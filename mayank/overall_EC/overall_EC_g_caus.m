clear all
close all
%%
load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\Code\WC_anal\general\')
addpath('F:\WC_Germany\Overall_EC\')
addpath('F:\Code\Chronux\spectral_analysis\continuous\')
addpath('F:\Code\bsmart\')

drive_letter = 'F';

%%
dsf = 32;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.06/niqf 10/niqf]);

params.Fs = Fsd;
params.err = [2 0.05];
params.tapers = [10 19];
W = 0.05;
params.fpass = [0 30];
winlength = 15;
winslide = 15;
movingwin = [winlength winslide];

f_i = linspace(W,10,250);

for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf8 lf3 wcv_minus_spike

    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = downsample(wcv_d,dsf);
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_d,dsf);
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = downsample(lf3_d,dsf);

    data = [wcv_d lf3_d];
    %compute AIC as a function of model order
    aic_vec = aic_test(data,15,20);
    %set model order to first value where the relative change is less than 2%
    model_order = 1+find(diff(aic_vec)./aic_vec(2:end) < 0.02,1,'first');
    if isempty(model_order)
        disp('AIC did not converge, setting model order to 10')
        model_order = 10;
    end
    
    freq_vec = 0:0.06:5; %frequencies to evaluate g-causality
    disp(sprintf('estimating moving window g-caus with model order %d',model_order))
    [Fxy,Fyx] = mov_bi_ga(data,1,round(length(wcv_d)/5),round(Fsd*winlength),model_order,Fsd,freq_vec);
    Fxy = squeeze(Fxy);
    Fyx = squeeze(Fyx);
    t_ar = linspace(0,round(length(wcv_d)/5)/Fsd,size(Fxy,2));
    subplot(2,1,1)
    pcolor(t_ar,freq_vec,Fxy);shading flat, caxis([0 0.4])
    title('MP -> LF3')
    subplot(2,1,2)
    pcolor(t_ar,freq_vec,Fyx);shading flat, caxis([0 0.4])
    title('LF3 -> MP')

    
end