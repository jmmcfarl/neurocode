clear all
close all
%%
load G:\WC_Germany\wc_database.mat
addpath('G:\Code\WC_anal\general\')
addpath('G:\WC_Germany\Overall_EC\')
addpath('G:\Code\Chronux\spectral_analysis\continuous\')

drive_letter = 'G';

%%
dsf = 32;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.075/niqf 10/niqf]);

params.Fs = Fsd;
params.err = [2 0.05];
params.tapers = [3 5];
W = 0.05;
params.fpass = [0 5];
winlength = 25;
winslide = 5;
movingwin = [winlength winslide];

cax = [-3 0];

for d = 95:length(wc_db)
    
    cdir = wc_db(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(wc_db(d).region,'_',wc_db(d).celltype,'_',wc_db(d).date);
    load used_data lf8 lf3 wcv_minus_spike

    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = zscore(downsample(wcv_d,dsf));
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = zscore(downsample(lf8_d,dsf));
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = zscore(downsample(lf3_d,dsf));
    
    t_axis = (1:length(wcv_d))/Fsd;

    [Cw3,phiw3,Sw3,Sww,S33,t,f] = cohgramc(wcv_d,lf3_d,movingwin,params);
    [Cw8,phiw8,Sw8,Sww,S88,t,f] = cohgramc(wcv_d,lf8_d,movingwin,params);
    [C83,phi83,S83,S88,S33,t,f] = cohgramc(lf8_d,lf3_d,movingwin,params);
 
    S38 = conj(S83);
    S8w = conj(Sw8);
    S3w = conj(Sw3);
    
    num = abs(S3w.*S88 - S38.*S8w).^2;
    denom = (S33.*S88 - abs(S38).^2).*(Sww.*S88 - abs(Sw8).^2);
    partial_w3 = num./denom; 
 
    numer = abs(S8w.*S33 - S83.*S3w).^2;
    denomer = (S88.*S33 - abs(S83).^2).*(Sww.*S33 - abs(Sw3).^2);
    partial_w8 = numer./denomer; 

    numer = abs(S38.*Sww - S3w.*Sw8).^2;
    denomer = (S33.*Sww - abs(S3w).^2).*(S88.*Sww - abs(S8w).^2);
    partial_83 = numer./denomer;

    Fig = figure('visible','off');
set(Fig,'PaperUnits','centimeters');
set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(3,2,1)
    pcolor(t,f,log10(Sww'));shading flat;
    xl = xlim();
    set(gca,'yscale','log')
    caxis(cax)
    title('MP')
    subplot(3,2,3)
    pcolor(t,f,log10(S88'));shading flat;
    caxis(cax)
    set(gca,'yscale','log')
    title('LF8')
    subplot(3,2,5)
    pcolor(t,f,log10(S33'));shading flat;
    caxis(cax)
    set(gca,'yscale','log')    
    title('LF3')
    subplot(3,2,2)
    pcolor(t,f,partial_w8');shading flat;
    caxis([0 1])
    set(gca,'yscale','log')
    title('Partial w8')
subplot(3,2,4)    
    pcolor(t,f,partial_w3');shading flat;
    caxis([0 1])
    set(gca,'yscale','log')
    title('Partial w3')
    subplot(3,2,6)
    pcolor(t,f,partial_83');shading flat;
    caxis([0 1])
   set(gca,'yscale','log')
   title('Partial 83')
    tname = ['G:\WC_Germany\overall_EC\hc_cort_partial_theta_grams_log\' s_name];
    print('-dpng',tname);
    close

    Fig = figure('visible','off');
set(Fig,'PaperUnits','centimeters');
set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(3,2,1)
    pcolor(t,f,log10(Sww'));shading flat;
    xl = xlim();
    caxis(cax)
    title('MP')
    subplot(3,2,3)
    pcolor(t,f,log10(S88'));shading flat;
    caxis(cax)
    title('LF8')
    subplot(3,2,5)
    pcolor(t,f,log10(S33'));shading flat;
    caxis(cax)
    title('LF3')
    subplot(3,2,2)
    pcolor(t,f,partial_w8');shading flat;
    caxis([0 1])
    title('Partial w8')
subplot(3,2,4)    
    pcolor(t,f,partial_w3');shading flat;
    caxis([0 1])
    title('Partial w3')
    subplot(3,2,6)
    pcolor(t,f,partial_83');shading flat;
    caxis([0 1])
    title('Partial 83')
    tname = ['G:\WC_Germany\overall_EC\hc_cort_partial_theta_grams\' s_name];
    print('-dpng',tname);
    close

    
end

%%
