% clear all
load G:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010.mat
addpath('G:\WC_Germany\persistent_2010\')
addpath('G:\Code\smoothing\software\')

dsf = 16;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.05/niqf 10/niqf]);
[b2,a2] = butter(2,[0.05/niqf 40/niqf]);

params.Fs = Fsd;
params.err = 0;
params.tapers = [2 3];
params.fpass = [0 30];
windowSize = 20;
windowSlide = 2;
movingwin = [windowSize windowSlide];

amp_range = linspace(-3,4,500);
amp_range8 = linspace(-3,5,500);
amp_range3 = linspace(-3,4,500);

for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

    load used_data lf8 lf3 wcv_minus_spike

    wcvf = filtfilt(b,a,wcv_minus_spike);
    wcvf = downsample(wcvf,dsf);
    wcvf = zscore(wcvf);
    lf8f = filtfilt(b,a,lf8);
    lf8f = downsample(lf8f,dsf);
    lf8f = zscore(lf8f);
    lf3f = filtfilt(b,a,lf3);
    lf3f = downsample(lf3f,dsf);
    lf3f = zscore(lf3f);
    wcvf2 = filtfilt(b2,a2,wcv_minus_spike);
    wcvf2 = downsample(wcvf2,dsf);
    wcvf2 = zscore(wcvf2);
    lf8f2 = filtfilt(b2,a2,lf8);
    lf8f2 = downsample(lf8f2,dsf);
    lf8f2 = zscore(lf8f2);
    lf3f2 = filtfilt(b2,a2,lf3);
    lf3f2 = downsample(lf3f2,dsf);
    lf3f2 = zscore(lf3f2);

    ac_time = (1:length(wcvf))/Fsd;
    
    total_dur = ac_time(end);
    numWins = floor((total_dur-windowSize)/windowSlide);
    t_axis = (0:numWins-1)*windowSlide+windowSize/2;
    data_dist = zeros(numWins,length(amp_range));
    data_dist8 = zeros(numWins,length(amp_range));
    data_dist3 = zeros(numWins,length(amp_range));
    for w = 1:numWins
        begT = (w-1)*windowSlide;
        endT = begT + windowSize;
        begInd = round(Fsd*begT+1);
        endInd = begInd + round(Fsd*windowSize);
        data_seg = wcvf(begInd:endInd);
        data_dist(w,:) = gpkde(data_seg(:),0.06,[-3 4 500]);
        data_seg = lf8f(begInd:endInd);
        data_dist8(w,:) = gpkde(data_seg(:),0.06,[-3 5 500]);
        data_seg = lf3f(begInd:endInd);
        data_dist3(w,:) = gpkde(data_seg(:),0.06,[-3 4 500]);
    end


    [Pw,t,f] = mtspecgramc(wcvf2,movingwin,params);
    [P8,t,f] = mtspecgramc(lf8f2,movingwin,params);
    [P3,t,f] = mtspecgramc(lf3f2,movingwin,params);

    f1 = figure;
    set(f1,'papersize',[4.0 2.0]);
    subplot(3,2,1)
    pcolor(t_axis,amp_range,log(data_dist'));shading flat
    caxis([-4 -0.5])
    subplot(3,2,3)
    pcolor(t_axis,amp_range8,log(data_dist8'));shading flat
    caxis([-4 -0.5])
    subplot(3,2,5)
    pcolor(t_axis,amp_range3,log(data_dist3'));shading flat
    caxis([-4 -0.5])
    subplot(3,2,2)
    pcolor(t,f,log(Pw'));shading flat
    caxis([-8 0])
    set(gca,'yscale','log')
    subplot(3,2,4)
    pcolor(t,f,log(P8'));shading flat
    caxis([-8 0])
    set(gca,'yscale','log')
    subplot(3,2,6)
    pcolor(t,f,log(P3'));shading flat
    caxis([-8 0])
    set(gca,'yscale','log')
    tname = ['G:\WC_Germany\parietal_cortical_2010\nearulamp_spec_comb\' s_name];
    print('-dpng',tname);
    close

end