% clear all
load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\WC_Germany\persistent_2010\')
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\Code\smoothing\software\')
addpath('F:\Code\general\')

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

rate_smooth_sig = round(Fsd*3);

for d = 1:100
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

    load used_data lf8 lf2 lf5 wcv_minus_spike
    load spike_time_jmm
    spkid = round(spkid/dsf);

    wcvf = filtfilt(b,a,wcv_minus_spike);
    wcvf = downsample(wcvf,dsf);
    wcvf = zscore(wcvf);
    lf8f = filtfilt(b,a,lf8);
    lf8f = downsample(lf8f,dsf);
    lf8f = zscore(lf8f);
    lf2 = lf2/sess_data(d).gains(2);
    lf5 = lf5/sess_data(d).gains(5);
    lf2_r = lf2-lf5;
    lf2_rf = filtfilt(b,a,lf2_r);
    lf2_rf = zscore(downsample(lf2_rf,dsf));
%     lf3f = filtfilt(b,a,lf3);
%     lf3f = downsample(lf3f,dsf);
%     lf3f = zscore(lf3f);
    wcvf2 = filtfilt(b2,a2,wcv_minus_spike);
    wcvf2 = downsample(wcvf2,dsf);
    wcvf2 = zscore(wcvf2);
    lf8f2 = filtfilt(b2,a2,lf8);
    lf8f2 = downsample(lf8f2,dsf);
    lf8f2 = zscore(lf8f2);
    lf2_rf2 = filtfilt(b2,a2,lf2_r);
    lf2_rf2 = zscore(downsample(lf2_rf2,dsf));
%     lf3f2 = filtfilt(b2,a2,lf3);
%     lf3f2 = downsample(lf3f2,dsf);
%     lf3f2 = zscore(lf3f2);
    rate_t = 1:length(lf8f);
    
    sm_rate = jmm_smooth_1d_cor(hist(spkid,rate_t),rate_smooth_sig);
    rate_t = rate_t/Fsd;
    
    ac_time = (1:length(wcvf))/Fsd;
    
    total_dur = ac_time(end);
    numWins = floor((total_dur-windowSize)/windowSlide);
    t_axis = (0:numWins-1)*windowSlide+windowSize/2;
    data_dist = zeros(numWins,length(amp_range));
    data_dist8 = zeros(numWins,length(amp_range));
    data_dist2r = zeros(numWins,length(amp_range));
    for w = 1:numWins
        begT = (w-1)*windowSlide;
        endT = begT + windowSize;
        begInd = round(Fsd*begT+1);
        endInd = begInd + round(Fsd*windowSize);
        data_seg = wcvf(begInd:endInd);
        data_dist(w,:) = gpkde(data_seg(:),0.06,[-3 4 500]);
        data_seg = lf8f(begInd:endInd);
        data_dist8(w,:) = gpkde(data_seg(:),0.06,[-3 5 500]);
        data_seg = lf2_rf(begInd:endInd);
        data_dist2r(w,:) = gpkde(data_seg(:),0.06,[-3 4 500]);
    end


    [Pw,t,f] = mtspecgramc(wcvf2,movingwin,params);
    [P8,t,f] = mtspecgramc(lf8f2,movingwin,params);
    [P2r,t,f] = mtspecgramc(lf2_rf2,movingwin,params);
    
%     f1 = figure;
%     set(f1,'papersize',[4.0 2.0]);
%     subplot(2,2,1)
%     pcolor(t_axis,amp_range,log(data_dist'));shading flat
%     caxis([-4 -0.5])
%     subplot(2,2,2)
%     plot(rate_t,sm_rate);
%     xlim([t_axis(1) t_axis(end)])
%     subplot(2,2,3)
%     pcolor(t,f,log(Pw'));shading flat
%     caxis([-8 0])
%     set(gca,'yscale','log')
%     subplot(2,2,4)
%     pcolor(t,f,log(P8'));shading flat
%     caxis([-8 0])
%     set(gca,'yscale','log')
%     set(gca,'yscale','log')
%     tname = ['G:\WC_Germany\overall_EC\nearulamp_spec_comb_rate\' s_name];
%     print('-dpng',tname);
%     close

    
    f1 = figure;
    set(f1,'papersize',[4.0 2.0]);
    subplot(3,2,1)
    pcolor(t_axis,amp_range,log(data_dist'));shading flat
    xlabel('Time (s)','Fontsize',12)
    ylabel('Amplitude (z)','fontsize',12)
    title('MP','fontsize',12)
    caxis([-4 -0.5])
    subplot(3,2,3)
    pcolor(t_axis,amp_range8,log(data_dist8'));shading flat
    xlabel('Time (s)','Fontsize',12)
    ylabel('Amplitude (z)','fontsize',12)
    title('Cortical LFP','fontsize',12)
    caxis([-4 -0.5])
    subplot(3,2,5)
    pcolor(t_axis,amp_range3,log(data_dist2r'));shading flat
    caxis([-4 -0.5])
    xlabel('Time (s)','Fontsize',12)
    ylabel('Amplitude (z)','fontsize',12)
    title('Hipp LFP','fontsize',12)
    subplot(3,2,2)
    pcolor(t,f,log(Pw'));shading flat
    caxis([-8 0])
    set(gca,'yscale','log')
    xlabel('Time (s)','Fontsize',12)
    ylabel('Frequency (Hz)','fontsize',12)
    title('MP','fontsize',12)
    subplot(3,2,4)
    pcolor(t,f,log(P8'));shading flat
    caxis([-8 0])
    set(gca,'yscale','log')
    xlabel('Time (s)','Fontsize',12)
    ylabel('Frequency (Hz)','fontsize',12)
    title('Cortical LFP','fontsize',12)
    subplot(3,2,6)
    pcolor(t,f,log(P2r'));shading flat
    caxis([-8 0])
    set(gca,'yscale','log')
    xlabel('Time (s)','Fontsize',12)
    ylabel('Frequency (Hz)','fontsize',12)
    title('Hipp LFP','fontsize',12)
    tname = ['F:\WC_Germany\overall_EC\nearulamp_spec_comb\' s_name];
    print('-dpng',tname);
    close

end