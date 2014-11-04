addpath(genpath('~/James_scripts_old/Nlx2Mat_relDec11/'))
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))
%%
clear all

% cd /Users/james/Analysis/Mayank/sleep/2012-7-5_Patch/2012-7-5_17-11-29 %actually some decent LFP UDS, MP pretty much unclear though. Maybe one place where a pdown is conceivable.
% cd /Users/james/Analysis/Mayank/sleep/2012-7-9_Sleep_Patch/2012-7-9_12-59-28  %very short, nothing really clear
% cd /Users/james/Analysis/Mayank/sleep/2012-7-11_sleep_WC/2012-7-11_13-42-46 %not much clear LFP UDS
% cd /Users/james/Analysis/Mayank/sleep/2012-7-8_2/2012-7-8_13-9-45 %very long rec, not clear MP uds though
% cd /Users/james/Analysis/Mayank/sleep/2012-7-9_Sleep_Patch/2012-7-9_12-59-28 %too short
% cd /Users/james/Analysis/Mayank/sleep/2012_07_03/2012-7-3_15-25-44 %and others from this date, not really usable MP
% cd /Users/james/Analysis/Mayank/sleep/Sleep_07-06/2012-7-6_16-12-14 %not really usable MP

% cd /Users/james/Analysis/Mayank/sleep/2012-7-10_Sleep_WC/2012-7-10_13-49-40 %MP has wierd activity with apparently long down states, and very short unstable ups. could be perceived as pers downs...
cd /Users/james/Analysis/Mayank/sleep/2012-7-13#1/2012-7-13_13-31-49 %nice long rec with best examples ive seen so far
% cd /Users/james/Analysis/Mayank/sleep/2012-7-13#2/2012-7-13_16-54-15 %long recording, some descent MP UDS at times, but no clear pers downs

%% LOAD RAW DATA

FieldSelectionFlags = [1 0 1 0 1];
ExtractMode = 1;
HeaderExtractionFlag = 1;

i = 1
fprintf('Loading CSC%d\n',i)
    Filename = sprintf('CSC%d.Ncs',i);
%     [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
%         Nlx2MatCSC(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
        Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    mp_t = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    Fs(i) = median(CSC_SampleFrequencies);
    mp = CSC_Samples(:)*str2num(Header{15}(13:end));
    clear CSC_Samples
mp_d = decimate(mp,16);
mp_t = downsample(mp_t,16);
mp_Fs = Fs/16;

for i = 2:8
    fprintf('Loading CSC%d\n',i)
    Filename = sprintf('CSC%d.Ncs',i);
    [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
        Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    csc_time{i} = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    Fs(i) = median(CSC_SampleFrequencies);
    csc{i} = CSC_Samples(:)*str2num(Header{15}(13:end));
    bad_samps = find(isnan(csc_time{i}));
    csc_time{i}(bad_samps) = [];
    csc{i}(bad_samps) = [];
    clear CSC_Samples
end
Fs = unique(Fs(2));


%% DOWN-SAMPLE

dsf = 1;
Fsd = Fs/dsf;
[b,a] = butter(2,[0.2 10]/(Fs/2));
[b2,a2] = butter(2,[30 100]/(Fs/2));
hf_smooth = round(0.025*Fs);
for i = 2:8
    csc_lf{i} = zscore(filtfilt(b,a,csc{i}));
    csc_hf{i} = zscore(filtfilt(b2,a2,csc{i}));
    csc_hf{i} = zscore(sqrt(jmm_smooth_1d_cor(csc_hf{i}.^2,hf_smooth)));
end

[b,a] = butter(2,[0.1 40]/(mp_Fs/2));
mp_lf = zscore(filtfilt(b,a,-mp_d));

%%
true_heka_fs = 1.99984e+004;
heka_dsf = 10;
heka_fsd = true_heka_fs/heka_dsf;

load ./aligned_heka.mat
dc_time = (1:length(heka_data))/true_heka_fs+dc_offset+mp_t(1)/1e6;
dc_time_d = downsample(dc_time,heka_dsf);
dc_data_d = zscore(decimate(heka_data,heka_dsf));
% use_points = find(dc_time_d > sp & dc_time_d < ep);
% std_scale = std(dc_data_d(use_points));
% dc_data_d = (dc_data_d - mean(dc_data_d(use_points)))/std(dc_data_d(use_points));

[b,a] = butter(1,[0.2 40]/true_heka_fs);
dc_data_lf = downsample(filtfilt(b,a,heka_data),heka_dsf);
% dc_data_lf = (dc_data_lf - mean(dc_data_lf(use_points)))/std(dc_data_lf(use_points));

%%
close all

figure
hold on
% [ff,h1,h2] = plotyy(dc_time_d,dc_data_d,csc_time{6}/1e6,csc_f{6})
% [ff,h1,h2] = plotyy(csc_time{6}/1e6,csc_f{6},dc_time_d,dc_data_d)
% plot(dc_time_d,zscore(dc_data_d)*1.5+3,'r','linewidth',1)
plot(mp_t/1e6,mp_lf+3,'r')
plot(csc_time{6}/1e6,csc_lf{6},'k','linewidth',1)
plot(csc_time{6}/1e6,zscore(csc_hf{6})-4,'b','linewidth',1)
% plot(csc_time{6}/1e6,csc_lf{6},'b','linewidth',2)
% plot(csc_time{6}(ctx_ds_locs)/1e6,csc_lf{6}(ctx_ds_locs),'ro','linewidth',2)
% xlim(examp_1)
xlabel('Time (s)','fontsize',18)
ylabel('Amplitude (z)','fontsize',18)
set(gca,'fontsize',16,'fontname','arial')

%%
win_size = 15;
n_bins = ceil(range(mp_t/1e6)/win_size);
for ii = 1:n_bins
    ep = (ii-1)*win_size + mp_t(1)/1e6;
    xlim([ep ep + win_size]);
    pause
end
%% LOAD MUA
amp_threshold = 25;
max_overlap = 0.5;
[mua,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
clear binned* sm_rate
dsf = 1;
down_ts = downsample(csc_time{6},dsf);
for ch = 1:8;
    cur_mua = mua{ch};
    cur_mua(cur_mua > down_ts(end) | cur_mua < down_ts(1)) = [];
    binned_spks(:,ch) = hist(cur_mua,down_ts)';
    sm_rate(:,ch) = zscore(jmm_smooth_1d_cor(binned_spks(:,ch),round(Fsd*0.05)));
end
binned_spks(end,:) = mean(binned_spks);

lfp_mua = sum(binned_spks(:,6:8),2);
lfp_rate = zscore(jmm_smooth_1d_cor(lfp_mua,round(Fsd*0.05)));

hpc_mua = sum(binned_spks(:,2),2);
hpc_rate = zscore(jmm_smooth_1d_cor(hpc_mua,round(Fsd*0.05)));

plot(csc_time{6}/1e6,lfp_rate-5,'g')
%%
maxlag = round(Fs*1);
maxlag_mp = round(heka_fsd*1);

ds_threshold = [1 1.5 2 2.5 3];
% ds_threshold = 1.5;
for dd = 1:length(ds_threshold)
    [ctx_ds_amps,ctx_ds_locs] = findpeaks(-csc_lf{6},'minpeakheight',ds_threshold(dd),'minpeakdistance',min_sep);
    ctx_ds_amps(~ismember(ctx_ds_locs,uds_inds)) = [];
    ctx_ds_locs(~ismember(ctx_ds_locs,uds_inds)) = [];
    
    [~,mec_ds_locs] = findpeaks(-dc_data_lf,'minpeakheight',ds_threshold(dd),'minpeakdistance',min_sep);
    mec_ds_locs(~ismember(mec_ds_locs,uds_inds_mec)) = [];

    ctx_ds_locs_interp = round(interp1(dc_time_d,1:length(dc_time_d),csc_time{6}(ctx_ds_locs)/1e6));
    bad_ds = find(isnan(ctx_ds_locs_interp));
    bad_ds = [bad_ds find(csc_time{6}(ctx_ds_locs)/1e6 < sp | csc_time{6}(ctx_ds_locs)/1e6 > ep)];
    bad_ds = [bad_ds find(ctx_ds_locs < maxlag | ctx_ds_locs > length(csc_lf{6})-maxlag)];
    bad_ds = [bad_ds find(ctx_ds_locs_interp < maxlag_mp | ctx_ds_locs_interp > length(dc_data_lf)-maxlag_mp)];
    ctx_ds_locs_interp(bad_ds) = [];
    ctx_ds_locs(bad_ds) = [];
    ctx_ds_amps(bad_ds) = [];
    
    inter_ds_intervals{dd} = diff(ctx_ds_locs)/Fs;
    mec_inter_ds_intervals{dd} = diff(mec_ds_locs)/heka_fsd;
    
    ds_trig_mp = zeros(length(ctx_ds_locs),2*maxlag_mp+1);
    ds_trig_lfp = zeros(length(ctx_ds_locs),2*maxlag+1);
%     ds_trig_mua = zeros(length(ctx_ds_locs),2*maxlag+1);
    for i = 1:length(ctx_ds_locs)
        uset1 = (ctx_ds_locs(i)-maxlag):(ctx_ds_locs(i)+maxlag);
        ds_trig_lfp(i,:) = csc_lf{6}(uset1);
%         ds_trig_mua(i,:) = lfp_rate(uset1);
        uset2 = (ctx_ds_locs_interp(i)-maxlag_mp):(ctx_ds_locs_interp(i)+maxlag_mp);
        ds_trig_mp(i,:) = dc_data_lf(uset2);
    end
    trig_lfp_avg(dd,:) = mean(ds_trig_lfp);
    trig_mp_avg(dd,:) = mean(ds_trig_mp);
%     trig_mua_avg(dd,:) = mean(ds_trig_mua);
end

lags = (-maxlag:maxlag)/Fsd;
heka_lags = (-maxlag_mp:maxlag_mp)/heka_fsd;

%%
cmap = jet(4);
xi = linspace(0.1,20,100);
for i = 1:4
    y = ksdensity(inter_ds_intervals{i},xi,'support','positive');
    subplot(2,1,1)
    plot(xi,y,'color',cmap(i,:))
    hold on
     set(gca,'yscale','log')
ylim([0.005 1])
xlim([0 20])
     y = ksdensity(mec_inter_ds_intervals{i},xi,'support','positive');
    subplot(2,1,2)
    plot(xi,y,'color',cmap(i,:))
    hold on
     set(gca,'yscale','log')
ylim([0.005 1])
xlim([0 20])
end    
   
%%
xi = linspace(-4,4,100);
mp_trig_dist = zeros(length(heka_lags),100);
for i = 1:length(heka_lags)
   mp_trig_dist(i,:) = ksdensity(ds_trig_mp(:,i),xi); 
end
    
%%
figure
subplot(3,1,1)
plot(lags,trig_lfp_avg)
xlabel('Time (s)','fontsize',18)
ylabel('Amplitude (z)','fontsize',18)
axis tight
xlim([-0.75 0.75])
title('Cortical LFP','fontsize',20)
set(gca,'fontsize',16)
grid on
subplot(3,1,2)
plot(heka_lags,trig_mp_avg)
xlabel('Time (s)','fontsize',18)
ylabel('Amplitude (z)','fontsize',18)
axis tight
xlim([-0.75 0.75])
set(gca,'fontsize',16)
grid on
title('MEC MP','fontsize',20)

subplot(3,1,3)
plot(lags,trig_mua_avg)
xlabel('Time (s)','fontsize',18)
ylabel('Amplitude (z)','fontsize',18)
axis tight
xlim([-0.75 0.75])
set(gca,'fontsize',16)
grid on
title('Cortical MUA','fontsize',20)


%%
% temp = pwd;
% sl = find(temp == '\',1,'last');
% dname = temp(sl+1:end);
% save(dname,'binned_spks','dsf', 'Fsd', 't', 'cscd');