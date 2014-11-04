clear all
close all

addpath('C:/Code/Nlx2Mat_new/')
addpath('C:/Code/general_functions/')

%7-10-13
%7-13-13

% cd C:\WC_Germany\sleep\2012-7-5_Patch\2012-7-5_17-11-29 %not too clear
% cd C:\WC_Germany\sleep\2012-7-9_Sleep_Patch\2012-7-9_12-59-28 %one short example. questionable i guess
cd C:\WC_Germany\sleep\2012-7-10_Sleep_WC\2012-7-10_13-49-40 %not really so clear, but might have pers downs
% cd C:\WC_Germany\sleep\2012-7-11_sleep_WC\2012-7-11_13-14-7 %no real clear examps
% cd C:\WC_Germany\sleep\2012-7-11_sleep_WC\2012-7-11_13-42-46 %couple OK examples with descent LFP and MUA
% cd C:\WC_Germany\sleep\2012-7-13#1\2012-7-13_12-47-19 %no clear examps
% cd C:\WC_Germany\sleep\2012-7-13#1\2012-7-13_13-31-49 %a number of good examples with good LFP, and OK MUA, might have some pers downs
% cd C:\WC_Germany\sleep\2012-7-13#2\2012-7-13_16-9-10 %no good examples
% cd C:\WC_Germany\sleep\2012-7-13#2\2012-7-13_16-54-15 %no clear examples, lack of UDS

%% LOAD RAW DATA

load ./aligned_heka


FieldSelectionFlags = [1 0 1 0 1];
ExtractMode = 1;
HeaderExtractionFlag = 1;

i = 1
fprintf('Loading CSC%d\n',i)
    Filename = sprintf('CSC%d.Ncs',i);
    [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
        Nlx2MatCSC(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    mp_t = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    Fs(i) = median(CSC_SampleFrequencies);
    mp = CSC_Samples(:)*str2num(Header{15}(13:end));
    clear CSC_Samples
mp_d = decimate(mp,32);
mp_t = downsample(mp_t,32);

for i = 2:8
    fprintf('Loading CSC%d\n',i)
    Filename = sprintf('CSC%d.Ncs',i);
    [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
        Nlx2MatCSC(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    csc_time{i} = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    Fs(i) = median(CSC_SampleFrequencies);
    csc{i} = CSC_Samples(:)*str2num(Header{15}(13:end));
%     bad_samps = find(isnan(csc_time{i}));
bad_samps = find(csc_time{i} < mp_t(1) | csc_time{i} > mp_t(end));
    csc_time{i}(bad_samps) = [];
    csc{i}(bad_samps) = [];
    clear CSC_Samples
csc{i} = decimate(csc{i},2);
csc_time{i} = downsample(csc_time{i},2);
Fs(i) = Fs(i)/2;
end
Fs = unique(Fs(2));


%% DOWN-SAMPLE
dsf = 1;
Fsd = Fs/dsf;
[b,a] = butter(2,[0.2 10]/(Fs/2));
for i = 2:8
    csc_f{i} = zscore(filtfilt(b,a,csc{i}));
end

%% LOAD MUA
amp_threshold = 25;
max_overlap = 0.5;
[mua,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
% load ./sync_times.mat
% synct_d = downsample(synct,8);
% counts = cellfun(@(x) length(x),mua_amps)/range(synct)*1e6
% figure
% plot(1:7,counts(1:7),'o-')

%% BIN MUA
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

ctx_mua = sum(binned_spks(:,6:8),2);
ctx_rate = zscore(jmm_smooth_1d_cor(ctx_mua,round(Fsd*0.05)));

hpc_mua = sum(binned_spks(:,4),2);
hpc_rate = zscore(jmm_smooth_1d_cor(hpc_mua,round(Fsd*0.05)));

figure
plot(mean(binned_spks),'o-');
%%
% ulist_hpc = [2:4];
% ulist_ctx = [5:8];
% master_spks_hpc = [];
% master_spks_ctx = [];
% for i = ulist_hpc
% master_spks_hpc = [master_spks_hpc mua{i}];
% end
% for i = ulist_ctx
% master_spks_ctx = [master_spks_ctx mua{i}];
% end

win = 20;
bt = 0;
ac_time_d = downsample(ac_time,2);
ac_time_d(length(csc_f{6})+1:end) = [];

dc_data = downsample(heka_data,10);
dc_data(dc_data > -40) = -40;


add_dsf = 4;

figure
plot(downsample(ac_time_d,add_dsf),downsample(csc_f{6},add_dsf))
hold on
plot(downsample(dc_time,10*add_dsf),downsample(dc_data,add_dsf)/3+23,'r')
plot(downsample(ac_time_d,add_dsf),downsample(ctx_rate-5,add_dsf),'k')
plot(downsample(ac_time_d,add_dsf),downsample(hpc_rate-8,add_dsf),'g')
legend('Ctx LFP','MEC MP','Ctx MUA','Hpc MUA');
xlabel('Time (s)')
ylabel('Amplitude')
for i = 1:100
    xlim([bt bt + win]);
    pause
    bt = bt + win;
%     clf
end

% raster(master_spks_hpc'/1e6,csc_time{2}([1 end])/1e6,'r',0);
% hold on
% raster(master_spks_ctx'/1e6,csc_time{2}([1 end])/1e6,'b',4);
% hold on
% plot(csc_time{6}/1e6,csc_f{6}*1.5+10,'b')
% plot(mp_t/1e6,-zscore(mp_d)*2+16,'k')
% plot(csc_time{6}/1e6,lfp_rate*2,'r')
% axis tight
% plot(csc_time{6}/1e6,hpc_rate*2,'g')


%%
% sp = 100*1e6;
% ep = 220*1e6;

% sp = 1175*1e6;
% ep = 1320*1e6;

sp = 550*1e6;
ep = 1000*1e6;

ctx_ch = 8;
sind_mp = find(mp_t > sp,1,'first');
eind_mp = find(mp_t > ep,1,'first');
use_mp = -zscore(mp_d(sind_mp:eind_mp));
[b,a] = butter(2,[0.1 4]/(2016/2));
use_mp_lp = zscore(filtfilt(b,a,use_mp));

sind_lfp = find(csc_time{ctx_ch} > sp,1,'first');
eind_lfp = find(csc_time{ctx_ch} > ep,1,'first');
use_lfp = csc{ctx_ch}(sind_lfp:eind_lfp);
use_mua = lfp_rate(sind_lfp:eind_lfp);
use_norm_lfp = csc_f{ctx_ch}(sind_lfp:eind_lfp);
use_lfp_t = csc_time{ctx_ch}(sind_lfp:eind_lfp);

% params.Fs = 2016;
% params.tapers = [4 7];
% win = 20;
% [Smp,f] = mtspectrumsegc(use_mp,win,params);
% [Slfp,f] = mtspectrumsegc(use_lfp,win,params);

maxlag = round(2016*1);
[ctx_ds_amps,ctx_ds_locs] = findpeaks(-use_norm_lfp,'minpeakheight',2);
% [ctx_ds_amps,ctx_ds_locs] = findpeaks(-use_mp_lp,'minpeakheight',2);
bad_ds = find(ctx_ds_locs < maxlag | ctx_ds_locs > length(use_lfp_t) - maxlag);
ctx_ds_locs(bad_ds) = [];
ctx_ds_locs_interp = round(interp1(mp_t(sind_mp:eind_mp),1:length(use_mp),use_lfp_t(ctx_ds_locs)));

ds_trig_mp = zeros(length(ctx_ds_locs),2*maxlag+1);
ds_trig_lfp = zeros(length(ctx_ds_locs),2*maxlag+1);
ds_trig_mua = zeros(length(ctx_ds_locs),2*maxlag+1);
for i = 1:length(ctx_ds_locs)
    uset1 = (ctx_ds_locs(i)-maxlag):(ctx_ds_locs(i)+maxlag);
    ds_trig_lfp(i,:) = use_norm_lfp(uset1);
    ds_trig_mua(i,:) = use_mua(uset1);
    uset2 = (ctx_ds_locs_interp(i)-maxlag):(ctx_ds_locs_interp(i)+maxlag);
    ds_trig_mp(i,:) = use_mp_lp(uset2);
end

lags = (-maxlag:maxlag)/Fsd;

%%
% temp = pwd;
% sl = find(temp == '\',1,'last');
% dname = temp(sl+1:end);
% save(dname,'binned_spks','dsf', 'Fsd', 't', 'cscd');