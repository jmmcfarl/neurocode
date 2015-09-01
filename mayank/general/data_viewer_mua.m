%% data viewer

clear all
close all

%% for hc_cort
load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cor_dir.mat')
load C:\WC_Germany\overall_calcs\HC_Cort\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\HC_Cort\UDS_synch_state_dur\UDS_synch_state_dur_data
% 
%% for EC
% load('C:\WC_Germany\overall_calcs\overall_dir.mat')
% load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
% load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data
% load C:\WC_Germany\overall_calcs\overall_info_file_data
% 
%%

d =1

cd(over_dir{d})
pwd

winSize = 10;

Fs = 2016;
dsf = 1;
niqf = Fs/2;
Fsd = Fs/dsf;

[b2,a2] = butter(2,150/niqf,'high');
[b,a] = butter(2,[0.01/niqf 200/niqf]);


load used_data lf8 lf5 lf7 wcv_minus_spike
if ~exist('used_data_lf6.mat')
    disp('creating lf6')
    get_lf6;
end
load used_data_lf6 lf6


down_trans8{d} = round(down_trans8{d}*8/dsf);
up_trans8{d} = round(up_trans8{d}*8/dsf);

lf8_lf = filter(b,a,lf7);
lf8_ld = downsample(lf8_lf,dsf);
lf8_ld = zscore(lf8_ld);

wcv_lf = filter(b,a,wcv_minus_spike);
wcv_ld = downsample(wcv_lf,dsf);
wcv_ld = zscore(wcv_ld);

lf8_f = filter(b2,a2,lf8);
lf8_d = downsample(lf8_f,dsf);
lf8_d = zscore(lf8_d);

lf7_f = filter(b2,a2,lf7);
lf7_d = downsample(lf7_f,dsf);
lf7_d = zscore(lf7_d);

% [mua_times7] = mua_extract(lf7,2);
% mua_times7 = round(mua_times7/dsf);
% mua_isis = diff(mua_times7)/Fsd;
% mua_bin = zeros(size(lf7));
% mua_bin(mua_times7) = 1;
% mua_locrate = jmm_smooth_1d(mua_bin,200)*Fsd;
% mua_up_times = zeros(size(mua_locrate));
% mua_up_times(mua_locrate > 20) = 1;
% mua_up_times = logical(mua_up_times);
% 
% mua_spikes_isup = mua_up_times(mua_times7);
% up_spikes = mua_times7(mua_spikes_isup);
% down_spikes = mua_times7(~mua_spikes_isup);
% 
% down_spikes_ld = [];
% for s = 1:(length(down_trans8{d})-1)
%     cur_spikes = mua_times7(mua_times7 > down_trans8{d}(s) & mua_times7 < up_trans8{d}(s+1));
%     down_spikes_ld = [down_spikes_ld;cur_spikes];
% end

% up_isis = diff(up_spikes)/Fsd;
% down_isis = diff(down_spikes)/Fsd;

lf6_f = filter(b2,a2,lf6);
lf6_d = downsample(lf6_f,dsf);
lf6_d = zscore(lf6_d);

lf5_f = filter(b2,a2,lf5);
lf5_d = downsample(lf5_f,dsf);
lf5_d = zscore(lf5_d);

% hist_range = linspace(0,1,100);
% hist(down_isis,hist_range)
% xlim([0 0.9])
% xlabel('ISI (s)')
% ylabel('Count')
% title('Down state ISI distribution')

t = (1:length(lf8_d))/Fsd;

numWins = floor(max(t)/winSize);
for i = 1:numWins
    seg_beg = (i-1)*winSize*Fsd+1;
    seg_end = i*winSize*Fsd;
    plot(t(seg_beg:seg_end),lf8_d(seg_beg:seg_end)/8+3,'k','linewidth',1)
    hold on
% cur_mua_times = mua_times7(find(mua_times7 > seg_beg & mua_times7 <= seg_end));
% cur_mua_isup = mua_up_times(cur_mua_times);
% cur_up_mua_times = cur_mua_times(cur_mua_isup);
% cur_down_mua_times = cur_mua_times(~cur_mua_isup);

% cur_down_mua_ld = down_spikes_ld(find(down_spikes_ld > seg_beg & down_spikes_ld <= seg_end));
% cur_isis = diff(cur_mua_times)/2016;
% loc_isis = zeros(1,length(cur_isis));
% for ci = 3:length(cur_isis)
%     loc_isis(ci) = mean(cur_isis(ci-2:ci));
% end
% isi_times = cur_mua_times(2:end);

    plot(t(seg_beg:seg_end),lf7_d(seg_beg:seg_end)/8+1,'m','linewidth',1)
        hold on
%         plot(t(cur_up_mua_times),lf7_d(cur_up_mua_times)/8+1,'g.')
%         plot(t(cur_down_mua_times),lf7_d(cur_down_mua_times)/8+1,'r.')
% plot(t(cur_down_mua_ld),lf7_d(cur_down_mua_ld),'r.')
%         plot(t(seg_beg:seg_end),mua_locrate(seg_beg:seg_end)/10,'g')
%     plot(t(seg_beg:seg_end),lf6_d(seg_beg:seg_end)/8-1,'c','linewidth',1)
%     plot(t(seg_beg:seg_end),lf5_d(seg_beg:seg_end)/8-3,'g','linewidth',1)
    plot(t(seg_beg:seg_end),lf8_ld(seg_beg:seg_end)+1,'linewidth',1)
    plot(t(seg_beg:seg_end),wcv_ld(seg_beg:seg_end)-2,'r','linewidth',1)
    ylim([-4 4])
    pause
    clf
    
end


