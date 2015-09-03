%%
clear all
addpath('~/James_scripts/mayank/parietal_cortical_2010/');
addpath('~/James_scripts/mayank/persistent_9_27_2010/');
addpath('~/James_scripts/mayank/hsmm_state_detection/');
addpath('~/James_scripts/mayank/persistent_downs/');
addpath('~/James_scripts/mayank/final_pdown_analysis/');


load ~/Analysis/Mayank/final_pdown_analysis/compiled_data.mat
% load('~/Analysis/Mayank/final_pdown_analysis/compiled_corticalMP_data.mat');


min_rec_dur = 500; %minimum total duration of recording (in sec)
data_ids = [data(:).id];
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);
used_dirs(ismember(data_ids(used_dirs),no_cell) | ismember(data_ids(used_dirs),unclear_uds)) = [];

data = data(used_dirs);

load ~/Analysis/Mayank/final_pdown_analysis/fin_pdown_core_analysis_fin.mat
% load('~/Analysis/Mayank/final_pdown_analysis/fin_pdown_core_analysis_corticalMP_fin.mat');
if length(core_data) ~= length(data)
    error('Data mismatch');
end

%% parameters
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
lcf = 0.05; hcf = 10; %bandpass filter for LF amp

%% get basic data 
d = 6; %example id

cd(data(d).new_dir)
pwd

load ./used_data lf7 wcv_minus_spike

%get LF amps
[lfp_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);

%get rid of data after the last good sample
end_time = min(data(d).ep,data(d).dp);
ep = find(t_axis >= end_time,1);
if ~isempty(ep)
    lfp_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; 
else
    ep = length(t_axis);
end

%% get state seq data
load ./pa_hsmm_state_seq_combined_fin_nd.mat
load ./pa_hsmm_state_seq7_combined_fin_nd.mat
hsmm_ctx = hsmm7;
lfp_state_seq = hsmm_bbstate_seq7;
mp_state_seq = hsmm_bbstate_seq;

[new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,mp_state_seq);
dur_uds = sum(diff(new_seg_inds,[],2))/Fsd; %total duration of UDS segs
    %get state transition times for lfp and MP
    [up_trans_inds_mp,down_trans_inds_mp] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds_lfp,down_trans_inds_lfp] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);

    %eliminate state transitions after the data end point
    bad_mp_states = find(up_trans_inds_mp > ep | down_trans_inds_mp > ep);
    up_trans_inds_mp(bad_mp_states) = []; down_trans_inds_mp(bad_mp_states) = [];
    bad_lfp_states = find(up_trans_inds_lfp > ep | down_trans_inds_lfp > ep);
    up_trans_inds_lfp(bad_lfp_states) = []; down_trans_inds_lfp(bad_lfp_states) = [];

%% get cortical period data
load ./allEC_ctx_period_data_hsmm.mat
lfp_period_vec = nan(size(wcv_lf));
for i = 1:size(new_seg_inds,1)
    cur_inds = new_seg_inds(i,1):new_seg_inds(i,2);
    cur_inds_used = find(cur_inds <= ep);
    lfp_period_vec(cur_inds(cur_inds_used)) = lf8_period_f{i}(cur_inds_used);
end

    %% get MP and LFP UDS state vectors
    mp_state_number = nan(length(wcv_lf),1);
    mp_state_vec = zeros(length(wcv_lf),1);
    lfp_state_number = nan(length(wcv_lf),1);
    lfp_state_vec = zeros(length(wcv_lf),1);
    for i = 1:length(up_trans_inds_mp)-1
        mp_state_vec(up_trans_inds_mp(i):down_trans_inds_mp(i)) = 1;
        mp_state_number(up_trans_inds_mp(i):down_trans_inds_mp(i)) = 2*(i-1)+1;
        mp_state_number(down_trans_inds_mp(i):up_trans_inds_mp(i+1)) = 2*(i-1) + 2;
    end
    for i = 1:length(up_trans_inds_lfp)-1
        lfp_state_vec(up_trans_inds_lfp(i):down_trans_inds_lfp(i)) = 1;
        lfp_state_number(up_trans_inds_lfp(i):down_trans_inds_lfp(i)) = 2*(i-1)+1;
        lfp_state_number(down_trans_inds_lfp(i):up_trans_inds_lfp(i+1)) = 2*(i-1) + 2;
    end
    mp_state_number(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;
    mp_state_vec(isnan(lfp_period_vec)) = nan;
    lfp_state_vec(isnan(lfp_period_vec)) = nan;

%%
% load ./aligned_heka.mat
% ac_time_d = downsample(ac_time,dsf);
% dc_fs = 1/median(diff(dc_time));
% dc_dsf = 20; %ddown-sample fac for DC data
% dc_fsd = dc_fs/dc_dsf;
% dc_data = decimate(dc_data,dc_dsf);
% dc_time = downsample(dc_time,dc_dsf);
% 
% %correct for junction potential
% jp = -7; %in mV
% dc_data = dc_data + jp;

%%
%%
close all
window_size = 20;
window_size_n = round(Fsd*window_size);
% window_size_dc = round(dc_fsd*window_size);
n_windows = ceil(max(t_axis)/window_size);
for i = 1:n_windows
    cur_set = (i-1)*window_size_n + (1:window_size_n);
    cur_set(cur_set > length(t_axis)) = [];
%     cur_set_dc = find(dc_time >= ac_time_d(cur_set(1)) & dc_time < ac_time_d(cur_set(end)));
    
    subplot(2,1,1); hold on
%     plot(dc_time(cur_set_dc),dc_data(cur_set_dc));hold on
plot(t_axis(cur_set),wcv_lf(cur_set));
plot(t_axis(cur_set),(mp_state_vec(cur_set)),'k')
    set(gca,'fontsize',10,'fontname','arial');
    box off
    
    subplot(2,1,2)
    hold on
    plot(t_axis(cur_set),lfp_lf(cur_set),'r','linewidth',1)
    plot(t_axis(cur_set),lfp_state_vec(cur_set)-1,'k')
    set(gca,'fontsize',10,'fontname','arial');
    box off
    
%     subplot(3,1,3)
%     hold on
%     plot(ac_time_d(cur_set),ctx_mua_rate(cur_set)-5,'r')
%     plot(ac_time_d(cur_set),hpc_mua_rate(cur_set)-8,'k')
%     set(gca,'fontsize',10,'fontname','arial');
%     box off
%     xlabel('Time (s)','fontsize',14)
    ylabel('Amplitude (V)','fontsize',14)
 
    pause
    clf
end

