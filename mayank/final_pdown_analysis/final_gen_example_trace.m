%%
clear all
load C:/WC_Germany/final_pdown_analysis/compiled_data.mat
load C:/WC_Germany/final_pdown_analysis/mua_classification.mat

fig_dir = 'C:\WC_Germany\final_pdown_analysis\figures\';

min_rec_dur = 500; %in sec
data_ids = [data(:).id];
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);
used_dirs(ismember(data_ids(used_dirs),unclear_uds)) = [];

data = data(used_dirs);

peak_hpcmua_loc = peak_hpcmua_loc(used_dirs);
peak_hpcmua_rate = peak_hpcmua_rate(used_dirs);
usable_mua = usable_mua(used_dirs,:);

load C:/WC_Germany/final_pdown_analysis/fin_pdown_core_analysis.mat
if length(core_data) ~= length(data)
    error('Data mismatch');
end

min_hpcrate = 1;
usable_hpc_mua = ~isnan(peak_hpcmua_loc) & peak_hpcmua_rate >= min_hpcrate;

%% parameters
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;

rate_sm = round(Fsd*0.05);

%%
d = 34;

cd(data(d).dir)
pwd

load ./used_data lf7 lf3 lf5 wcv_minus_spike
if data(d).hpc_lfp == 3
    load ./used_data lf3 lf5
    if data(d).is_old_type
        lf3 = lf3 + lf5; %redefine LF3 wrt gnd
    end
    hpc_lfp = lf3;
else
    load ./used_data lf2
    hpc_lfp = lf2;
end

[lfp_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
hpc_lf = get_lf_features(hpc_lfp,raw_Fs,Fsd,[lcf hcf]);

load ./sync_times.mat
synct_d = downsample(synct,dsf);

load ./spike_time_jmm

end_time = min(data(d).ep,data(d).dp);
ep = find(t_axis >= end_time,1);
if ~isempty(ep)
    synct_d(ep+1:end) = []; lfp_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; lfp_hf(ep+1:end) = []; %hpc_hf(ep+1:end) = [];
    hpc_lf(ep+1:end) = []; %hpc_hf(ep+1:end) = []; wcv_hf(ep+1:end) = [];
else
    ep = length(t_axis);
end

mp_spike_rate = hist(synct(spkid),synct_d)*Fsd;
mp_spike_rate(end) = 0;

%%
load ./pa_hsmm_state_seq_combined_fin_nd.mat
load ./pa_hsmm_state_seq7_combined_fin_nd.mat
hsmm8 = hsmm7;
lfp_state_seq = hsmm_bbstate_seq7;
mp_state_seq = hsmm_bbstate_seq;

[new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,mp_state_seq);
dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
[up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
[up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);

bad_mp_states = find(up_trans_inds > ep | down_trans_inds > ep);
up_trans_inds(bad_mp_states) = []; down_trans_inds(bad_mp_states) = [];
bad_lfp_states = find(up_trans_inds8 > ep | down_trans_inds8 > ep);
up_trans_inds8(bad_lfp_states) = []; down_trans_inds8(bad_lfp_states) = [];

%%
load ./allEC_ctx_period_data_hsmm.mat
lfp_period_vec = nan(size(wcv_lf));
for i = 1:size(new_seg_inds,1)
    cur_inds = new_seg_inds(i,1):new_seg_inds(i,2);
    cur_inds_used = find(cur_inds <= ep);
    lfp_period_vec(cur_inds(cur_inds_used)) = lf8_period_f{i}(cur_inds_used);
end

%%
mp_state_number = nan(length(wcv_lf),1);
mp_state_vec = zeros(length(wcv_lf),1);
lfp_state_number = nan(length(wcv_lf),1);
lfp_state_vec = zeros(length(wcv_lf),1);
for i = 1:length(up_trans_inds)-1
    mp_state_vec(up_trans_inds(i):down_trans_inds(i)) = 1;
    mp_state_number(up_trans_inds(i):down_trans_inds(i)) = 2*(i-1)+1;
    mp_state_number(down_trans_inds(i):up_trans_inds(i+1)) = 2*(i-1) + 2;
end
for i = 1:length(up_trans_inds8)-1
    lfp_state_vec(up_trans_inds8(i):down_trans_inds8(i)) = 1;
    lfp_state_number(up_trans_inds8(i):down_trans_inds8(i)) = 2*(i-1)+1;
    lfp_state_number(down_trans_inds8(i):up_trans_inds8(i+1)) = 2*(i-1) + 2;
end
mp_state_number(isnan(lfp_period_vec)) = nan;
lfp_state_number(isnan(lfp_period_vec)) = nan;
mp_state_vec(isnan(lfp_period_vec)) = nan;
lfp_state_number(isnan(lfp_period_vec)) = nan;

%%
if exist('./mua_data3.mat','file')
    load ./mua_data3
    if usable_hpc_mua(d)
        hpc_mua_times = mua_times{peak_hpcmua_loc(d)};
        hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
        hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
        hpc_mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
    else
        hpc_mua_rate = nan(size(synct_d));
    end
    usable_ctx_mua = find(usable_mua(d,6:end));
    used_ctx_mua_chs{d} = usable_ctx_mua;
    ctx_mua_times = [];
    for ii = 1:length(usable_ctx_mua)
        ctx_mua_times = [ctx_mua_times mua_times{usable_ctx_mua(ii)+5}];
    end
    ctx_mua_times = sort(ctx_mua_times);
    ctx_mua_times(ctx_mua_times > synct_d(end) | ctx_mua_times < synct_d(1)) = [];
    ctx_mua_rate = hist(ctx_mua_times,synct_d)*Fsd;
    ctx_mua_rate = jmm_smooth_1d_cor(ctx_mua_rate,rate_sm);
    
    if isempty(usable_ctx_mua)
        ctx_mua_rate = nan(size(synct_d));
    end
    
    if length(hpc_mua_rate) > length(synct_d)
        hpc_mua_rate = hpc_mua_rate(1:length(synct_d));
        ctx_mua_rate = ctx_mua_rate(1:length(synct_d));
    end
else
    hpc_mua_rate = nan(size(synct_d));
    ctx_mua_rate = nan(size(synct_d));
end

mp_spike_rate = jmm_smooth_1d_cor(mp_spike_rate,rate_sm);

hpc_mua_rate = zscore(hpc_mua_rate);
ctx_mua_rate = zscore(ctx_mua_rate);
mp_spike_rate = zscore(mp_spike_rate);

%%
load ./aligned_heka.mat
ac_time_d = downsample(ac_time,dsf);
dc_fs = 1/median(diff(dc_time));
dc_dsf = 20;
dc_fsd = dc_fs/dc_dsf;
dc_data = decimate(dc_data,dc_dsf);
dc_time = downsample(dc_time,dc_dsf);

%correct for junction potential
jp = -7; %in mV
dc_data = dc_data + jp;
%%
wcv_up_log = nan(size(t_axis));
lfp_up_log = nan(size(t_axis));

for ns = 1:hmm.Nsegs
    cur_seg = new_seg_inds(ns,1):new_seg_inds(ns,2);
    wcv_up_log(cur_seg) = logical(mp_state_seq{ns}-1);
    lfp_up_log(cur_seg) = logical(lfp_state_seq{ns}-1);
end

%%
%%
close all
window_size = 20;
window_size_n = round(Fsd*window_size);
window_size_dc = round(dc_fsd*window_size);
n_windows = ceil(max(t_axis)/window_size);
for i = 1:n_windows
    cur_set = (i-1)*window_size_n + (1:window_size_n);
    cur_set(cur_set > length(t_axis)) = [];
    cur_set_dc = find(dc_time >= ac_time_d(cur_set(1)) & dc_time < ac_time_d(cur_set(end)));
    
    subplot(3,1,1)
    plot(dc_time(cur_set_dc),dc_data(cur_set_dc));hold on
    plot(ac_time_d(cur_set),wcv_up_log(cur_set)*20-50,'k')
    set(gca,'fontsize',10,'fontname','arial');
    box off
    
    subplot(3,1,2)
    hold on
    plot(ac_time_d(cur_set),lfp_lf(cur_set),'r','linewidth',1)
    plot(ac_time_d(cur_set),lfp_up_log(cur_set),'k')
    set(gca,'fontsize',10,'fontname','arial');
    box off
    
    subplot(3,1,3)
    hold on
    plot(ac_time_d(cur_set),ctx_mua_rate(cur_set)-5,'r')
    plot(ac_time_d(cur_set),hpc_mua_rate(cur_set)-8,'k')
    set(gca,'fontsize',10,'fontname','arial');
    box off
    xlabel('Time (s)','fontsize',14)
    ylabel('Amplitude (V)','fontsize',14)
 
    pause
    clf
end

%%
close all
bt = 385;
et = 405;
cur_set = find(ac_time_d >= bt & ac_time_d <= et);
cur_set_dc = find(dc_time >= bt & dc_time <= et);
    
    h1 = figure();
    subplot(2,1,1)
    plot(dc_time(cur_set_dc),dc_data(cur_set_dc));hold on
    plot(ac_time_d(cur_set),wcv_up_log(cur_set)*20-50,'k')
    ylabel('MP (mV)','fontsize',14);
    
    subplot(2,1,2)
    hold on
    plot(ac_time_d(cur_set),lfp_lf(cur_set),'r','linewidth',1)
    plot(ac_time_d(cur_set),lfp_up_log(cur_set),'k')
    ylabel('LFP amplitude (z)','fontsize',14);

%     subplot(3,1,3)
%     hold on
%     plot(ac_time_d(cur_set),ctx_mua_rate(cur_set)-5,'r')
%     plot(ac_time_d(cur_set),hpc_mua_rate(cur_set)-8,'k')
%     xlabel('Time (s)','fontsize',14)
%     ylabel('Firing rate (z)','fontsize',14);
 
figufy(h1);
fig_width = 4; rel_height = 1.6;
figufy(h1);
fname = [fig_dir 'examptrace_mplfp.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

