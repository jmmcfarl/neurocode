clear all
close all
addpath('C:\WC_Germany\hsmm_state_detection\');
addpath('C:\WC_Germany\hsmm_uds_code\');
addpath('C:\WC_Germany\parietal_cortical_2010\\');


clear all
data_dir = cell(0);
data_type = cell(0);
data_exp = cell(0);
data_ep = [];
data_dp = [];
data_hpc_lfp = [];

cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
for ii = 1:length(l3mec)
    data_dir = cat(2,data_dir,combined_dir{l3mec(ii)});
    data_type = cat(2,data_type,{'L3MEC'});
    if ismember(l3mec(ii),old_data_inds)
        data_exp = cat(2,data_exp,{'T'});
    else
        data_exp = cat(2,data_exp,{'S'});
    end
    data_ep = cat(2,data_ep,Inf);
    data_dp = cat(2,data_dp,nan);
    data_hpc_lfp = cat(2,data_hpc_lfp,hpc_lfp(l3mec(ii)));
end
for ii = 1:length(l3mec_np)
    data_dir = cat(2,data_dir,combined_dir{l3mec_np(ii)});
    data_type = cat(2,data_type,{'L3MEC'});
    if ismember(l3mec_np(ii),old_data_inds)
        data_exp = cat(2,data_exp,{'T'});
    else
        data_exp = cat(2,data_exp,{'S'});
    end
    data_ep = cat(2,data_ep,Inf);
    data_dp = cat(2,data_dp,nan);
    data_hpc_lfp = cat(2,data_hpc_lfp,hpc_lfp(l3mec_np(ii)));
end
for ii = 1:length(l3lec)
    data_dir = cat(2,data_dir,combined_dir{l3lec(ii)});
    data_type = cat(2,data_type,{'L3LEC'});
    if ismember(l3lec(ii),old_data_inds)
        data_exp = cat(2,data_exp,{'T'});
    else
        data_exp = cat(2,data_exp,{'S'});
    end
    data_ep = cat(2,data_ep,Inf);
    data_dp = cat(2,data_dp,nan);
    data_hpc_lfp = cat(2,data_hpc_lfp,hpc_lfp(l3lec(ii)));
end
for ii = 1:length(l3lec_np)
    data_dir = cat(2,data_dir,combined_dir{l3lec_np(ii)});
    data_type = cat(2,data_type,{'L3LEC'});
    if ismember(l3lec_np(ii),old_data_inds)
        data_exp = cat(2,data_exp,{'T'});
    else
        data_exp = cat(2,data_exp,{'S'});
    end
    data_ep = cat(2,data_ep,Inf);
    data_dp = cat(2,data_dp,nan);
    data_hpc_lfp = cat(2,data_hpc_lfp,hpc_lfp(l3lec_np(ii)));
end

cd C:\WC_Germany\persistent_downs\
load ./new_pdown_dir
cur_uset = find(new_pdown_use == 1);
for ii = 1:length(cur_uset)
    data_dir = cat(2,data_dir,new_pdown_dir{cur_uset(ii)});
    data_type = cat(2,data_type,{'L3MEC'});
    data_exp = cat(2,data_exp,{'S'});
    data_ep = cat(2,data_ep,new_pdown_ep(cur_uset(ii)));
    data_dp = cat(2,data_ep,new_pdown_dp(cur_uset(ii)));
    data_hpc_lfp = cat(2,data_hpc_lfp,2);
end
pdown_inds = nan(67,1);
pdown_inds(68:length(data_dir)) = cur_uset;
%%
load ./new_down_core_analysis
% cd C:\WC_Germany\sven_thomas_combined\
% load ./combined_dir.mat


min_rec_dur = 500; %in sec
pdown_inds = pdown_inds(data_ep > min_rec_dur);
used_dirs = find(data_ep > min_rec_dur);
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 4;
lcf = 0.05;
to_zscore = 1;

window_size = 30;

%%
close all
%33 good example. good trace at ~1300

d = 5; 
fprintf('Testing cell %d of %d\n',d,length(new_pdown_dir));

cor_d = find(pdown_inds==d);
cur_rt2_downs = rt2_downs{cor_d};
cur_rt2_ups = rt2_ups{cor_d};

cd(new_pdown_dir{d})
% cd(combined_dir{d})

%
load ./all_combined_mp_uds.mat
load ./all_combined_lf7_uds.mat
lfp_state_seq = hmm_bbstate_seq7;
mp_state_seq = hmm_bbstate_seq;
[new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
[up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
[up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);

load ./used_data lf7 wcv_minus_spike
[lf7_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf],to_zscore);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);

load ./aligned_heka.mat
ac_time_d = downsample(ac_time,dsf);
dc_fs = 1/median(diff(dc_time));
dc_dsf = 20;
dc_fsd = dc_fs/dc_dsf;
dc_data = decimate(dc_data,dc_dsf);
dc_time = downsample(dc_time,dc_dsf);

wcv_up_log = nan(size(t_axis));
lfp_up_log = nan(size(t_axis));

for ns = 1:hmm.Nsegs
    cur_seg = new_seg_inds(ns,1):new_seg_inds(ns,2);
    wcv_up_log(cur_seg) = logical(mp_state_seq{ns}-1);
    lfp_up_log(cur_seg) = logical(lfp_state_seq{ns}-1);
end

window_size_n = round(Fsd*window_size);
window_size_dc = round(dc_fsd*window_size);
n_windows = ceil(max(t_axis)/window_size);
for i = 1:n_windows
    cur_set = (i-1)*window_size_n + (1:window_size_n);
    cur_set(cur_set > length(t_axis)) = [];
    cur_set_dc = find(dc_time >= ac_time_d(cur_set(1)) & dc_time < ac_time_d(cur_set(end)));
    subplot(2,1,1)
    plot(dc_time(cur_set_dc),dc_data(cur_set_dc))
    hold on
    plot(ac_time_d(cur_set),wcv_up_log(cur_set)*20-50,'k')
   
    plot(ac_time_d(down_trans_inds(cur_rt2_downs)),-50*ones(size(cur_rt2_downs)),'r*');
    plot(ac_time_d(up_trans_inds(cur_rt2_ups)),-50*ones(size(cur_rt2_ups)),'g*');
    xlim(dc_time(cur_set_dc([1 end])));

    
    %     subplot(3,1,[2 3])
    subplot(2,1,[2])
%     plot(ac_time_d(cur_set),wcv_lf(cur_set)+5)
    hold on
    plot(ac_time_d(cur_set),lf7_lf(cur_set),'r','linewidth',1)
%     plot(ac_time_d(cur_set),wcv_up_log(cur_set)+5,'k')
    plot(ac_time_d(cur_set),lfp_up_log(cur_set),'k')
    
    xlabel('Time (s)','fontsize',16)
    ylabel('Amplitude (V)','fontsize',16)
    
    pause
    clf
end

