
clear all
data_dir = cell(0);
data_type = cell(0);
data_exp = cell(0);
data_ep = [];
data_dp = [];
data_hpc_lfp = [];
data_hpc_mua = [];

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
    data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3mec(ii)));
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
    data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3mec_np(ii)));
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
    data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3lec(ii)));
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
    data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3lec_np(ii)));
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
    data_hpc_mua = cat(2,data_hpc_mua,new_pdown_hpcmua(cur_uset(ii))-1);
end

addpath('C:\WC_Germany\parietal_cortical_2010\');
addpath('C:\WC_Germany\persistent_9_27_2010\');
addpath('C:\WC_Germany\new_mec\');
addpath('C:\WC_Germany\overall_EC');
addpath('C:\WC_Germany\hsmm_state_detection');
addpath('C:\WC_Germany\sven_thomas_combined');
addpath('C:\Code\general_functions\');
addpath('C:\WC_Germany\persistent_downs\');
%%

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 40;
hcf_hf = 80;
hcf_sm = 0.025;
rate_sm = round(Fsd*0.05);

backlag = 6*Fsd;
forwardlag = 8*Fsd;
lags = (-backlag:forwardlag);

params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 25;

%robust persistence requires MP skips states longer than this
thresh_lf8_updur = 0.5;
thresh_lf8_downdur = 0.5;


%%
min_rec_dur = 500; %in sec
used_dirs = find(data_ep > min_rec_dur);
for d = 1:length(used_dirs)
    cd(data_dir{used_dirs(d)})
    pwd
    
    load ./used_data lf7 wcv_minus_spike
    [lf8_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    
    if data_hpc_lfp(d) == 3
        load ./used_data lf3 lf5
        if ismember(d,old_data_inds)
            lf3 = lf3 + lf5; %redefine LF3 wrt gnd
        end
        hpc_hf = get_hf_features(lf3,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    else
        load ./used_data lf2
        hpc_hf = get_hf_features(lf2,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    end

    lf8_hf = get_hf_features(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    
    load ./sync_times.mat
    synct_d = downsample(synct,dsf);

    end_time = min(data_ep(used_dirs(d)),data_dp(used_dirs(d)));
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        lf8_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; synct_d(ep+1:end) = [];
        hpc_hf(ep+1:end) = []; lf8_hf(ep+1:end) = [];
    else
        ep = length(t_axis);
    end
    
    %% LOAD IN MUA
    if exist('./mua_data3.mat','file')
        load ./mua_data3
        if ~isnan(data_hpc_mua(used_dirs(d)))
            hpc_mua_times = mua_times{data_hpc_mua(used_dirs(d))};
            hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
            hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
            hpc_mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
        else
            hpc_mua_rate = nan(size(synct_d));
        end
        ctx_mua_times = sort([mua_times{7} mua_times{8}]);
        ctx_mua_rate = hist(ctx_mua_times,synct_d)*Fsd;
        ctx_mua_rate = jmm_smooth_1d_cor(ctx_mua_rate,rate_sm);
        if length(hpc_mua_rate) > length(synct_d)
            hpc_mua_rate = hpc_mua_rate(1:length(synct_d));
            ctx_mua_rate = ctx_mua_rate(1:length(synct_d));
        end
    else
        hpc_mua_rate = nan(size(synct_d));
        ctx_mua_rate = nan(size(synct_d));
    end
    avg_hpc_mua_rate(d) = nanmean(hpc_mua_rate);
    avg_ctx_mua_rate(d) = nanmean(ctx_mua_rate);
    std_hpc_mua_rate(d) = nanstd(hpc_mua_rate);
    std_ctx_mua_rate(d) = nanstd(ctx_mua_rate);
    hpc_mua_rate = zscore(hpc_mua_rate);
    ctx_mua_rate = zscore(ctx_mua_rate);
    
    %%
    if exist('./all_combined_mp_uds.mat','file')
        load ./all_combined_mp_uds.mat
        load ./all_combined_lf7_uds.mat
    else
        load ./pa_hsmm_state_seq7_combined_fin_nd.mat
        load ./pa_hsmm_state_seq_combined_fin_nd.mat
    end
    hmm8 = hmm7;
    lfp_state_seq = hmm_bbstate_seq7;
    mp_state_seq = hmm_bbstate_seq;
    
    %         if you want to flip the roles of the MP and LFP
    %         temp = hsmm_bbstate_seq8;
    %         hsmm_bbstate_seq8 = hsmm_bbstate_seq;
    %         hsmm_bbstate_seq = temp;
    
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    bad_mp_states = find(up_trans_inds > ep | down_trans_inds > ep);
    up_trans_inds(bad_mp_states) = []; down_trans_inds(bad_mp_states) = [];
    bad_lfp_states = find(up_trans_inds8 > ep | down_trans_inds8 > ep);
    up_trans_inds8(bad_lfp_states) = []; down_trans_inds8(bad_lfp_states) = [];
    
    fract_uds_dur(d) = dur_uds/t_axis(end);
    
    %%
    if exist('./allEC_ctx_period_data.mat','file')
        load ./allEC_ctx_period_data.mat
    else
        load ./combined_lf7_period_data_fin_nd
    end
    lf8_period_vec = nan(size(wcv_lf));
    for i = 1:size(new_seg_inds,1)
        cur_inds = new_seg_inds(i,1):new_seg_inds(i,2);
        cur_inds_used = find(cur_inds <= ep);
        lf8_period_vec(cur_inds(cur_inds_used)) = lf8_period_f{i}(cur_inds_used);
    end
    
    %%
    mp_state_number = nan(length(wcv_lf),1);
    mp_state_vec = zeros(length(wcv_lf),1);
    lfp_state_number = nan(length(wcv_lf),1);
    for i = 1:length(up_trans_inds)-1
       mp_state_vec(up_trans_inds(i):down_trans_inds(i)) = 1;
       mp_state_number(up_trans_inds(i):down_trans_inds(i)) = 2*(i-1)+1;
       mp_state_number(down_trans_inds(i):up_trans_inds(i+1)) = 2*(i-1) + 2;
    end
    for i = 1:length(up_trans_inds8)-1
       lfp_state_number(up_trans_inds8(i):down_trans_inds8(i)) = 2*(i-1)+1;
       lfp_state_number(down_trans_inds8(i):up_trans_inds8(i+1)) = 2*(i-1) + 2;
    end
    mp_state_number(isnan(lf8_period_vec)) = nan;
    lfp_state_number(isnan(lf8_period_vec)) = nan;
    mp_state_vec(isnan(lf8_period_vec)) = nan;
    
    %% compute state duration distributions
    [mp_state_durations{d}] = compute_state_durations_seg(mp_state_seq,Fsd);
    [lfp_state_durations{d}] = compute_state_durations_seg(lfp_state_seq,Fsd);
    mp_state_durations{d}{1}(bad_mp_states) = []; mp_state_durations{d}{2}(bad_mp_states) = [];
    lfp_state_durations{d}{1}(bad_lfp_states) = []; lfp_state_durations{d}{2}(bad_lfp_states) = [];
    
    mp_updurs = mp_state_durations{d}{2};
    mp_downdurs = mp_state_durations{d}{1};
    lf8_updurs= lfp_state_durations{d}{2};
    lf8_downdurs = lfp_state_durations{d}{1};
    
    %% compute corresponding state transitions and transition lags
    [corresp_lf8_upinds{d},corresp_lf8_downinds{d}] = greedy_find_corresponding_ncx_state_transitions_simp(...
        up_trans_inds,down_trans_inds,up_trans_inds8,down_trans_inds8);
    
    %find non-skipped states
    non_skipped_mp_up_states{d} = find(~isnan(corresp_lf8_upinds{d}));
    non_skipped_mp_down_states{d} = find(~isnan(corresp_lf8_downinds{d}));
    
    %compute transition lags for non-skipped states
    mp_uplags{d} = nan(size(up_trans_inds));
    mp_downlags{d} = nan(size(up_trans_inds));
    mp_uplags{d}(non_skipped_mp_up_states{d}) = up_trans_inds(non_skipped_mp_up_states{d}) - up_trans_inds8(corresp_lf8_upinds{d}(non_skipped_mp_up_states{d}));
    mp_downlags{d}(non_skipped_mp_down_states{d}) = down_trans_inds(non_skipped_mp_down_states{d}) - down_trans_inds8(corresp_lf8_downinds{d}(non_skipped_mp_down_states{d}));
    
    %% compute persistence
%     [mp_upskipped{d},mp_downskipped{d}] = greedy_find_skipped_ncx_states(...
%         corresp_lf8_upinds{d},corresp_lf8_downinds{d},lfp_state_durations{d}{2},lfp_state_durations{d}{1},thresh_lf8_downdur,thresh_lf8_updur);
    [mp_upskipped{d},mp_downskipped{d}] = greedy_find_skipped_ncx_states_v2(...
        corresp_lf8_upinds{d},corresp_lf8_downinds{d},lfp_state_durations{d}{2},lfp_state_durations{d}{1},...
        up_trans_inds8,down_trans_inds8,mp_state_vec,thresh_lf8_downdur,thresh_lf8_updur);
    
    rt2_ups{d} = find(mp_upskipped{d}.rnum_skipped > 0);
    t2_ups{d} = find(mp_upskipped{d}.num_skipped > 0);
    nrt2_ups{d} = find(mp_upskipped{d}.num_skipped == 0);
    fract_rt2_ups(d) = length(rt2_ups{d})/(length(rt2_ups{d}) + length(nrt2_ups{d}));
    fract_t2_ups(d) = length(t2_ups{d})/(length(t2_ups{d}) + length(nrt2_ups{d}));
    
    rt2_downs{d} = find(mp_downskipped{d}.rnum_skipped > 0);
    nrt2_downs{d} = find(mp_downskipped{d}.num_skipped==0); %number of down states that didn't skip any LFP up states
    fract_rt2_downs(d) = length(rt2_downs{d})/(length(rt2_downs{d}) + length(nrt2_downs{d}));
    
    robust_non_skipped_mp_ups{d} = [rt2_ups{d}; nrt2_ups{d}];
    robust_non_skipped_mp_downs{d} = [rt2_downs{d}; nrt2_downs{d}];
    
    skipped_lfp_trans = [];
    for j = 1:length(up_trans_inds)
        skipped_lfp_trans = [skipped_lfp_trans mp_upskipped{d}.inds{j}];
    end
    non_skipped_lfp_trans = setdiff(1:length(down_trans_inds8),skipped_lfp_trans);
    n_uskip(d) = length(skipped_lfp_trans);
    
    dskipped_lfp_trans = [];
    for j = 1:length(up_trans_inds)
        dskipped_lfp_trans = [dskipped_lfp_trans mp_downskipped{d}.inds{j}];
    end
    non_dskipped_lfp_trans = setdiff(1:length(up_trans_inds8),dskipped_lfp_trans);
    n_dskip(d) = length(dskipped_lfp_trans);
    
    %% CREATE INTERPOLATED MP DATA
    if exist('./aligned_heka.mat','file')
        load ./aligned_heka.mat
        dc_fs = 1/median(diff(dc_time));
        dc_dsf = 20;
        dc_fsd = dc_fs/dc_dsf;
        [dc_fb,dc_fa] = butter(4,20/(dc_fsd/2),'low');
        dc_data = decimate(dc_data,dc_dsf);
        dc_data = filtfilt(dc_fb,dc_fa,dc_data);
        dc_time = downsample(dc_time,dc_dsf);
        dc_time = dc_time(:);
        
        dc_interp_data = interp1(dc_time,dc_data,t_axis);
        ind_interp = ceil(interp1(dc_time,1:length(dc_time),t_axis));
        used = find(~isnan(ind_interp));
        t_error = dc_time(ind_interp(used)) - t_axis(used)';
        used = used(t_error <= 1/Fsd);
        
        dc_interp_data(setdiff(1:length(dc_interp_data),used)) = nan;
        if nanstd(dc_interp_data) < 1
            dc_interp_data = dc_interp_data*100;
        end
    else
        dc_interp_data = nan(size(wcv_lf));
    end
    
    %COMPUTE AVG UP- AND DOWN-STATE CONDITIONAL DC MP VALUES FOR
    %NORMALIZATION
    dc_upstate_amp(d) = nanmean(dc_interp_data(mp_state_vec == 1));
    dc_downstate_amp(d) = nanmean(dc_interp_data(mp_state_vec == 0));
    dc_avg_amp(d) = nanmean(dc_interp_data(~isnan(mp_state_vec)));
    
    %%
    min_nstates = 5;
    lf8_utrig_lf8(d,:) = get_event_trig_avg(lf8_lf,up_trans_inds8,backlag,forwardlag);
    lf8_utrig_hpc_hf(d,:) = get_event_trig_avg(hpc_hf,up_trans_inds8,backlag,forwardlag);
    lf8_utrig_ctx_hf(d,:) = get_event_trig_avg(lf8_hf,up_trans_inds8,backlag,forwardlag);
    lf8_utrig_mp(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds8,backlag,forwardlag);
    lf8_utrig_DCmp(d,:) = get_event_trig_avg(dc_interp_data,up_trans_inds8,backlag,forwardlag);
    lf8_utrig_hpc_mua(d,:) = get_event_trig_avg(hpc_mua_rate,up_trans_inds8,backlag,forwardlag);
    lf8_utrig_ctx_mua(d,:) = get_event_trig_avg(ctx_mua_rate,up_trans_inds8,backlag,forwardlag);

    mp_utrig_hpc_mua(d,:) = get_event_trig_avg(hpc_mua_rate,up_trans_inds,backlag,forwardlag);
    mp_utrig_ctx_mua(d,:) = get_event_trig_avg(ctx_mua_rate,up_trans_inds,backlag,forwardlag);

    lf8_dtrig_lf8(d,:) = get_event_trig_avg(lf8_lf,down_trans_inds8,backlag,forwardlag);
    lf8_dtrig_hpc_hf(d,:) = get_event_trig_avg(hpc_hf,down_trans_inds8,backlag,forwardlag);
    lf8_dtrig_ctx_hf(d,:) = get_event_trig_avg(lf8_hf,down_trans_inds8,backlag,forwardlag);
    lf8_dtrig_mp(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds8,backlag,forwardlag);
    lf8_dtrig_DCmp(d,:) = get_event_trig_avg(dc_interp_data,down_trans_inds8,backlag,forwardlag);
    lf8_dtrig_hpc_mua(d,:) = get_event_trig_avg(hpc_mua_rate,down_trans_inds8,backlag,forwardlag);
    lf8_dtrig_ctx_mua(d,:) = get_event_trig_avg(ctx_mua_rate,down_trans_inds8,backlag,forwardlag);

    mp_dtrig_hpc_mua(d,:) = get_event_trig_avg(hpc_mua_rate,down_trans_inds,backlag,forwardlag);
    mp_dtrig_ctx_mua(d,:) = get_event_trig_avg(ctx_mua_rate,down_trans_inds,backlag,forwardlag);

    if n_uskip(d) > min_nstates        
        lf8_utrig_lf8_pup(d,:) = get_event_trig_avg(lf8_lf,up_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_utrig_mp_pup(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_utrig_DCmp_pup(d,:) = get_event_trig_avg(dc_interp_data,up_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_lf8_pup(d,:) = get_event_trig_avg(lf8_lf,down_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_hpc_hf_pup(d,:) = get_event_trig_avg(hpc_hf,down_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_ctx_hf_pup(d,:) = get_event_trig_avg(lf8_hf,down_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_mp_pup(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_DCmp_pup(d,:) = get_event_trig_avg(dc_interp_data,down_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_hpc_mua_pup(d,:) = get_event_trig_avg(hpc_mua_rate,down_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_ctx_mua_pup(d,:) = get_event_trig_avg(ctx_mua_rate,down_trans_inds8(skipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_mp_pup_trial(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds8(skipped_lfp_trans),backlag,forwardlag,[],mp_state_number,1);
        lf8_dtrig_DCmp_pup_trial(d,:) = get_event_trig_avg(dc_interp_data,down_trans_inds8(skipped_lfp_trans),backlag,forwardlag,[],mp_state_number,1);
    end
    lf8_utrig_lf8_npup(d,:) = get_event_trig_avg(lf8_lf,up_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
    lf8_utrig_mp_npup(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
    lf8_utrig_DCmp_npup(d,:) = get_event_trig_avg(dc_interp_data,up_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
    lf8_dtrig_lf8_npup(d,:) = get_event_trig_avg(lf8_lf,down_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
    lf8_dtrig_hpc_hf_npup(d,:) = get_event_trig_avg(hpc_hf,down_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
    lf8_dtrig_ctx_hf_npup(d,:) = get_event_trig_avg(lf8_hf,down_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
    lf8_dtrig_hpc_mua_npup(d,:) = get_event_trig_avg(hpc_mua_rate,down_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
    lf8_dtrig_ctx_mua_npup(d,:) = get_event_trig_avg(ctx_mua_rate,down_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
    lf8_dtrig_mp_npup(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
    lf8_dtrig_DCmp_npup(d,:) = get_event_trig_avg(dc_interp_data,down_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag);
     lf8_dtrig_mp_npup_trial(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag,[],mp_state_number,1);
    lf8_dtrig_DCmp_npup_trial(d,:) = get_event_trig_avg(dc_interp_data,down_trans_inds8(non_skipped_lfp_trans),backlag,forwardlag,[],mp_state_number,1);
   
    lfp_uptrans_mpdown = up_trans_inds8(mp_state_vec(up_trans_inds8) == 0);
    lfp_uptrans_mpdown(ismember(lfp_uptrans_mpdown,up_trans_inds8(dskipped_lfp_trans))) = [];
    if n_dskip(d) > min_nstates
        lf8_utrig_lf8_pdown(d,:) = get_event_trig_avg(lf8_lf,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
        lf8_utrig_hpc_hf_pdown(d,:) = get_event_trig_avg(hpc_hf,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
        lf8_utrig_ctx_hf_pdown(d,:) = get_event_trig_avg(lf8_hf,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
        lf8_utrig_hpc_mua_pdown(d,:) = get_event_trig_avg(hpc_mua_rate,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
        lf8_utrig_ctx_mua_pdown(d,:) = get_event_trig_avg(ctx_mua_rate,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
        lf8_utrig_mp_pdown(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
        lf8_utrig_DCmp_pdown(d,:) = get_event_trig_avg(dc_interp_data,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
        lf8_utrig_mp_pdown_trial(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag,[],mp_state_number,1);
        lf8_utrig_DCmp_pdown_trial(d,:) = get_event_trig_avg(dc_interp_data,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag,[],mp_state_number,1);
%         [lf8_utrig_DCmp_pdown_trial(d,:),ss,nn] = get_event_trig_avg(dc_interp_data,up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag,5,mp_state_number,1);
        lf8_dtrig_lf8_pdown(d,:) = get_event_trig_avg(lf8_lf,down_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_mp_pdown(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
        lf8_dtrig_DCmp_pdown(d,:) = get_event_trig_avg(dc_interp_data,down_trans_inds8(dskipped_lfp_trans),backlag,forwardlag);
    end
    lf8_utrig_lf8_npdown(d,:) = get_event_trig_avg(lf8_lf,up_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag);
    lf8_utrig_hpc_hf_npdown(d,:) = get_event_trig_avg(hpc_hf,up_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag);
    lf8_utrig_ctx_hf_npdown(d,:) = get_event_trig_avg(lf8_hf,up_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag);
    lf8_utrig_hpc_mua_npdown(d,:) = get_event_trig_avg(hpc_mua_rate,up_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag);
    lf8_utrig_ctx_mua_npdown(d,:) = get_event_trig_avg(ctx_mua_rate,up_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag);
    lf8_utrig_mp_npdown(d,:) = get_event_trig_avg(wcv_lf,lfp_uptrans_mpdown,backlag,forwardlag);
    lf8_utrig_DCmp_npdown(d,:) = get_event_trig_avg(dc_interp_data,lfp_uptrans_mpdown,backlag,forwardlag);
%     lf8_utrig_mp_npdown_trial(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag,[],mp_state_number,1);
%     lf8_utrig_DCmp_npdown_trial(d,:) = get_event_trig_avg(dc_interp_data,up_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag,[],mp_state_number,1);
    lf8_dtrig_lf8_npdown(d,:) = get_event_trig_avg(lf8_lf,down_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag);
    lf8_dtrig_mp_npdown(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag);
    lf8_dtrig_DCmp_npdown(d,:) = get_event_trig_avg(dc_interp_data,down_trans_inds8(non_dskipped_lfp_trans),backlag,forwardlag);
    
end

cd C:\WC_Germany\persistent_downs\
% save allEC_core_analysis
save new_down_trig_avgs3 lags *trig* *skip dc_* avg_* std_*

%%
% close all
min_nstates = 15;
l3mec = find(strcmp(data_type,'L3MEC'));
pup_used = find(n_uskip > min_nstates & ismember(used_dirs,l3mec));
pdown_used = find(n_dskip > min_nstates & ismember(used_dirs,l3mec));

bad_set = [19 23];
pdown_used(ismember(pdown_used,bad_set)) = [];
pup_used(ismember(pup_used,bad_set)) = [];
% pup_used(pup_used > 67) = [];
% pdown_used(pdown_used > 67) = [];
% pup_used(pup_used <= 67) = [];
% pdown_used(pdown_used <= 67) = [];

%%
figure; 
subplot(2,1,1);hold on
h1=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_hpc_mua(pup_used,:)),nanstd(lf8_utrig_hpc_mua(pup_used,:))/sqrt(sum(~isnan(lf8_utrig_hpc_mua(pup_used,1)))),{'color','k'});
h2=shadedErrorBar(lags/Fsd,nanmean(mp_utrig_hpc_mua(pup_used,:)),nanstd(mp_utrig_hpc_mua(pup_used,:))/sqrt(sum(~isnan(mp_utrig_hpc_mua(pup_used,1)))),{'color','r'});
xlim([-2 3]);
xlabel('Time (s)','fontsize',12);
ylabel('HPC MUA rate (z)','fontsize',12);
title('Up-transition triggered','fontsize',12);
legend([h1.mainLine h2.mainLine],'Cortical state-transitions','MEC state-transitions');

subplot(2,1,2);hold on
h1=shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_hpc_mua(pup_used,:)),nanstd(lf8_dtrig_hpc_mua(pup_used,:))/sqrt(sum(~isnan(lf8_dtrig_hpc_mua(pup_used,1)))),{'color','k'});
h2=shadedErrorBar(lags/Fsd,nanmean(mp_dtrig_hpc_mua(pup_used,:)),nanstd(mp_dtrig_hpc_mua(pup_used,:))/sqrt(sum(~isnan(mp_dtrig_hpc_mua(pup_used,1)))),{'color','r'});
xlim([-2 3])
xlabel('Time (s)','fontsize',12);
ylabel('HPC MUA rate (z)','fontsize',12);
title('Down-transition triggered','fontsize',12);
legend([h1.mainLine h2.mainLine],'Cortical state-transitions','MEC state-transitions');

%% PERS UPS
figure;
% subplot(3,1,1);hold on
hold on;
h1 = shadedErrorBar(lags/Fsd,mean(lf8_dtrig_lf8_pup(pup_used,:)),std(lf8_dtrig_lf8_pup(pup_used,:))/sqrt(length(pup_used)),{'color','b'});
h2 = shadedErrorBar(lags/Fsd,mean(lf8_dtrig_lf8_npup(pup_used,:)),std(lf8_dtrig_lf8_npup(pup_used,:))/sqrt(length(pup_used)),{'color','r'});
h3 = shadedErrorBar(lags/Fsd,mean(lf8_dtrig_mp_pup(pup_used,:)),std(lf8_dtrig_mp_pup(pup_used,:))/sqrt(length(pup_used)),{'color','k'});
h4 = shadedErrorBar(lags/Fsd,mean(lf8_dtrig_mp_npup(pup_used,:)),std(lf8_dtrig_mp_npup(pup_used,:))/sqrt(length(pup_used)),{'color','g'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Ctx LFP Pers','Ctx LFP non-pers','MP pers','MP non-pers'});
xlim([-2 3])
yl = ylim(); 
line([0 0],yl,'color','k','linestyle','--');
line([-2 3],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',14);
ylabel('Amplitude (z)','fontsize',14);
title('Ctx LFP and MEC MP','fontsize',14);
fillPage(gcf,'papersize',[6 5]);
print('pers_ups_dtrig_avgs_LFPMP','-dtiff');
close

% subplot(3,1,2);hold on
figure; hold on;
h1=shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_ctx_mua_pup(pup_used,:)),nanstd(lf8_dtrig_ctx_mua_pup(pup_used,:))/sqrt(sum(~isnan(lf8_dtrig_ctx_mua_pup(pup_used,1)))),{'color','b'});
h2=shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_ctx_mua_npup(pup_used,:)),nanstd(lf8_dtrig_ctx_mua_npup(pup_used,:))/sqrt(sum(~isnan(lf8_dtrig_ctx_mua_npup(pup_used,1)))),{'color','r'});
h3=shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_hpc_mua_pup(pup_used,:)),nanstd(lf8_dtrig_hpc_mua_pup(pup_used,:))/sqrt(sum(~isnan(lf8_dtrig_hpc_mua_pup(pup_used,1)))),{'color','k'});
h4=shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_hpc_mua_npup(pup_used,:)),nanstd(lf8_dtrig_hpc_mua_npup(pup_used,:))/sqrt(sum(~isnan(lf8_dtrig_hpc_mua_npup(pup_used,1)))),{'color','g'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Ctx MUA Pers','Ctx MUA non-pers','Hpc MUA pers','Hpc MUA non-pers'});
xlim([-2 3])
yl = ylim(); 
line([0 0],yl,'color','k','linestyle','--');
line([-2 3],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',14);
ylabel('Amplitude (z)','fontsize',14);
title('Ctx MUA and Hpc MUA','fontsize',14);
fillPage(gcf,'papersize',[6 5]);
print('pers_ups_dtrig_avgs_mua','-dtiff');
close

% subplot(3,1,3);hold on
figure; hold on
h1=shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_ctx_hf_pup(pup_used,:)),nanstd(lf8_dtrig_ctx_hf_pup(pup_used,:))/sqrt(length(pup_used)),{'color','b'});
h2=shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_ctx_hf_npup(pup_used,:)),nanstd(lf8_dtrig_ctx_hf_npup(pup_used,:))/sqrt(length(pup_used)),{'color','r'});
h3=shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_hpc_hf_pup(pup_used,:)),nanstd(lf8_dtrig_hpc_hf_pup(pup_used,:))/sqrt(length(pup_used)),{'color','k'});
h4=shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_hpc_hf_npup(pup_used,:)),nanstd(lf8_dtrig_hpc_hf_npup(pup_used,:))/sqrt(length(pup_used)),{'color','g'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Ctx HF Pers','Ctx HF non-pers','Hpc HF pers','Hpc HF non-pers'});
xlim([-2 3])
yl = ylim(); 
line([0 0],yl,'color','k','linestyle','--');
line([-2 3],[0 0],'color','k','linestyle','--');
title('Ctx HFpow and Hpc HFpow','fontsize',14);
xlabel('Time lag (s)','fontsize',14);
ylabel('Amplitude (z)','fontsize',14);
fillPage(gcf,'papersize',[6 5]);
print('pers_ups_dtrig_avgs_pow','-dtiff');
close

% fillPage(gcf,'papersize',[6 14]);
% print('pers_ups_dtrig_avgs','-dtiff');
% close
%% PERS DOWNS
figure; hold on
% subplot(3,1,1);hold on
h1=shadedErrorBar(lags/Fsd,mean(lf8_utrig_lf8_pdown(pdown_used,:)),std(lf8_utrig_lf8_pdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','b'});
h2=shadedErrorBar(lags/Fsd,mean(lf8_utrig_lf8_npdown(pdown_used,:)),std(lf8_utrig_lf8_npdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','r'});
h3=shadedErrorBar(lags/Fsd,mean(lf8_utrig_mp_pdown(pdown_used,:)),std(lf8_utrig_mp_pdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','k'});
h4=shadedErrorBar(lags/Fsd,mean(lf8_utrig_mp_npdown(pdown_used,:)),std(lf8_utrig_mp_npdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','g'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Ctx LFP Pers','Ctx LFP non-pers','MP pers','MP non-pers'});
xlim([-2 3])
yl = ylim(); 
line([0 0],yl,'color','k','linestyle','--');
line([-2 3],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',14);
ylabel('Amplitude (z)','fontsize',14);
title('Ctx LFP and MEC MP','fontsize',14);
fillPage(gcf,'papersize',[6 5]);
print('pers_downs_utrig_avgs_LFPMP','-dtiff');
close

% subplot(3,1,2);hold on
figure; hold on
h1=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_ctx_mua_pdown(pdown_used,:)),nanstd(lf8_utrig_ctx_mua_pdown(pdown_used,:))/sqrt(sum(~isnan(lf8_utrig_ctx_mua_pdown(pdown_used,1)))),{'color','b'});
h2=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_ctx_mua_npdown(pdown_used,:)),nanstd(lf8_utrig_ctx_mua_npdown(pdown_used,:))/sqrt(sum(~isnan(lf8_utrig_ctx_mua_npdown(pdown_used,1)))),{'color','r'});
h3=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_hpc_mua_pdown(pdown_used,:)),nanstd(lf8_utrig_hpc_mua_pdown(pdown_used,:))/sqrt(sum(~isnan(lf8_utrig_hpc_mua_pdown(pdown_used,1)))),{'color','k'});
h4=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_hpc_mua_npdown(pdown_used,:)),nanstd(lf8_utrig_hpc_mua_npdown(pdown_used,:))/sqrt(sum(~isnan(lf8_utrig_hpc_mua_npdown(pdown_used,1)))),{'color','g'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Ctx MUA Pers','Ctx MUA non-pers','Hpc MUA pers','Hpc MUA non-pers'});
xlim([-2 3])
yl = ylim(); 
line([0 0],yl,'color','k','linestyle','--');
line([-2 3],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',14);
ylabel('Amplitude (z)','fontsize',14);
title('Ctx MUA and Hpc MUA','fontsize',14);
fillPage(gcf,'papersize',[6 5]);
print('pers_downs_utrig_avgs_mua','-dtiff');
close

% subplot(3,1,3);hold on
figure; hold on
h1=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_ctx_hf_pdown(pdown_used,:)),nanstd(lf8_utrig_ctx_hf_pdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','b'});
h2=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_ctx_hf_npdown(pdown_used,:)),nanstd(lf8_utrig_ctx_hf_npdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','r'});
h3=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_hpc_hf_pdown(pdown_used,:)),nanstd(lf8_utrig_hpc_hf_pdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','k'});
h4=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_hpc_hf_npdown(pdown_used,:)),nanstd(lf8_utrig_hpc_hf_npdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','g'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Ctx HF Pers','Ctx HF non-pers','Hpc HF pers','Hpc HF non-pers'});
xlim([-2 3])
yl = ylim(); 
line([0 0],yl,'color','k','linestyle','--');
line([-2 3],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',14);
ylabel('Amplitude (z)','fontsize',14);
title('Ctx HFpow and Hpc HFpow','fontsize',14);
fillPage(gcf,'papersize',[6 5]);
print('pers_downs_utrig_avgs_pow','-dtiff');
close

% fillPage(gcf,'papersize',[6 14]);
% print('pers_downs_utrig_avgs','-dpdf');
% close

%%
lf8_dtrig_DCmp_npup_trial_norm = bsxfun(@minus,lf8_dtrig_DCmp_npup_trial,dc_downstate_amp');
lf8_dtrig_DCmp_pup_trial_norm = bsxfun(@minus,lf8_dtrig_DCmp_pup_trial,dc_downstate_amp');
figure; hold on
% shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_DCmp_pup(pup_used,:)),nanstd(lf8_dtrig_DCmp_pup(pup_used,:))/sqrt(length(pup_used)),{'color','b'});
shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_DCmp_npup_trial_norm(pup_used,:)),nanstd(lf8_dtrig_DCmp_npup_trial_norm(pup_used,:))/sqrt(length(pup_used)),{'color','r'});
shadedErrorBar(lags/Fsd,nanmean(lf8_dtrig_DCmp_pup_trial_norm(pup_used,:)),nanstd(lf8_dtrig_DCmp_pup_trial_norm(pup_used,:))/sqrt(length(pup_used)),{'color','g'});
% plot(lags/Fsd,mean(lf8_dtrig_DCmp_npup_trial(pup_used,:)),'k');
%%
% lf8_utrig_DCmp_npdown_norm = bsxfun(@minus,lf8_utrig_DCmp_npdown,dc_downstate_amp');
lf8_utrig_DCmp_pdown_trial_norm = bsxfun(@minus,lf8_utrig_DCmp_pdown_trial,dc_downstate_amp');
% lf8_utrig_DCmp_npdown_norm = lf8_utrig_mp_npdown;
% lf8_utrig_DCmp_pdown_trial_norm = lf8_utrig_mp_pdown;

figure;hold on
% plot(lags/Fsd,lf8_utrig_DCmp_npdown_norm(pdown_used,:),'r','linewidth',0.1)
% plot(lags/Fsd,lf8_utrig_DCmp_pdown_trial_norm(pdown_used,:),'k','linewidth',0.1)
h1=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_DCmp_npdown_norm(pdown_used,:)),nanstd(lf8_utrig_DCmp_npdown_norm(pdown_used,:))/sqrt(length(pdown_used)),{'color','r'});
h2=shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_DCmp_pdown_trial_norm(pdown_used,:)),nanstd(lf8_utrig_DCmp_pdown_trial_norm(pdown_used,:))/sqrt(length(pdown_used)),{'color','k'});
% legend([h1.mainLine h2.mainLine],{'Non-persistent','Persistent'});
xlim([-2 3])
yl = ylim(); 
line([0 0],yl,'color','k','linestyle','--');
line([-2 3],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',14);
ylabel('Amplitude relative to DS (mV)','fontsize',14);

%%
for i = 1:length(pdown_used)
    i
    plot(lags/Fsd,lf8_utrig_DCmp_pdown_trial_norm(pdown_used(i),:),'k','linewidth',0.1)
hold on
plot(lags/Fsd,lf8_utrig_DCmp_npdown_trial_norm(pdown_used(i),:),'r','linewidth',0.1);
xlim([-2 3])
yl = ylim(); 
legend({'Non-persistent','Persistent'});
line([0 0],yl,'color','k','linestyle','--');
line([-2 3],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',14);
ylabel('Amplitude relative to DS (mV)','fontsize',14);

pause
clf
end
