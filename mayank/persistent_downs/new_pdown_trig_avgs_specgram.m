
clear all
data_dir = cell(0);
data_type = cell(0);
data_exp = cell(0);
data_ep = [];
data_dp = [];
data_hpc_lfp = [];
% data_hpc_mua = [];

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
%     data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3mec(ii)));
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
%     data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3mec_np(ii)));
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
%     data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3lec(ii)));
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
%     data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3lec_np(ii)));
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
%     data_hpc_mua = cat(2,data_hpc_mua,new_pdown_hpcmua(cur_uset(ii))-1);
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
load ./mua_classification
min_hpcrate = 1;
usable_hpc_mua = ~isnan(peak_hpcmua_loc) & peak_hpcmua_rate >= min_hpcrate;

%%

raw_Fs = 2016;
dsf = 4;
usfac = 8/dsf;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 40;
hcf_hf = 80;
hcf_sm = 0.025;
rate_sm = round(Fsd*0.05);

wo = 50/(Fsd/2);
bw = wo/50;
[b_notch,a_notch] = iirnotch(wo,bw);
wo = 100/(Fsd/2);
[b_notch2,a_notch2] = iirnotch(wo,bw);
% wo = 150/(Fsd/2);
% [b_notch3,a_notch3] = iirnotch(wo,bw);

backlag = 3*Fsd;
forwardlag = 5*Fsd;
lags = (-backlag:forwardlag);

params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 25;

%robust persistence requires MP skips states longer than this
thresh_lf8_updur = 0.5;
thresh_lf8_downdur = 0.5;

min_rec_dur = 500; %in sec
used_dirs = find(data_ep > min_rec_dur);

%%
n_freqs = 35;
min_freq = 5;
max_freq = 120;
desired_wfreqs = linspace(min_freq,max_freq,n_freqs);
scales = Fsd./(desired_wfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

%%
for d = 1:length(used_dirs)
    cd(data_dir{used_dirs(d)})
    pwd
    
    load ./used_data lf7 wcv_minus_spike
    [lf8_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    wcv_ds = zscore(decimate(wcv_minus_spike,dsf));
    
    if data_hpc_lfp(d) == 3
        load ./used_data lf3 lf5
        if ismember(d,old_data_inds)
            lf3 = lf3 + lf5; %redefine LF3 wrt gnd
        end
%         hpc_hf = get_hf_features(lf3,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
        hpc_ds = zscore(decimate(lf3,dsf));
    else
        load ./used_data lf2
%         hpc_hf = get_hf_features(lf2,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
        hpc_ds = zscore(decimate(lf2,dsf));
        hpc_ds_notch = filtfilt(b_notch,a_notch,hpc_ds);
        hpc_ds_notch = filtfilt(b_notch2,a_notch2,hpc_ds_notch);
%         hpc_ds_notch = filtfilt(b_notch3,a_notch3,hpc_ds_notch);
    end

%     lf8_hf = get_hf_features(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);    
    lf8_ds = zscore(decimate(lf7,dsf));
    lf8_ds_notch = filtfilt(b_notch,a_notch,lf8_ds);
    lf8_ds_notch = filtfilt(b_notch2,a_notch2,lf8_ds_notch);
%     lf8_ds_notch = filtfilt(b_notch3,a_notch3,lf8_ds_notch);
    
    load ./sync_times.mat
    synct_d = downsample(synct,dsf);

    end_time = min(data_ep(used_dirs(d)),data_dp(used_dirs(d)));
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        lf8_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; synct_d(ep+1:end) = [];
        hpc_ds(ep+1:end) = []; lf8_ds(ep+1:end) = []; wcv_ds(ep+1:end) = [];
    else
        ep = length(t_axis);
    end
        
    %%
    lf8_specgram = abs(cwt(lf8_ds_notch,scales,'cmor1-1'));
    hpc_specgram = abs(cwt(hpc_ds_notch,scales,'cmor1-1'));
%     mp_specgram = abs(cwt(wcv_ds,scales,'cmor1-1'));
    
    lf8_specgram = zscore(lf8_specgram');
    hpc_specgram = zscore(hpc_specgram');
%     mp_specgram = zscore(mp_specgram');
    
    z_thresh = 5;
    lf8_specgram(lf8_specgram > z_thresh) = z_thresh;
    lf8_specgram(lf8_specgram < -z_thresh) = -z_thresh;
    lf8_specgram = zscore(lf8_specgram);
    
    hpc_specgram(hpc_specgram > z_thresh) = z_thresh;
    hpc_specgram(hpc_specgram < -z_thresh) = -z_thresh;
    hpc_specgram = zscore(hpc_specgram);
    
%     mp_specgram(mp_specgram > z_thresh) = z_thresh;
%     mp_specgram(mp_specgram < -z_thresh) = -z_thresh;
%     mp_specgram = zscore(mp_specgram);
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
% load ./pa_hsmm_state_seq_combined_fin_nd.mat
% load ./pa_hsmm_state_seq7_combined_fin_nd.mat
%     hmm8 = hsmm7;
%     lfp_state_seq = hsmm_bbstate_seq7;
%     mp_state_seq = hsmm_bbstate_seq;
    for i = 1:length(lfp_state_seq)
        temp = ceil((1/usfac):(1/usfac):length(lfp_state_seq{i}));
        lfp_state_seq{i} = lfp_state_seq{i}(temp);
        mp_state_seq{i} = mp_state_seq{i}(temp);
    end
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
% load ./allEC_ctx_period_data_hsmm.mat
    for i = 1:length(lf8_period_f)
        temp = ceil((1/usfac):(1/usfac):length(lf8_period_f{i}));
        lf8_period_f{i} = lf8_period_f{i}(temp);
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
%     [corresp_lf8_upinds{d},corresp_lf8_downinds{d}] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         up_trans_inds,down_trans_inds,up_trans_inds8,down_trans_inds8);
    [corresp_lf8_upinds{d},corresp_lf8_downinds{d}] = find_corresponding_state_transitions_lookback(...
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
    
    skipped_lfp_dtrans = [];
    for j = 1:length(up_trans_inds)
        skipped_lfp_dtrans = [skipped_lfp_dtrans mp_upskipped{d}.inds{j}];
    end
    non_skipped_lfp_dtrans = setdiff(1:length(down_trans_inds8),skipped_lfp_dtrans);
    
    %don't count LFP downtrans for down states that are too short
    too_short_lfp_downs = find(lf8_downdurs < 0.5);
    skipped_lfp_dtrans(ismember(skipped_lfp_dtrans,too_short_lfp_downs)) = [];
    non_skipped_lfp_dtrans(ismember(non_skipped_lfp_dtrans,too_short_lfp_downs)) = [];
    
    n_uskip(d) = length(skipped_lfp_dtrans);
    
    skipped_lfp_utrans = [];
    for j = 1:length(up_trans_inds)
        skipped_lfp_utrans = [skipped_lfp_utrans mp_downskipped{d}.inds{j}];
    end
    non_skipped_lfp_utrans = setdiff(1:length(up_trans_inds8),skipped_lfp_utrans);
 
    %don't count LFP uptrans for up states that are too short
    too_short_lfp_ups = find(lf8_updurs < 0.5);
    skipped_lfp_utrans(ismember(skipped_lfp_utrans,too_short_lfp_ups)) = [];
    non_skipped_lfp_utrans(ismember(non_skipped_lfp_utrans,too_short_lfp_ups)) = [];

    n_dskip(d) = length(skipped_lfp_utrans);
        
    %%
    min_nstates = 5;
    lf8_utrig_lf8(d,:) = get_event_trig_avg(lf8_lf,up_trans_inds8,backlag,forwardlag);
    lf8_utrig_mp(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds8,backlag,forwardlag);

    if n_uskip(d) > min_nstates
        lf8_dtrig_lf8_pup(d,:) = get_event_trig_avg(lf8_ds,down_trans_inds8(skipped_lfp_dtrans),backlag,forwardlag);
        lf8_dtrig_lf8_spec_pup(d,:,:) = get_event_trig_avg(lf8_specgram,down_trans_inds8(skipped_lfp_dtrans),backlag,forwardlag);
        lf8_dtrig_hpc_pup(d,:) = get_event_trig_avg(hpc_ds,down_trans_inds8(skipped_lfp_dtrans),backlag,forwardlag);
        lf8_dtrig_hpc_spec_pup(d,:,:) = get_event_trig_avg(hpc_specgram,down_trans_inds8(skipped_lfp_dtrans),backlag,forwardlag);
%         lf8_dtrig_mp_pup(d,:) = get_event_trig_avg(wcv_ds,down_trans_inds8(skipped_lfp_dtrans),backlag,forwardlag);
%         lf8_dtrig_mp_spec_pup(d,:,:) = get_event_trig_avg(mp_specgram,down_trans_inds8(skipped_lfp_dtrans),backlag,forwardlag);
    end
    lf8_dtrig_lf8_npup(d,:) = get_event_trig_avg(lf8_ds,down_trans_inds8(non_skipped_lfp_dtrans),backlag,forwardlag);
    lf8_dtrig_lf8_spec_npup(d,:,:) = get_event_trig_avg(lf8_specgram,down_trans_inds8(non_skipped_lfp_dtrans),backlag,forwardlag);
    lf8_dtrig_hpc_npup(d,:) = get_event_trig_avg(hpc_ds,down_trans_inds8(non_skipped_lfp_dtrans),backlag,forwardlag);
    lf8_dtrig_hpc_spec_npup(d,:,:) = get_event_trig_avg(hpc_specgram,down_trans_inds8(non_skipped_lfp_dtrans),backlag,forwardlag);
%     lf8_dtrig_mp_npup(d,:) = get_event_trig_avg(wcv_ds,up_trans_inds8(non_skipped_lfp_dtrans),backlag,forwardlag);
%     lf8_dtrig_mp_spec_npup(d,:,:) = get_event_trig_avg(mp_specgram,up_trans_inds8(non_skipped_lfp_dtrans),backlag,forwardlag);
   
    if n_dskip(d) > min_nstates
         lf8_utrig_lf8_pdown(d,:) = get_event_trig_avg(lf8_ds,up_trans_inds8(skipped_lfp_utrans),backlag,forwardlag);
        lf8_utrig_lf8_spec_pdown(d,:,:) = get_event_trig_avg(lf8_specgram,up_trans_inds8(skipped_lfp_utrans),backlag,forwardlag);
        lf8_utrig_hpc_pdown(d,:) = get_event_trig_avg(hpc_ds,up_trans_inds8(skipped_lfp_utrans),backlag,forwardlag);
        lf8_utrig_hpc_spec_pdown(d,:,:) = get_event_trig_avg(hpc_specgram,up_trans_inds8(skipped_lfp_utrans),backlag,forwardlag);
%         lf8_utrig_mp_pdown(d,:) = get_event_trig_avg(wcv_ds,up_trans_inds8(skipped_lfp_utrans),backlag,forwardlag);
%         lf8_utrig_mp_spec_pdown(d,:,:) = get_event_trig_avg(mp_specgram,up_trans_inds8(skipped_lfp_utrans),backlag,forwardlag);
    end
    lf8_utrig_lf8_npdown(d,:) = get_event_trig_avg(lf8_ds,up_trans_inds8(non_skipped_lfp_utrans),backlag,forwardlag);
    lf8_utrig_lf8_spec_npdown(d,:,:) = get_event_trig_avg(lf8_specgram,up_trans_inds8(non_skipped_lfp_utrans),backlag,forwardlag);
    lf8_utrig_hpc_npdown(d,:) = get_event_trig_avg(hpc_ds,up_trans_inds8(non_skipped_lfp_utrans),backlag,forwardlag);
    lf8_utrig_hpc_spec_npdown(d,:,:) = get_event_trig_avg(hpc_specgram,up_trans_inds8(non_skipped_lfp_utrans),backlag,forwardlag);
%     lf8_utrig_mp_npdown(d,:) = get_event_trig_avg(wcv_ds,up_trans_inds8(non_skipped_lfp_utrans),backlag,forwardlag);
%     lf8_utrig_mp_spec_npdown(d,:,:) = get_event_trig_avg(mp_specgram,up_trans_inds8(non_skipped_lfp_utrans),backlag,forwardlag);
 
end

cd C:\WC_Germany\persistent_downs\
% save allEC_core_analysis
save new_down_trig_avgs_spec3 lags *trig* *skip avg_* wfreqs

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
cd C:\persDowns_paper\Figs\
close all

figure; 
subplot(2,1,1);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean((lf8_utrig_lf8_spec_npdown(pdown_used,:,:))))');shading interp
xlim([-0.25 1.]);
ylim([5 100]);
xlabel('Time (s)','fontsize',12);
ylabel('Frequency (Hz)','fontsize',12);
title('No Persistent down','fontsize',12);
line([-2 3],[0 0],'color','k');
ca = caxis();
colorbar
subplot(2,1,2);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean((lf8_utrig_lf8_spec_pdown(pdown_used,:,:))))');shading interp
xlim([-0.25 1.]);
ylim([5 100]);
xlabel('Time (s)','fontsize',12);
ylabel('Frequency (Hz)','fontsize',12);
title('Persistent down','fontsize',12);
line([-2 3],[0 0],'color','k');
caxis(ca);
colorbar
fillPage(gcf,'papersize',[6 10]);
% print('utrig_avg_ctx_spec','-dpdf');                                                                                                                                                                                                                                                                                                                                                                             ','-dpdf');

% figure
% hold on
% shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_lf8_npdown(pdown_used,:)),nanstd(lf8_utrig_lf8_npdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','r'});
% shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_lf8_pdown(pdown_used,:)),nanstd(lf8_utrig_lf8_pdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','b'});
% xlim([-0.5 1.25]);
% xlabel('Time (s)','fontsize',12);
% ylabel('Amplitude (z)','fontsize',12);
% % lf8_npdown = squeeze(nanmean(lf8_utrig_lf8_spec_npdown(pdown_used,:,:)))';
% % lf8_pdown = squeeze(nanmean(lf8_utrig_lf8_spec_pdown(pdown_used,:,:)))';
% % lf8_dom = (lf8_npdown - lf8_pdown);
% % figure
% % pcolor(lags/Fsd,wfreqs,lf8_dom);shading interp
% % xlim([-0.5 1.25]);
% % colorbar
%%
figure; 
subplot(2,1,1);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean(lf8_utrig_hpc_spec_npdown(pdown_used,:,:)))');shading interp
xlim([-0.25 1.]);
ylim([5 100]);
xlabel('Time (s)','fontsize',12);
ylabel('HPC MUA rate (z)','fontsize',12);
title('Up-transition triggered','fontsize',12);
line([-2 3],[0 0],'color','k');
ca = caxis();
colorbar
subplot(2,1,2);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean(lf8_utrig_hpc_spec_pdown(pdown_used,:,:)))');shading interp
xlim([-0.25 1.]);
ylim([5 100]);
xlabel('Time (s)','fontsize',12);
ylabel('HPC MUA rate (z)','fontsize',12);
title('Up-transition triggered','fontsize',12);
line([-2 3],[0 0],'color','k');
caxis(ca);
colorbar
fillPage(gcf,'papersize',[6 10]);
print('utrig_avg_hpc_spec','-dpdf');                                                                                                                                                                                                                                                                                                                                                                            

% % hpc_npdown = squeeze(nanmean(lf8_utrig_hpc_spec_npdown(pdown_used,:,:)))';
% % hpc_pdown = squeeze(nanmean(lf8_utrig_hpc_spec_pdown(pdown_used,:,:)))';
% % hpc_dom = (hpc_npdown - hpc_pdown);
% % figure
% % pcolor(lags/Fsd,wfreqs,hpc_dom);shading interp
% % xlim([-0.5 1.25]);
% % colorbar
% 
% figure
% hold on
% shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_hpc_npdown(pdown_used,:)),nanstd(lf8_utrig_hpc_npdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','r'});
% shadedErrorBar(lags/Fsd,nanmean(lf8_utrig_hpc_pdown(pdown_used,:)),nanstd(lf8_utrig_hpc_pdown(pdown_used,:))/sqrt(length(pdown_used)),{'color','b'});
% xlim([-0.5 1.25]);

%%
figure; 
subplot(2,1,1);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean((lf8_dtrig_lf8_spec_npup(pup_used,:,:))))');shading interp
xlim([-0.25 0.75]);
ylim([5 100]);
xlabel('Time (s)','fontsize',12);
ylabel('Frequency (Hz)','fontsize',12);
title('No Persistent down','fontsize',12);
line([-2 3],[0 0],'color','k');
caxis([-0.6 0.4])
colorbar
subplot(2,1,2);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean((lf8_dtrig_lf8_spec_pup(pup_used,:,:))))');shading interp
xlim([-0.25 0.75]);
ylim([5 100]);
xlabel('Time (s)','fontsize',12);
ylabel('Frequency (Hz)','fontsize',12);
title('Persistent down','fontsize',12);
line([-2 3],[0 0],'color','k');
caxis([-0.6 0.4])
colorbar
fillPage(gcf,'papersize',[6 10]);
% print('dtrig_avg_ctx_spec','-dpdf');

%%
figure; 
subplot(2,1,1);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean((lf8_dtrig_hpc_spec_npup(pup_used,:,:))))');shading interp
xlim([-0.25 0.75]);
ylim([5 100]);
xlabel('Time (s)','fontsize',12);
ylabel('Frequency (Hz)','fontsize',12);
title('No Persistent down','fontsize',12);
line([-2 3],[0 0],'color','k');
caxis([-0.4 0.4])
colorbar
subplot(2,1,2);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean((lf8_dtrig_hpc_spec_pup(pup_used,:,:))))');shading interp
xlim([-0.25 0.75]);
ylim([5 100]);
xlabel('Time (s)','fontsize',12);
ylabel('Frequency (Hz)','fontsize',12);
title('Persistent down','fontsize',12);
line([-2 3],[0 0],'color','k');
caxis([-0.4 0.4])
colorbar
fillPage(gcf,'papersize',[6 10]);
% print('dtrig_avg_hpc_spec','-dpdf');

%%
figure; 
subplot(2,1,1);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean(lf8_utrig_mp_spec_npdown(pdown_used,:,:)))');shading interp
xlim([-0.5 1.25]);
xlabel('Time (s)','fontsize',12);
ylabel('HPC MUA rate (z)','fontsize',12);
title('Up-transition triggered','fontsize',12);
line([-2 3],[0 0],'color','k');
ylim([5 200]);
ca = caxis();
subplot(2,1,2);hold on
pcolor(lags/Fsd,wfreqs,squeeze(nanmean(lf8_utrig_mp_spec_pdown(pdown_used,:,:)))');shading interp
xlim([-0.5 1.25]);
xlabel('Time (s)','fontsize',12);
ylabel('HPC MUA rate (z)','fontsize',12);
title('Up-transition triggered','fontsize',12);
line([-2 3],[0 0],'color','k');
caxis(ca);
ylim([5 200]);

