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
load ./mua_classification
min_hpcrate = 5;
usable_hpc_mua = ~isnan(peak_hpcmua_loc) & peak_hpcmua_rate >= min_hpcrate;
%%

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.01;
lcf_hf = 40;
hcf_hf = 100;
hcf_sm = 0.025;
rate_sm = round(Fsd*0.05);

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
for d = 1 :length(used_dirs)
    cd(data_dir{used_dirs(d)})
    pwd
    
    load ./used_data lf7 wcv_minus_spike
    [lf8_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    lf8_hf = get_hf_features(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);

    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    
    end_time = min(data_ep(used_dirs(d)),data_dp(used_dirs(d)));
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        lf8_hf(ep+1:end) = [];
        lf8_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; synct_d(ep+1:end) = [];
    else
        ep = length(t_axis);
    end
    
    %%
%     if exist('./all_combined_mp_uds.mat','file')
%         load ./all_combined_mp_uds.mat
%         load ./all_combined_lf7_uds.mat
%     else
%         load ./pa_hsmm_state_seq7_combined_fin_nd.mat
%         load ./pa_hsmm_state_seq_combined_fin_nd.mat
%     end
load ./pa_hsmm_state_seq_combined_fin_nd.mat
load ./pa_hsmm_state_seq7_combined_fin_nd.mat
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
%     if exist('./allEC_ctx_period_data.mat','file')
%         load ./allEC_ctx_period_data.mat
%     else
%         load ./combined_lf7_period_data_fin_nd
%     end
load ./allEC_ctx_period_data_hsmm.mat
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
    mp_state_number(isnan(lf8_period_vec)) = nan;
    lfp_state_number(isnan(lf8_period_vec)) = nan;
    mp_state_vec(isnan(lf8_period_vec)) = nan;
    lfp_state_number(isnan(lf8_period_vec)) = nan;
    
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
%     [corresp_mp_upinds{d},corresp_mp_downinds{d}] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         up_trans_inds8,down_trans_inds8,up_trans_inds,down_trans_inds);
    [corresp_lf8_upinds{d},corresp_lf8_downinds{d}] = find_corresponding_state_transitions_lookback(...
        up_trans_inds,down_trans_inds,up_trans_inds8,down_trans_inds8);
    [corresp_mp_upinds{d},corresp_mp_downinds{d}] = find_corresponding_state_transitions_lookback(...
        up_trans_inds8,down_trans_inds8,up_trans_inds,down_trans_inds);
    
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
    t2_downs{d} = find(mp_downskipped{d}.num_skipped > 0);
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
if exist('./mua_data3.mat','file')
    load ./mua_data3
    
    usable_ctx_mua = find(usable_mua(used_dirs(d),6:end));
    used_ctx_mua_chs{d} = usable_ctx_mua;
    ctx_mua_times = [];
    for ii = 1:length(usable_ctx_mua)
        ctx_mua_times = [ctx_mua_times mua_times{usable_ctx_mua(ii)+5}];
    end
    ctx_mua_times = sort(ctx_mua_times);
    ctx_mua_rate = hist(ctx_mua_times,synct_d)*Fsd;
    ctx_mua_rate = jmm_smooth_1d_cor(ctx_mua_rate,rate_sm);
    
    if isempty(usable_ctx_mua)
        ctx_mua_rate = nan(size(synct_d));
    end  
    
%     ctx_mua_times = sort([mua_times{6} mua_times{7} mua_times{8}]);
%     ctx_mua_rate = hist(ctx_mua_times,synct_d)*Fsd;
%     ctx_mua_rate = jmm_smooth_1d_cor(ctx_mua_rate,rate_sm);
    if length(ctx_mua_rate) > length(synct_d)
        ctx_mua_rate = ctx_mua_rate(1:length(synct_d));
    end
else
    ctx_mua_rate = nan(size(synct_d));
end
avg_ctx_mua_rate(d) = nanmean(ctx_mua_rate);
std_ctx_mua_rate(d) = nanstd(ctx_mua_rate);
ctx_mua_rate = zscore(ctx_mua_rate);

    %%
%     mp_up_amps = nan(length(up_trans_inds),1);
%     mp_down_amps = nan(length(up_trans_inds),1);
%     for ii = 1:length(up_trans_inds)-1
%         cur_up_inds = up_trans_inds(ii):down_trans_inds(ii);
%         cur_down_inds = down_trans_inds(ii):up_trans_inds(ii+1);
%         mp_up_amps(ii) = nanmean(dc_interp_data_normd(cur_up_inds));
%         mp_down_amps(ii) = nanmean(dc_interp_data_normd(cur_down_inds));
% %         mp_up_amps(ii) = prctile(dc_interp_data_normd(cur_up_inds),75);
% %         mp_down_amps(ii) = prctile(dc_interp_data_normd(cur_down_inds),75);
%     end
    
    [ctx_up_amps{d},down_amps] = get_state_amplitudes(lf8_lf,up_trans_inds8,down_trans_inds8);
    [ctx_up_amps_hf{d},down_amps_hf] = get_state_amplitudes(lf8_hf,up_trans_inds8,down_trans_inds8);
    [ctx_up_amps_mua{d},down_amps_mua] = get_state_amplitudes(ctx_mua_rate,up_trans_inds8,down_trans_inds8);

    %normalize up-state amplitudes
    ctx_up_amps{d} = nanzscore(ctx_up_amps{d});
    ctx_up_amps_hf{d} = nanzscore(ctx_up_amps_hf{d});
    ctx_up_amps_mua{d} = nanzscore(ctx_up_amps_mua{d});
    
    ctx_down_amps{d} = nanzscore(down_amps);
    ctx_down_amps_hf{d} = nanzscore(down_amps_hf);
    ctx_down_amps_mua{d} = nanzscore(down_amps_mua);
    
    ctx_up_durs{d} = lfp_state_durations{d}{2};
    
    uset = find(~isnan(corresp_mp_upinds{d}));
    corresp_mp_uplag{d} = nan(length(up_trans_inds8),1);
    corresp_mp_uplag{d}(uset) = mp_uplags{d}(corresp_mp_upinds{d}(uset))/Fsd;
    
    uset = find(~isnan(corresp_mp_upinds{d}));
    corresp_mp_downlag{d} = nan(length(up_trans_inds8),1);
    corresp_mp_downlag{d}(uset) = mp_downlags{d}(corresp_mp_upinds{d}(uset))/Fsd;

    corresp_mp_amps{d} = nan(length(up_trans_inds8),1);

    prev_ctx_downdur{d} = nan(length(up_trans_inds8),1);
    prev_ctx_downdur{d}(2:end) = lfp_state_durations{d}{1}(1:end-1);
    
    median_uplag = round(nanmedian(mp_uplags{d}));
    ctx_up_isskipped{d} = nan(length(up_trans_inds8),1);
    %     for ii = 1:length(rt2_downs{d})
    %         cur_skipped = mp_downskipped{d}.inds{rt2_downs{d}(ii)};
    %         ctx_up_isskipped{d}(cur_skipped) = 1;
    %     end
    %     for ii = 1:length(up_trans_inds)
    %         if ~isnan(corresp_lf8_upinds{d}(ii))
    %             ctx_up_isskipped{d}(corresp_lf8_upinds{d}(ii)) = 0;
    %         end
    %     end
    ctx_up_isskipped{d}(skipped_lfp_utrans) = 1;
    ctx_up_isskipped{d}(non_skipped_lfp_utrans) = 0;
    ctx_ups_mp_up = find(mp_state_vec(up_trans_inds8) == 1);
    ctx_up_isskipped{d}(ctx_ups_mp_up) = nan;
    
   ctx_down_isskipped{d} = nan(length(down_trans_inds8),1);
%     for ii = 1:length(rt2_ups{d})
%         cur_skipped = mp_upskipped{d}.inds{rt2_ups{d}(ii)};
%         ctx_down_isskipped{d}(cur_skipped) = 1;
%     end
%     for ii = 1:length(down_trans_inds)
%         if ~isnan(corresp_lf8_downinds{d}(ii))
%             ctx_down_isskipped{d}(corresp_lf8_downinds{d}(ii)) = 0;
%         end
%     end
    ctx_down_isskipped{d}(skipped_lfp_dtrans) = 1;
    ctx_down_isskipped{d}(non_skipped_lfp_dtrans) = 0;
    ctx_downs_mp_down = find(mp_state_vec(down_trans_inds8) == 0);
    ctx_down_isskipped{d}(ctx_ups_mp_up) = nan;
    
    ctx_up_ispers{d} = nan(length(up_trans_inds8),1);
    for i = 1:length(up_trans_inds)
        if ~isnan(corresp_lf8_upinds{d}(i))
            if ismember(i,rt2_ups{d})
                ctx_up_ispers{d}(corresp_lf8_upinds{d}(i)) = 1;
            else
                ctx_up_ispers{d}(corresp_lf8_upinds{d}(i)) = 0;
            end
        end
    end

    time_since_last_mp_down{d} = nan(size(up_trans_inds8));
    for i = 1:length(time_since_last_mp_down{d})
%         if ~ismember(i,skipped_lfp_trans)
            prev_mp_down = find(down_trans_inds < up_trans_inds8(i),1,'last');
            if ~isempty(prev_mp_down)
                time_since_last_mp_down{d}(i) = (up_trans_inds8(i) - down_trans_inds(prev_mp_down))/Fsd;
            end
%         end
    end

    time_since_last_mp_up{d} = nan(size(up_trans_inds8));
    for i = 1:length(time_since_last_mp_up{d})
        %         if ~ismember(i,skipped_lfp_trans)
        prev_mp_up = find(up_trans_inds < down_trans_inds8(i),1,'last');
        if ~isempty(prev_mp_up)
            time_since_last_mp_up{d}(i) = (down_trans_inds8(i) - up_trans_inds(prev_mp_up))/Fsd;
        end
        %         end
    end

    %%
    not_skipped = ctx_up_isskipped{d} == 0;
    
    [uplag_ctxamp_c(d),uplag_ctxamp_p(d)] = nancorr(corresp_mp_uplag{d}(not_skipped),ctx_up_amps{d}(not_skipped),'spearman');
    [uplag_ctxamphf_c(d),uplag_ctxamphf_p(d)] = nancorr(corresp_mp_uplag{d}(not_skipped),ctx_up_amps_hf{d}(not_skipped),'spearman');
    [downlag_ctxamp_c(d),downlag_ctxamp_p(d)] = nancorr(corresp_mp_downlag{d}(not_skipped),ctx_up_amps{d}(not_skipped),'spearman');
    [downlag_ctxamphf_c(d),downlag_ctxamphf_p(d)] = nancorr(corresp_mp_downlag{d}(not_skipped),ctx_up_amps_hf{d}(not_skipped),'spearman');
    [prvdowndur_ctxamp_c(d),prvdowndur_ctxamp_p(d)] = nancorr(prev_ctx_downdur{d},ctx_up_amps{d},'spearman');
    [prvdowndur_ctxamphf_c(d),prvdowndur_ctxamphf_p(d)] = nancorr(prev_ctx_downdur{d},ctx_up_amps_hf{d},'spearman');
   
    
    [b,dev,stats] = glmfit(ctx_up_amps_hf{d},ctx_up_isskipped{d},'binomial');
    pdown_ctxamphf_p(d) = stats.p(2);
    pdown_ctxamphf_t(d) = stats.t(2);
    [b,dev,stats] = glmfit(ctx_up_amps{d},ctx_up_isskipped{d},'binomial');
    pdown_ctxamp_p(d) = stats.p(2);
    pdown_ctxamp_t(d) = stats.t(2);
%     [b,dev,stats] = glmfit(corresp_mp_uplag{d},ctx_up_isskipped{d},'binomial');
%     pdown_uplag_p(d) = stats.p(2);
%     pdown_uplag_t(d) = stats.t(2);
    [b,dev,stats] = glmfit(time_since_last_mp_down{d},ctx_up_isskipped{d},'binomial');
    pdown_tslu_p(d) = stats.p(2);
    pdown_tslu_t(d) = stats.t(2);
    [b,dev,stats] = glmfit(prev_ctx_downdur{d},ctx_up_isskipped{d},'binomial');
    pdown_prvctxdown_p(d) = stats.p(2);
    pdown_prvctxdown_t(d) = stats.t(2);
 
    [b,dev,stats] = glmfit([time_since_last_mp_down{d} ctx_up_amps_hf{d}],ctx_up_isskipped{d},'binomial');
    pdown_tslu_joint_p(d) = stats.p(2);
    pdown_upamphf_joint_p(d) = stats.p(3);
    pdown_tslu_joint_t(d) = stats.t(2);
    pdown_upamphf_joint_t(d) = stats.t(3);
    
    if ~isnan(avg_ctx_mua_rate(d))
    [b,dev,stats] = glmfit(ctx_up_amps_mua{d},ctx_up_isskipped{d},'binomial');
    pdown_ctxampmua_p(d) = stats.p(2);
    [uplag_ctxampmua_c(d),uplag_ctxampmua_p(d)] = nancorr(corresp_mp_uplag{d},ctx_up_amps_mua{d},'spearman');
    [downlag_ctxampmua_c(d),downlag_ctxampmua_p(d)] = nancorr(corresp_mp_downlag{d},ctx_up_amps_mua{d},'spearman');
   [prvdowndur_ctxampmua_c(d),prvdowndur_ctxampmua_p(d)] = nancorr(prev_ctx_downdur{d},ctx_up_amps_mua{d},'spearman');
    else
       pdown_ctxampmua_p(d) = nan;
       uplag_ctxampmua_c(d) = nan;
       uplag_ctxampmua_p(d) = nan;
       downlag_ctxampmua_c(d) = nan;
       downlag_ctxampmua_p(d) = nan;
       prvdowndur_ctxampmua_c(d) = nan;
       prvdowndur_ctxampmua_p(d) = nan;
    end
    
    %%
        [b,dev,stats] = glmfit(ctx_up_amps_hf{d},ctx_up_ispers{d},'binomial');
    pup_ctxamphf_p(d) = stats.p(2);
    pup_ctxamphf_t(d) = stats.t(2);
    [b,dev,stats] = glmfit(ctx_up_amps{d},ctx_up_ispers{d},'binomial');
    pup_ctxamp_p(d) = stats.p(2);
    pup_ctxamp_t(d) = stats.t(2);
        [b,dev,stats] = glmfit(ctx_down_amps_hf{d},ctx_down_isskipped{d},'binomial');
    pup_damphf_p(d) = stats.p(2);
    pup_damphf_t(d) = stats.t(2);
    [b,dev,stats] = glmfit(ctx_down_amps{d},ctx_down_isskipped{d},'binomial');
    pup_damp_p(d) = stats.p(2);
    pup_damp_t(d) = stats.t(2);
    [b,dev,stats] = glmfit(corresp_mp_uplag{d},ctx_up_ispers{d},'binomial');
    pup_uplag_p(d) = stats.p(2);
    pup_uplag_t(d) = stats.t(2);
    [b,dev,stats] = glmfit(time_since_last_mp_up{d},ctx_up_ispers{d},'binomial');
    pup_tslu_p(d) = stats.p(2);
    pup_tslu_t(d) = stats.t(2);
    [b,dev,stats] = glmfit(prev_ctx_downdur{d},ctx_up_ispers{d},'binomial');
    pup_prvctxdown_p(d) = stats.p(2);
     pup_prvctxdown_t(d) = stats.t(2);

    [b,dev,stats] = glmfit([ctx_down_amps_hf{d} ctx_down_amps{d}],ctx_down_isskipped{d},'binomial');
    pup_damp_joint_p(d) = stats.p(2);
    pup_damphf_joint_p(d) = stats.p(3);
     pup_damp_joint_t(d) = stats.t(2);
    pup_damphf_joint_t(d) = stats.t(3);
   
    if ~isnan(avg_ctx_mua_rate(d))
        [b,dev,stats] = glmfit(ctx_up_amps_mua{d},ctx_up_ispers{d},'binomial');
        pup_ctxampmua_p(d) = stats.p(2);
        pup_ctxampmua_t(d) = stats.t(2);
        [b,dev,stats] = glmfit(ctx_down_amps_mua{d},ctx_down_isskipped{d},'binomial');
        pup_dampmua_p(d) = stats.p(2);
        pup_dampmua_t(d) = stats.t(2);
    else
        pup_ctxampmua_p(d) = nan;
         pup_ctxampmua_t(d) = nan;
        pup_dampmua_p(d) = nan;
         pup_dampmua_t(d) = nan;
    end
    
   %%
   
   
end

%%
cd C:\WC_Germany\persistent_downs\
save new_down_trig_transfun_norm_newdefs_lf2 Fsd ctx_up_isskipped pdown_* prvdown* pup* down* *skip uplag_ctx* corresp_mp_* ctx_up_* ctx_down_* time_since_last_mp_* fract* prev_*

%%
 min_nstates = 15;
l3mec = find(strcmp(data_type(used_dirs),'L3MEC'));
l3lec = find(strcmp(data_type(used_dirs),'L3LEC'));
pup_used = find(n_uskip > min_nstates & ismember(used_dirs,l3mec));
pdown_used = find(n_dskip > min_nstates & ismember(used_dirs,l3mec));

%%
% load new_down_core_analysis
%%
cd C:\WC_Germany\persistent_downs\
load ./new_down_trig_transfun_norm_newdefs

close all
min_pdown = 0.05;
used_pdowns = l3mec(fract_rt2_downs(l3mec) >= min_pdown);
% used_pdowns = used_pdowns(avg_ctx_mua_rate(used_pdowns) > 0.1);

nbins = 20;
pbins = linspace(0,100,nbins+1);
ctxhf_cond_pdown = nan(length(used_pdowns),nbins);
ctx_cond_pdown = nan(length(used_pdowns),nbins);
ctxdur_cond_pdown = nan(length(used_pdowns),nbins);
ctxmua_cond_pdown = nan(length(used_pdowns),nbins);
ctx_cond_uplag = nan(length(used_pdowns),nbins);
ctxmua_cond_uplag = nan(length(used_pdowns),nbins);
ctxtsu_cond_pdown = nan(length(used_pdowns),nbins);
ctxtsu_cond_uplag = nan(length(used_pdowns),nbins);
ctxtsu_cond_upamp = nan(length(used_pdowns),nbins);
ctxtsu_cond_upamphf = nan(length(used_pdowns),nbins);
ctxtsu_cond_upampmua = nan(length(used_pdowns),nbins);
ctxpdd_cond_pdown = nan(length(used_pdowns),nbins);
ctxpdd_cond_uplag = nan(length(used_pdowns),nbins);
ctxpdd_cond_upamp = nan(length(used_pdowns),nbins);
ctxpdd_cond_upamphf = nan(length(used_pdowns),nbins);
ctxpdd_cond_upampmua = nan(length(used_pdowns),nbins);
for ii = 1:length(used_pdowns)
    cur_bin_edges = prctile(ctx_up_amps_hf{used_pdowns(ii)},pbins);
    [nb,ni] = histc(ctx_up_amps_hf{used_pdowns(ii)},cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
           uset = find(ni==jj);
           cur_pdown = nansum(ctx_up_isskipped{used_pdowns(ii)}(uset)==1);
           cur_npdown = nansum(ctx_up_isskipped{used_pdowns(ii)}(uset)==0);           
%           cur_tot = sum(~isnan(ctx_up_isskipped{used_pdowns(ii)}(ni==jj)));
          ctxhf_cond_pdown(ii,jj) = cur_pdown/(cur_pdown+cur_npdown);
          
          uset = uset(ctx_up_isskipped{used_pdowns(ii)}(uset)==0);
          ctxhf_cond_uplag(ii,jj) = nanmean(corresp_mp_uplag{used_pdowns(ii)}(uset));
       end
    end
    cur_npdown = nansum(ctx_up_isskipped{used_pdowns(ii)}==1);
    cur_tot = sum(~isnan(ctx_up_isskipped{used_pdowns(ii)}));
    test(ii) = cur_npdown/cur_tot;
    
    cur_bin_edges = prctile(ctx_up_amps{used_pdowns(ii)},pbins);
    [nb,ni] = histc(ctx_up_amps{used_pdowns(ii)},cur_bin_edges);
    for jj = 1:nbins
        if nb(jj) > 0
            uset = find(ni==jj);
            cur_pdown = nansum(ctx_up_isskipped{used_pdowns(ii)}(uset)==1);
            cur_npdown = nansum(ctx_up_isskipped{used_pdowns(ii)}(uset)==0);
            %%           cur_tot = sum(~isnan(ctx_up_isskipped{used_pdowns(ii)}(ni==jj)));
            ctx_cond_pdown(ii,jj) = cur_pdown/(cur_pdown+cur_npdown);
            
            uset = uset(ctx_up_isskipped{used_pdowns(ii)}(uset)==0);
            ctx_cond_uplag(ii,jj) = nanmean(corresp_mp_uplag{used_pdowns(ii)}(uset));
        end
    end
    cur_bin_edges = prctile(ctx_up_durs{used_pdowns(ii)},pbins);
    [nb,ni] = histc(ctx_up_durs{used_pdowns(ii)},cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
           uset = find(ni==jj);
           cur_pdown = nansum(ctx_up_isskipped{used_pdowns(ii)}(uset)==1);
           cur_npdown = nansum(ctx_up_isskipped{used_pdowns(ii)}(uset)==0);           
%           cur_tot = sum(~isnan(ctx_up_isskipped{used_pdowns(ii)}(ni==jj)));
          ctxdur_cond_pdown(ii,jj) = cur_pdown/(cur_pdown+cur_npdown);
       end
    end
    cur_bin_edges = prctile(ctx_up_amps_mua{used_pdowns(ii)},pbins);
    [nb,ni] = histc(ctx_up_amps_mua{used_pdowns(ii)},cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
           uset = find(ni==jj);
           cur_pdown = nansum(ctx_up_isskipped{used_pdowns(ii)}(uset)==1);
           cur_npdown = nansum(ctx_up_isskipped{used_pdowns(ii)}(uset)==0);           
%           cur_tot = sum(~isnan(ctx_up_isskipped{used_pdowns(ii)}(ni==jj)));
          ctxmua_cond_pdown(ii,jj) = cur_pdown/(cur_pdown+cur_npdown);
         
          uset = uset(ctx_up_isskipped{used_pdowns(ii)}(uset)==0);
          ctxmua_cond_uplag(ii,jj) = nanmean(corresp_mp_uplag{used_pdowns(ii)}(uset));
       end
    end
    time_since_last_mp_down{used_pdowns(ii)}(time_since_last_mp_down{used_pdowns(ii)} > 10) = 10;
    
    tsu_cur_bin_edges = prctile(time_since_last_mp_down{used_pdowns(ii)},pbins);
    [nb,ni] = histc(time_since_last_mp_down{used_pdowns(ii)},tsu_cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
          cur_npdown = nansum(ctx_up_isskipped{used_pdowns(ii)}((ni==jj)==1));
          cur_tot = sum(~isnan(ctx_up_isskipped{used_pdowns(ii)}(ni==jj)));
          ctxtsu_cond_pdown(ii,jj) = cur_npdown/cur_tot;
          ctxtsu_cond_uplag(ii,jj) = nanmean(corresp_mp_uplag{used_pdowns(ii)}(ni==jj));
          ctxtsu_cond_upamp(ii,jj) = nanmean(ctx_up_amps{used_pdowns(ii)}(ni==jj));
          ctxtsu_cond_upamphf(ii,jj) = nanmean(ctx_up_amps_hf{used_pdowns(ii)}(ni==jj));
          ctxtsu_cond_upampmua(ii,jj) = nanmean(ctx_up_amps_mua{used_pdowns(ii)}(ni==jj));
       end
    end
    
    prvdowndur_cur_bin_edges = prctile(prev_ctx_downdur{used_pdowns(ii)},pbins);
    [nb,ni] = histc(prev_ctx_downdur{used_pdowns(ii)},tsu_cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
          cur_npdown = nansum(ctx_up_isskipped{used_pdowns(ii)}((ni==jj)==1));
          cur_tot = sum(~isnan(ctx_up_isskipped{used_pdowns(ii)}(ni==jj)));
          ctxpdd_cond_pdown(ii,jj) = cur_npdown/cur_tot;
          ctxpdd_cond_uplag(ii,jj) = nanmean(corresp_mp_uplag{used_pdowns(ii)}(ni==jj));
          ctxpdd_cond_upamp(ii,jj) = nanmean(ctx_up_amps{used_pdowns(ii)}(ni==jj));
          ctxpdd_cond_upamphf(ii,jj) = nanmean(ctx_up_amps_hf{used_pdowns(ii)}(ni==jj));
          ctxpdd_cond_upampmua(ii,jj) = nanmean(ctx_up_amps_mua{used_pdowns(ii)}(ni==jj));
       end
    end
end

%%
ctxhf_cond_uplag = nan(length(l3mec),nbins);
ctxhf_cond_downlag = nan(length(l3mec),nbins);
for ii = 1:length(l3mec)
    cur_bin_edges = prctile(ctx_up_amps_hf{l3mec(ii)},pbins);
    [nb,ni] = histc(ctx_up_amps_hf{l3mec(ii)},cur_bin_edges);
    for jj = 1:nbins
             uset = find(ni==jj);
     if nb(jj) > 0          
          uset = uset(ctx_up_isskipped{l3mec(ii)}(uset)==0);
          ctxhf_cond_uplag(ii,jj) = nanmean(corresp_mp_uplag{l3mec(ii)}(uset));
          uset = uset(ctx_up_isskipped{l3mec(ii)}(uset)==0);
          ctxhf_cond_downlag(ii,jj) = nanmean(corresp_mp_downlag{l3mec(ii)}(uset));
       end
    end
end
%%
cd C:\persDowns_paper\Figs\

pbin_cents = 0.5*pbins(1:end-1) + 0.5*pbins(2:end);
figure;
hold on
h2 = shadedErrorBar(pbin_cents,nanmean(ctx_cond_pdown),nanstd(ctx_cond_pdown)/sqrt(length(used_pdowns)),{'color','k'});
h3 = shadedErrorBar(pbin_cents,nanmean(ctxmua_cond_pdown),nanstd(ctxmua_cond_pdown)/sqrt(length(used_pdowns)),{'color','g'});
h1 = shadedErrorBar(pbin_cents,nanmean(ctxhf_cond_pdown),nanstd(ctxhf_cond_pdown)/sqrt(length(used_pdowns)),{'color','r'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'HF power','LF amplitude','MUA rate'});
xlabel('Ctx up-state amplitude (Percentile)','fontsize',12);
ylabel('Probability persistent down','fontsize',12);
box off
set(gca,'fontsize',10,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);
fname = 'pers_downs_vs_ctx';

% figure;
% hold on
% h2 = shadedErrorBar(pbin_cents,nanmean(ctx_cond_uplag),nanstd(ctx_cond_uplag)/sqrt(length(used_pdowns)),{'color','k'});
% h3 = shadedErrorBar(pbin_cents,nanmean(ctxmua_cond_uplag),nanstd(ctxmua_cond_uplag)/sqrt(length(used_pdowns)),{'color','g'});
% h1 = shadedErrorBar(pbin_cents,nanmean(ctxhf_cond_uplag),nanstd(ctxhf_cond_uplag)/sqrt(length(used_pdowns)),{'color','r'});
% legend([h1.mainLine h2.mainLine h3.mainLine],{'HF power','LF amplitude','MUA rate'});
% xlabel('Ctx up-state amplitude (Percentile)','fontsize',12);
% ylabel('Up-transition lag (s)','fontsize',12);
% box off
% set(gca,'fontsize',10,'fontname','arial');
% fillPage(gcf,'papersize',[5 5]);
% fname = 'Uplag_vs_ctxamp';

figure;
hold on
h1 = shadedErrorBar(pbin_cents,nanmean(ctxhf_cond_uplag),nanstd(ctxhf_cond_uplag)/sqrt(length(l3mec)),{'color','k'});
h2 = shadedErrorBar(pbin_cents,nanmean(ctxhf_cond_downlag),nanstd(ctxhf_cond_downlag)/sqrt(length(l3mec)),{'color','r'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'HF power','LF amplitude','MUA rate'});
xlabel('Ctx up-state amplitude (Percentile)','fontsize',12);
ylabel('Up-transition lag (s)','fontsize',12);
box off
set(gca,'fontsize',10,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);
fname = 'Uplag_vs_ctxamp';


tsu_bincents = 0.5*tsu_cur_bin_edges(1:end-1) + 0.5*tsu_cur_bin_edges(2:end);
pdd_bincents = 0.5*prvdowndur_cur_bin_edges(1:end-1) + 0.5*prvdowndur_cur_bin_edges(2:end);

figure; hold on
h1=shadedErrorBar(pbin_cents,nanmean(ctxpdd_cond_pdown),nanstd(ctxpdd_cond_pdown)/sqrt(length(used_pdowns)),{'color','b'});
h2=shadedErrorBar(pbin_cents,nanmean(ctxtsu_cond_pdown),nanstd(ctxtsu_cond_pdown)/sqrt(length(used_pdowns)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Ctx state','MP state'});
xlabel('Time since down-transition','fontsize',12);
ylabel('Probability persistent down','fontsize',12);
box off
set(gca,'fontsize',10,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);
fname = 'Pdown_vs_tsdown';
% figure; hold on
% h1=shadedErrorBar(tsu_bincents,nanmean(ctxtsu_cond_upamp),nanstd(ctxtsu_cond_upamp)/sqrt(length(used_pdowns)),{'color','b'});
% h2=shadedErrorBar(tsu_bincents,nanmean(ctxtsu_cond_upamphf),nanstd(ctxtsu_cond_upamphf)/sqrt(length(used_pdowns)),{'color','r'});
% legend([h1.mainLine h2.mainLine],{'Up state duration','Time since MP down-transition'});
% xlabel('Time since MP down-transition','fontsize',12);
% ylabel('Probability persistent down','fontsize',12);
% box off
% set(gca,'fontsize',10,'fontname','arial');
% fillPage(gcf,'papersize',[5 5]);

% figure; hold on
% h1=shadedErrorBar(pbin_cents,nanmean(ctxpdd_cond_upamphf),nanstd(ctxpdd_cond_upamphf)/sqrt(length(used_pdowns)),{'color','b'});
% % h2=shadedErrorBar(pbin_cents,nanmean(ctxpdd_cond_upamp),nanstd(ctxpdd_cond_upamp)/sqrt(length(used_pdowns)),{'color','r'});
% h3=shadedErrorBar(pbin_cents,nanmean(ctxpdd_cond_upampmua),nanstd(ctxpdd_cond_upampmua)/sqrt(length(used_pdowns)),{'color','g'});
% xlabel('Previous ctx down-state duration (Percentile)','fontsize',12);
% ylabel('Average ctx up-state amplitude (z)','fontsize',12);
% legend([h1.mainLine h3.mainLine],{'HF power','MUA rate'});
% box off
% set(gca,'fontsize',10,'fontname','arial');
% fillPage(gcf,'papersize',[5 5]);

%%
min_pup = 0.05;
used_pups = l3mec(fract_rt2_ups(l3mec) >= min_pup);
% used_pups = used_pups(avg_ctx_mua_rate(used_pups) > 0.1);

nbins = 20;
pbins = linspace(0,100,nbins+1);
% ctxhf_cond_mp_amp = nan(length(used_pups),nbins);
ctxhf_cond_pup = nan(length(used_pups),nbins);
ctx_cond_pup = nan(length(used_pups),nbins);
ctxdur_cond_pup = nan(length(used_pups),nbins);
ctxmua_cond_pup = nan(length(used_pups),nbins);
ctxtsu_cond_pup = nan(length(used_pups),nbins);
for ii = 1:length(used_pups)
%     cur_bin_edges = prctile(ctx_up_amps_hf{used_pups(ii)},pbins);
%     [nb,ni] = histc(ctx_up_amps_hf{used_pups(ii)},cur_bin_edges);
    cur_bin_edges = prctile(ctx_down_amps_hf{used_pups(ii)},pbins);
    [nb,ni] = histc(ctx_down_amps_hf{used_pups(ii)},cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
          cur_npup = nansum(ctx_up_ispers{used_pups(ii)}(ni==jj)==1);
          cur_tot = sum(~isnan(ctx_up_ispers{used_pups(ii)}(ni==jj)));
          ctxhf_cond_pup(ii,jj) = cur_npup/cur_tot;
       end
    end    
        cur_npdown = nansum(ctx_up_ispers{used_pups(ii)}==1);
    cur_tot = sum(~isnan(ctx_up_ispers{used_pups(ii)}));
    test(ii) = cur_npdown/cur_tot;

%     cur_bin_edges = prctile(ctx_up_amps{used_pups(ii)},pbins);
%     [nb,ni] = histc(ctx_up_amps{used_pups(ii)},cur_bin_edges);
    cur_bin_edges = prctile(-ctx_down_amps{used_pups(ii)},pbins);
    [nb,ni] = histc(-ctx_down_amps{used_pups(ii)},cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
          cur_npup = nansum(ctx_up_ispers{used_pups(ii)}(ni==jj)==1);
          cur_tot = sum(~isnan(ctx_up_ispers{used_pups(ii)}(ni==jj)));
          ctx_cond_pup(ii,jj) = cur_npup/cur_tot;
       end
    end
    cur_bin_edges = prctile(ctx_up_durs{used_pups(ii)},pbins);
    [nb,ni] = histc(ctx_up_durs{used_pups(ii)},cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
          cur_npup = nansum(ctx_up_ispers{used_pups(ii)}(ni==jj)==1);
          cur_tot = sum(~isnan(ctx_up_ispers{used_pups(ii)}(ni==jj)));
          ctxdur_cond_pup(ii,jj) = cur_npup/cur_tot;
       end
    end
%     cur_bin_edges = prctile(ctx_up_amps_mua{used_pups(ii)},pbins);
%     [nb,ni] = histc(ctx_up_amps_mua{used_pups(ii)},cur_bin_edges);
    cur_bin_edges = prctile(ctx_down_amps_mua{used_pups(ii)},pbins);
    [nb,ni] = histc(ctx_down_amps_mua{used_pups(ii)},cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
          cur_npup = nansum(ctx_up_ispers{used_pups(ii)}(ni==jj)==1);
          cur_tot = sum(~isnan(ctx_up_ispers{used_pups(ii)}(ni==jj)));
          ctxmua_cond_pup(ii,jj) = cur_npup/cur_tot;
       end
    end
    cur_bin_edges = prctile(time_since_last_mp_up{used_pups(ii)},pbins);
    [nb,ni] = histc(time_since_last_mp_up{used_pups(ii)},cur_bin_edges);
    for jj = 1:nbins
       if nb(jj) > 0
          cur_npup = nansum(ctx_up_ispers{used_pups(ii)}((ni==jj)==1));
          cur_tot = sum(~isnan(ctx_up_ispers{used_pups(ii)}(ni==jj)));
          ctxtsu_cond_pup(ii,jj) = cur_npup/cur_tot;
       end
    end
end

%%
pbin_cents = 0.5*pbins(1:end-1) + 0.5*pbins(2:end);
figure;
hold on
h1=shadedErrorBar(pbin_cents,nanmean(ctxhf_cond_pup),nanstd(ctxhf_cond_pup)/sqrt(length(used_pups)),{'color','r'});
h2=shadedErrorBar(pbin_cents,nanmean(ctx_cond_pup),nanstd(ctx_cond_pup)/sqrt(length(used_pups)),{'color','k'});
h3=shadedErrorBar(pbin_cents,nanmean(ctxmua_cond_pup),nanstd(ctxmua_cond_pup)/sqrt(length(used_pups)),{'color','g'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'HF','LF','MUA'});
xlabel('Cortical down state amplitude (prctile)','fontsize',12);
ylabel('Probability persistent up','fontsize',12);
box off
set(gca,'fontsize',10,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);
fname = 'Pup_vs_ctxamp';

figure; hold on
h1=shadedErrorBar(pbin_cents,nanmean(ctxdur_cond_pup),nanstd(ctxdur_cond_pup)/sqrt(length(used_pups)),{'color','b'});
h2=shadedErrorBar(pbin_cents,nanmean(ctxtsu_cond_pup),nanstd(ctxtsu_cond_pup)/sqrt(length(used_pups)),{'color','r'});
legend([h1.mainLine h2.mainLine ],{'Ctx','MP'});
xlabel('Time since up-transition','fontsize',12);
ylabel('Probability persistent up','fontsize',12);
box off
set(gca,'fontsize',10,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);
fname = 'Pup_vs_tsup';


%%
% close all
% min_pup = 0.05;
% used_pups = l3mec(fract_rt2_ups(l3mec) >= min_pup);
% % used_pups = used_pups(avg_ctx_mua_rate(used_pups) > 0.1);
% 
% nbins = 20;
% pbins = linspace(0,100,nbins+1);
% ctxhf_cond_pup = nan(length(used_pups),nbins);
% ctx_cond_pup = nan(length(used_pups),nbins);
% ctxdur_cond_pup = nan(length(used_pups),nbins);
% ctxmua_cond_pup = nan(length(used_pups),nbins);
% ctxtsu_cond_pup = nan(length(used_pups),nbins);
% ctxpdd_cond_pup = nan(length(used_pups),nbins);
% for ii = 1:length(used_pups)
%     cur_bin_edges = prctile(ctx_up_amps_hf{used_pups(ii)},pbins);
%     [nb,ni] = histc(ctx_up_amps_hf{used_pups(ii)},cur_bin_edges);
%     for jj = 1:nbins
%        if nb(jj) > 0
%            uset = find(ni==jj);
%            cur_pup = nansum(ctx_up_ispers{used_pups(ii)}(uset)==1);
%            cur_npup = nansum(ctx_up_ispers{used_pups(ii)}(uset)==0);           
%           ctxhf_cond_pup(ii,jj) = cur_pup/(cur_pup+cur_npup);
%        end
%     end
%     cur_npup = nansum(ctx_up_ispers{used_pups(ii)}==1);
%     cur_tot = sum(~isnan(ctx_up_ispers{used_pups(ii)}));
%     
%     cur_bin_edges = prctile(ctx_up_amps{used_pups(ii)},pbins);
%     [nb,ni] = histc(ctx_up_amps{used_pups(ii)},cur_bin_edges);
%     for jj = 1:nbins
%         if nb(jj) > 0
%             uset = find(ni==jj);
%             cur_pup = nansum(ctx_up_ispers{used_pups(ii)}(uset)==1);
%             cur_npup = nansum(ctx_up_ispers{used_pups(ii)}(uset)==0);
%             ctx_cond_pup(ii,jj) = cur_pup/(cur_pup+cur_npup);
%          end
%     end
%     cur_bin_edges = prctile(ctx_up_durs{used_pups(ii)},pbins);
%     [nb,ni] = histc(ctx_up_durs{used_pups(ii)},cur_bin_edges);
%     for jj = 1:nbins
%        if nb(jj) > 0
%            uset = find(ni==jj);
%            cur_pup = nansum(ctx_up_ispers{used_pups(ii)}(uset)==1);
%            cur_npup = nansum(ctx_up_ispers{used_pups(ii)}(uset)==0);           
%           ctxdur_cond_pup(ii,jj) = cur_pup/(cur_pup+cur_npup);
%        end
%     end
%     cur_bin_edges = prctile(ctx_up_amps_mua{used_pups(ii)},pbins);
%     [nb,ni] = histc(ctx_up_amps_mua{used_pups(ii)},cur_bin_edges);
%     for jj = 1:nbins
%        if nb(jj) > 0
%            uset = find(ni==jj);
%            cur_pup = nansum(ctx_up_ispers{used_pups(ii)}(uset)==1);
%            cur_npup = nansum(ctx_up_ispers{used_pups(ii)}(uset)==0);           
%           ctxmua_cond_pup(ii,jj) = cur_pup/(cur_pup+cur_npup);
%        end
%     end
%     time_since_last_mp_up{used_pups(ii)}(time_since_last_mp_up{used_pups(ii)} > 10) = 10;
%     
%     tsu_cur_bin_edges = prctile(time_since_last_mp_up{used_pups(ii)},pbins);
%     [nb,ni] = histc(time_since_last_mp_up{used_pups(ii)},tsu_cur_bin_edges);
%     for jj = 1:nbins
%        if nb(jj) > 0
%           cur_npup = nansum(ctx_up_ispers{used_pups(ii)}((ni==jj)==1));
%           cur_tot = sum(~isnan(ctx_up_ispers{used_pups(ii)}(ni==jj)));
%           ctxtsu_cond_pup(ii,jj) = cur_npup/cur_tot;
%        end
%     end
%     
%     prvupdur_cur_bin_edges = prctile(prev_ctx_downdur{used_pups(ii)},pbins);
%     [nb,ni] = histc(prev_ctx_downdur{used_pups(ii)},tsu_cur_bin_edges);
%     for jj = 1:nbins
%        if nb(jj) > 0
%           cur_npup = nansum(ctx_up_ispers{used_pups(ii)}((ni==jj)==1));
%           cur_tot = sum(~isnan(ctx_up_ispers{used_pups(ii)}(ni==jj)));
%           ctxpdd_cond_pup(ii,jj) = cur_npup/cur_tot;
%         end
%     end
% end

%%
cd C:\persdowns_paper\Figs\

pbin_cents = 0.5*pbins(1:end-1) + 0.5*pbins(2:end);
figure;
hold on
h2 = shadedErrorBar(pbin_cents,nanmean(ctx_cond_pup),nanstd(ctx_cond_pup)/sqrt(length(used_pups)),{'color','k'});
h3 = shadedErrorBar(pbin_cents,nanmean(ctxmua_cond_pup),nanstd(ctxmua_cond_pup)/sqrt(length(used_pups)),{'color','g'});
h1 = shadedErrorBar(pbin_cents,nanmean(ctxhf_cond_pup),nanstd(ctxhf_cond_pup)/sqrt(length(used_pups)),{'color','r'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'HF power','LF amplitude','MUA rate'});
xlabel('Ctx up-state amplitude (Percentile)','fontsize',12);
ylabel('Probability persistent up','fontsize',12);
box off
set(gca,'fontsize',10,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);


tsu_bincents = 0.5*tsu_cur_bin_edges(1:end-1) + 0.5*tsu_cur_bin_edges(2:end);
pdd_bincents = 0.5*prvupdur_cur_bin_edges(1:end-1) + 0.5*prvupdur_cur_bin_edges(2:end);

figure; hold on
h1=shadedErrorBar(pbin_cents,nanmean(ctxpdd_cond_pup),nanstd(ctxpdd_cond_pup)/sqrt(length(used_pups)),{'color','b'});
h2=shadedErrorBar(pbin_cents,nanmean(ctxtsu_cond_pup),nanstd(ctxtsu_cond_pup)/sqrt(length(used_pups)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Ctx state','MP state'});
xlabel('Time since up-transition','fontsize',12);
ylabel('Probability persistent up','fontsize',12);
box off
set(gca,'fontsize',10,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);



