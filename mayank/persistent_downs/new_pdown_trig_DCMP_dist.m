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
d=3
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
    mp_state_vec = nan(length(wcv_lf),1);
    lfp_state_number = nan(length(wcv_lf),1);
    for i = 1:length(up_trans_inds)-1
        mp_state_vec(up_trans_inds(i):down_trans_inds(i)) = 1;
        mp_state_vec(down_trans_inds(i):up_trans_inds(i+1)) = 0;
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
    dc_uds_amp(d) = dc_upstate_amp(d) - dc_downstate_amp(d);
    
    dc_interp_data_reld = dc_interp_data - dc_downstate_amp(d);
    dc_interp_data_relu = dc_interp_data - dc_upstate_amp(d);
    dc_interp_data_normd = dc_interp_data_reld/dc_uds_amp(d);
    dc_interp_data_normu = dc_interp_data_relu/dc_uds_amp(d);
    
        
    %%
     dtrans_buffer = round(Fsd*0.1);
    down_trans_inds8 = down_trans_inds8 + dtrans_buffer;
    down_trans_inds8(down_trans_inds8 > length(lfp_state_number)) = length(lfp_state_number);
    down_trans_inds = down_trans_inds + dtrans_buffer;
    down_trans_inds(down_trans_inds > length(lfp_state_number)) = length(lfp_state_number);

    mp_state_number = nan(length(wcv_lf),1);
    lfp_state_number = nan(length(wcv_lf),1);
    for i = 1:length(up_trans_inds)-1
        mp_state_number(up_trans_inds(i):down_trans_inds(i)) = 2*(i-1)+1;
        mp_state_number(down_trans_inds(i):up_trans_inds(i+1)) = 2*(i-1) + 2;
        mp_state_vec(up_trans_inds(i):down_trans_inds(i)) = 1;
        mp_state_vec(down_trans_inds(i):up_trans_inds(i+1)) = 0;
    end
    
    for i = 1:length(up_trans_inds8)-1
        lfp_state_number(up_trans_inds8(i):down_trans_inds8(i)) = 2*(i-1)+1;
        lfp_state_number(down_trans_inds8(i):up_trans_inds8(i+1)) = 2*(i-1) + 2;
    end
    mp_state_number(isnan(lf8_period_vec)) = nan;
    lfp_state_number(isnan(lf8_period_vec)) = nan;
    mp_state_vec(isnan(lf8_period_vec)) = nan;

        ctx_up_mp_down = up_trans_inds8(mp_state_vec(up_trans_inds8) == 0);
    ctx_down_mp_up = down_trans_inds8(mp_state_vec(down_trans_inds8) == 1);
    mp_down_nonskipped = ctx_up_mp_down(ismember(ctx_up_mp_down,up_trans_inds8(non_dskipped_lfp_trans)));

    %%
    xi_drel = linspace(-10,40,200);
    xi_urel = linspace(-40,20,200);
    xi_drel_norm = linspace(-0.4,1.5,200);
    xi_urel_norm = linspace(-1.5,0.75,200);
    backlag = round(Fsd*0.5);
    forwardlag = round(Fsd*1.25);
%     [ev_avg,lags,ev_std,n_events,utrig_mat] = get_event_trig_avg(dc_interp_data_reld(:),ctx_up_mp_down,backlag,forwardlag,1,mp_state_number,-1);
    [ev_avg,lags,ev_std,n_events,utrig_matN] = get_event_trig_avg(dc_interp_data_normd(:),ctx_up_mp_down,backlag,forwardlag,1,mp_state_number,-1);
    [ev_avg,lags,ev_std,n_events,utrig_matN_npdown] = get_event_trig_avg(dc_interp_data_normd(:),mp_down_nonskipped,backlag,forwardlag,1,mp_state_number,-1);
%     [ev_avg,lags,ev_std,n_events,mutrig_matN] = get_event_trig_avg(dc_interp_data_normd(:),up_trans_inds,backlag,forwardlag,1,mp_state_number,-2);
    if length(dskipped_lfp_trans) > 5
%         [ev_avg,lags,ev_std,n_events,utrig_mat_pdown] = get_event_trig_avg(dc_interp_data_reld(:),up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag,1,mp_state_number,-1);
        [ev_avg,lags,ev_std,n_events,utrig_mat_pdownN] = get_event_trig_avg(dc_interp_data_normd(:),up_trans_inds8(dskipped_lfp_trans),backlag,forwardlag,1,mp_state_number,-1);
    end
%     cur_amp_gram = nan(length(lags),length(xi_drel));
%     cur_amp_gram_pdown = nan(length(lags),length(xi_drel));
    cur_amp_gramN = nan(length(lags),length(xi_drel));
%     cur_amp_mgramN = nan(length(lags),length(xi_drel));
    cur_amp_gram_npdownN = nan(length(lags),length(xi_drel));
    cur_amp_gram_pdownN = nan(length(lags),length(xi_drel));
    for ii = 1:length(lags)
        uset = ~isnan(utrig_matN(:,ii));
        if sum(uset) > 5
%             cur_amp_gram(ii,:) = ksdensity(utrig_mat(:,ii),xi_drel);
            cur_amp_gramN(ii,:) = ksdensity(utrig_matN(:,ii),xi_drel_norm);
            cur_amp_gramN_npdown(ii,:) = ksdensity(utrig_matN_npdown(:,ii),xi_drel_norm);
        else
%             cur_amp_gram(ii,:) = nan(1,length(xi_drel));
            cur_amp_gramN(ii,:) = nan(1,length(xi_drel));
            cur_amp_gramN_npdown(ii,:) = nan(1,length(xi_drel));
        end
%         uset = ~isnan(mutrig_matN(:,ii));
%         if sum(uset) > 5
%             cur_amp_mgramN(ii,:) = ksdensity(mutrig_matN(:,ii),xi_drel_norm);
%         else
%             cur_amp_mgramN(ii,:) = nan(1,length(xi_drel_norm));
%         end
        if length(dskipped_lfp_trans) > 5
            uset = ~isnan(utrig_mat_pdownN(:,ii));
            if sum(uset) > 5
%                 cur_amp_gram_pdown(ii,:) = ksdensity(utrig_mat_pdown(uset,ii),xi_drel);
                cur_amp_gram_pdownN(ii,:) = ksdensity(utrig_mat_pdownN(uset,ii),xi_drel_norm);
            end
        end
    end
%     ctx_utrig_ampdist(d,:,:) = cur_amp_gram;
    ctx_utrig_ampdistN(d,:,:) = cur_amp_gramN;
    ctx_utrig_ampdistN_npdown(d,:,:) = cur_amp_gramN_npdown;
%     ctx_mutrig_ampdistN(d,:,:) = cur_amp_mgramN;
    if length(dskipped_lfp_trans) > 5
%         ctx_utrig_ampdist_pdown(d,:,:) = cur_amp_gram_pdown;
        ctx_utrig_ampdist_pdownN(d,:,:) = cur_amp_gram_pdownN;
    else
%         ctx_utrig_ampdist_pdown(d,:,:) = nan(length(lags),length(xi_drel));
        ctx_utrig_ampdist_pdownN(d,:,:) = nan(length(lags),length(xi_drel));
    end
    
    non_skipped_dtrans = setdiff(ctx_down_mp_up,down_trans_inds8(skipped_lfp_trans));
%     [ev_avg,lags,ev_std,n_events,dtrig_mat] = get_event_trig_avg(dc_interp_data_relu(:),ctx_down_mp_up,backlag,forwardlag,1,mp_state_number,-1);
    [ev_avg,lags,ev_std,n_events,dtrig_matN] = get_event_trig_avg(dc_interp_data_normu(:),ctx_down_mp_up,backlag,forwardlag,1,mp_state_number,-1);
    [ev_avg,lags,ev_std,n_events,dtrig_matN_npup] = get_event_trig_avg(dc_interp_data_normu(:),non_skipped_dtrans,backlag,forwardlag,1,mp_state_number,-1);
%     [ev_avg,lags,ev_std,n_events,mdtrig_matN] = get_event_trig_avg(dc_interp_data_normu(:),down_trans_inds,backlag,forwardlag,1,mp_state_number,-2);
    if length(skipped_lfp_trans) > 5
%         [ev_avg,lags,ev_std,n_events,dtrig_mat_pup] = get_event_trig_avg(dc_interp_data_relu(:),down_trans_inds8(skipped_lfp_trans),backlag,forwardlag,1,mp_state_number,-1);
        [ev_avg,lags,ev_std,n_events,dtrig_mat_pupN] = get_event_trig_avg(dc_interp_data_normu(:),down_trans_inds8(skipped_lfp_trans),backlag,forwardlag,1,mp_state_number,-1);
    end
%     cur_amp_gram = nan(length(lags),length(xi_urel));
%     cur_amp_gram_pup = nan(length(lags),length(xi_urel));
    cur_amp_gramN = nan(length(lags),length(xi_urel));
    cur_amp_gramN_npup = nan(length(lags),length(xi_urel));
%     cur_amp_mgramN = nan(length(lags),length(xi_urel));
    cur_amp_gram_pupN = nan(length(lags),length(xi_urel));
    for ii = 1:length(lags)
        uset = ~isnan(dtrig_matN(:,ii));
        if sum(uset) > 5
%             cur_amp_gram(ii,:) = ksdensity(dtrig_mat(uset,ii),xi_urel);
            cur_amp_gramN(ii,:) = ksdensity(dtrig_matN(uset,ii),xi_urel_norm);
        else
%             cur_amp_gram(ii,:) = nan(1,length(xi_urel));
            cur_amp_gramN(ii,:) = nan(1,length(xi_urel));
        end
        uset = ~isnan(dtrig_matN_npup(:,ii));
        if sum(uset) > 5
            cur_amp_gramN_npup(ii,:) = ksdensity(dtrig_matN_npup(uset,ii),xi_urel_norm);
        else
            cur_amp_gramN_npup(ii,:) = nan(1,length(xi_urel));
        end
%         uset = ~isnan(mdtrig_matN(:,ii));
%         if sum(uset) > 5
%             cur_amp_mgramN(ii,:) = ksdensity(mdtrig_matN(uset,ii),xi_urel_norm);
%         else
%             cur_amp_mgramN(ii,:) = nan(1,length(xi_urel));
%         end
        if length(skipped_lfp_trans) > 5
            uset = ~isnan(dtrig_mat_pupN(:,ii));
            if sum(uset) > 5
%                 cur_amp_gram_pup(ii,:) = ksdensity(dtrig_mat_pup(uset,ii),xi_urel);
                cur_amp_gram_pupN(ii,:) = ksdensity(dtrig_mat_pupN(uset,ii),xi_urel_norm);
            end
        end
    end
%     ctx_dtrig_ampdist(d,:,:) = cur_amp_gram;
    ctx_dtrig_ampdistN(d,:,:) = cur_amp_gramN;
    ctx_dtrig_ampdistN_npup(d,:,:) = cur_amp_gramN;
%     ctx_mdtrig_ampdistN(d,:,:) = cur_amp_mgramN;
    if length(skipped_lfp_trans) > 5
%         ctx_dtrig_ampdist_pup(d,:,:) = cur_amp_gram_pup;
        ctx_dtrig_ampdist_pupN(d,:,:) = cur_amp_gram_pupN;
    else
%         ctx_dtrig_ampdist_pup(d,:,:) = nan(length(lags),length(xi_urel));
        ctx_dtrig_ampdist_pupN(d,:,:) = nan(length(lags),length(xi_urel));
    end
    
end

%%
cd C:\WC_Germany\persistent_downs\
save new_down_trig_DCMPdist lags xi* ctx_*

%%
load ./new_down_core_analysis fract*

min_nstates = 15;
l3mec = find(strcmp(data_type(used_dirs),'L3MEC'));
l3lec = find(strcmp(data_type(used_dirs),'L3LEC'));
% pup_used = find(n_uskip > min_nstates & ismember(used_dirs,l3mec));
% pdown_used = find(n_dskip > min_nstates & ismember(used_dirs,l3mec));
pup_used = l3mec(fract_rt2_ups(l3mec) > 0.05);
pdown_used = l3mec(fract_rt2_downs(l3mec) > 0.05);

%this cell is too dominated by the up state to get a reliable reading of
%the trig avg distributions
pup_used(pup_used == 70) = [];
pdown_used(pdown_used==70) = [];

%no heka data for this cell
pup_used(pup_used == 44) = [];
pdown_used(pdown_used==44) = [];

%no heka data for this cell
pup_used(pup_used == 86) = [];
pdown_used(pdown_used==86) = [];

%%
% figure
% imagesc(lags/Fsd,xi_drel,squeeze(nanmean(sqrt(ctx_utrig_ampdist(pdown_used,:,:))))');
% set(gca,'ydir','normal');
% 
% figure
% imagesc(lags/Fsd,xi_drel,squeeze(nanmean(sqrt(ctx_utrig_ampdistN_npdown(pdown_used,:,:))))');
% set(gca,'ydir','normal');
% 
% % figure
% % imagesc(lags/Fsd,xi_drel,squeeze(nanmean(sqrt(ctx_utrig_ampdist(l3lec,:,:))))');
% % set(gca,'ydir','normal');
% 
% figure
% imagesc(lags/Fsd,xi_drel,squeeze(nanmean(sqrt(ctx_utrig_ampdist_pdown(pdown_used,:,:))))');
% set(gca,'ydir','normal');
%%
[~,ctx_utrig_maxlocs_npdown] = nanmax(ctx_utrig_ampdistN_npdown(pdown_used,:,:),[],3);
uu = find(~isnan(ctx_utrig_maxlocs_npdown));
ctx_utrig_maxlocs_npdown(uu) = xi_drel_norm(ctx_utrig_maxlocs_npdown(uu));

[~,ctx_utrig_maxlocs] = nanmax(ctx_utrig_ampdist_pdownN(pdown_used,:,:),[],3);
uu = find(~isnan(ctx_utrig_maxlocs));
ctx_utrig_maxlocs(uu) = xi_drel_norm(ctx_utrig_maxlocs(uu));
figure;hold on
shadedErrorBar(lags/Fsd,nanmedian(ctx_utrig_maxlocs_npdown),nanstd(ctx_utrig_maxlocs_npdown)/sqrt(length(pdown_used)));
shadedErrorBar(lags/Fsd,nanmedian(ctx_utrig_maxlocs),nanstd(ctx_utrig_maxlocs)/sqrt(length(pdown_used)),{'color','r'});

modefuns = nan(length(pdown_used),length(lags));
for ii = 1:length(pdown_used)
    mode_deriv = [nan diff(ctx_utrig_maxlocs_npdown(ii,:))];
    [~,temploc(ii)] = nanmax(mode_deriv);
    modefuns(ii,1:temploc(ii)) = ctx_utrig_maxlocs_npdown(ii,1:temploc(ii));
end

%%
figure
imagesc(lags/Fsd,xi_drel_norm,squeeze(nanmedian(sqrt(ctx_utrig_ampdistN_npdown(pdown_used,:,:))))');
set(gca,'ydir','normal');
xlabel('Time since ctx up-trans (s)','fontsize',12);
ylabel('MP relative amplitude','fontsize',12);
yl = ylim();
line([0 0],yl,'color','k');
hold on
shadedErrorBar(lags/Fsd,nanmedian(ctx_utrig_maxlocs_npdown),nanstd(ctx_utrig_maxlocs_npdown)/sqrt(length(pdown_used)));

% figure
% imagesc(lags/Fsd,xi_drel_norm,squeeze(nanmean(sqrt(ctx_utrig_ampdistN(l3lec,:,:))))');
% set(gca,'ydir','normal');

figure
imagesc(lags/Fsd,xi_drel_norm,squeeze(nanmedian(sqrt(ctx_utrig_ampdist_pdownN(pdown_used,:,:))))');
set(gca,'ydir','normal');
xlabel('Time since ctx up-trans (s)','fontsize',12);
ylabel('MP relative amplitude','fontsize',12);
yl = ylim();
line([0 0],yl,'color','k');
hold on
shadedErrorBar(lags/Fsd,nanmedian(ctx_utrig_maxlocs),nanstd(ctx_utrig_maxlocs)/sqrt(length(pdown_used)));

% figure
% imagesc(lags/Fsd,xi_drel_norm,squeeze(nanmean(sqrt(ctx_mutrig_ampdistN(pdown_used,:,:))))');
% set(gca,'ydir','normal');


%%
load new_down_core_analysis mp_uplags nrt2_downs
for ii = 1:length(pdown_used)
%     ii = 16
ii
    subplot(3,1,1);
    imagesc(lags/Fsd,xi_drel_norm,squeeze(sqrt(ctx_utrig_ampdistN_npdown(pdown_used(ii),:,:)))');
    set(gca,'ydir','normal');
    xlabel('Time since ctx up-trans (s)','fontsize',12);
    ylabel('MP relative amplitude','fontsize',12);
    hold on
    plot(lags/Fsd,ctx_utrig_maxlocs_npdown(ii,:),'w','linewidth',2);
    yl = ylim();
    line([0 0],yl,'color','k');
    subplot(3,1,2);
    imagesc(lags/Fsd,xi_drel_norm,squeeze(sqrt(ctx_utrig_ampdist_pdownN(pdown_used(ii),:,:)))');
    set(gca,'ydir','normal');
    xlabel('Time since ctx up-trans (s)','fontsize',12);
    ylabel('MP relative amplitude','fontsize',12);
    hold on
    plot(lags/Fsd,ctx_utrig_maxlocs(ii,:),'w','linewidth',2);
    yl = ylim();
    line([0 0],yl,'color','k');
    subplot(3,1,3)
    hist(mp_uplags{pdown_used(ii)}(nrt2_downs{pdown_used(ii)})/Fsd,50);
    xlim(lags([1 end])/Fsd);
    pause
    clf
end


%%
cur_avg = squeeze(nanmean(sqrt(ctx_utrig_ampdistN(pdown_used,:,:))))';
cur_avg = bsxfun(@rdivide,cur_avg,max(cur_avg));
figure
imagesc(lags/Fsd,xi_drel_norm,cur_avg);
set(gca,'ydir','normal');

cur_avg = squeeze(nanmean(sqrt(ctx_utrig_ampdistN(l3lec,:,:))))';
cur_avg = bsxfun(@rdivide,cur_avg,max(cur_avg));
figure
imagesc(lags/Fsd,xi_drel_norm,cur_avg);
set(gca,'ydir','normal');

cur_avg = squeeze(nanmean(sqrt(ctx_utrig_ampdist_pdownN(pdown_used,:,:))))';
cur_avg = bsxfun(@rdivide,cur_avg,max(cur_avg));
figure
imagesc(lags/Fsd,xi_drel_norm,cur_avg);
set(gca,'ydir','normal');


%%
figure
imagesc(lags,xi_urel,squeeze(nanmean(sqrt(ctx_dtrig_ampdist(pup_used,:,:))))');
set(gca,'ydir','normal');

figure
imagesc(lags,xi_urel,squeeze(nanmean(sqrt(ctx_dtrig_ampdist(l3lec,:,:))))');
set(gca,'ydir','normal');

figure
imagesc(lags,xi_urel,squeeze(nanmean(sqrt(ctx_dtrig_ampdist_pup(pup_used,:,:))))');
set(gca,'ydir','normal');

%%
figure
imagesc(lags/Fsd,xi_urel_norm,squeeze(nanmean(sqrt(ctx_dtrig_ampdistN(pup_used,:,:))))');
set(gca,'ydir','normal');
xlabel('Time since ctx down-trans (s)','fontsize',12);
ylabel('MP relative amplitude','fontsize',12);
yl = ylim();
line([0 0],yl,'color','k');

figure
imagesc(lags/Fsd,xi_urel_norm,squeeze(nanmean(sqrt(ctx_dtrig_ampdistN(l3lec,:,:))))');
set(gca,'ydir','normal');

figure
imagesc(lags/Fsd,xi_urel_norm,squeeze(nanmean(sqrt(ctx_dtrig_ampdist_pupN(pup_used,:,:))))');
set(gca,'ydir','normal');
xlabel('Time since ctx down-trans (s)','fontsize',12);
ylabel('MP relative amplitude','fontsize',12);
yl = ylim();
line([0 0],yl,'color','k');

figure
imagesc(lags/Fsd,xi_urel_norm,squeeze(nanmean(sqrt(ctx_mdtrig_ampdistN(pup_used,:,:))))');
set(gca,'ydir','normal');

%%
cur_avg = squeeze(nanmean(sqrt(ctx_dtrig_ampdistN(pup_used,:,:))))';
cur_avg = bsxfun(@rdivide,cur_avg,max(cur_avg));
figure
imagesc(lags,xi_urel_norm,cur_avg);
set(gca,'ydir','normal');

cur_avg = squeeze(nanmean(sqrt(ctx_dtrig_ampdist_pupN(pup_used,:,:))))';
cur_avg = bsxfun(@rdivide,cur_avg,max(cur_avg));
figure
imagesc(lags,xi_urel_norm,cur_avg);
set(gca,'ydir','normal');
