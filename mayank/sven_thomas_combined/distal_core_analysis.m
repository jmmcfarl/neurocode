
clear all
close all

cur_dir_letter = 'C';
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
cd C:\WC_Germany\sven_thomas_combined\

load .\distal_dir
distal_dir = distal_dir(distal_usable);
mctx_lfp = ctx_lfp(distal_usable);
mhpc_lfp = hpc_lfp(distal_usable);
mhpc_mua = hpc_mua(distal_usable);

load ./distal_lec_dir.mat
distal_dir = [distal_dir distal_lec(usable_distal_lec)];
ctx_lfp = [mctx_lfp ctx_lfp(usable_distal_lec)];
hpc_lfp = [mhpc_lfp hpc_lfp(usable_distal_lec)];
hpc_mua = [mhpc_mua hpc_mua(usable_distal_lec)];
distal_mec = (1:length(distal_usable));
distal_lec = (length(distal_usable)+1):(length(distal_usable)+length(usable_distal_lec));

%%
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;

maxlag = round(4*Fsd);
backlag = 4*Fsd;
forwardlag = 4*Fsd;

params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 25;

%robust persistence requires MP skips states longer than this
thresh_lf8_updur = 0.5;
thresh_lf8_downdur = 0.5;

%for state duration distribution calculations
up_range = [0.1 15];
down_range = [0.1 15];
% up_range = [0.1 3];
% down_range = [0.1 3];
numBins = 60;
lin_dur_grid = linspace(up_range(1),up_range(2),numBins+1);
log_dur_grid = logspace(log10(up_range(1)),log10(up_range(2)),numBins+1);

%%
% uset = [1 4 5 9 10];
% distal_dir = distal_dir(uset);
% hpc_mua = hpc_mua(uset);
% hpc_lfp = hpc_lfp(uset);
% ctx_lfp = ctx_lfp(uset);
for d = 1:length(distal_dir)
cur_dir = distal_dir{d};
    cd(cur_dir)
    pwd
    
    load ./used_data lf8 lf7 lf6 lf5 wcv_minus_spike
    if ctx_lfp(d) == 7
        lf8 = lf7;
    elseif ctx_lfp(d) == 6
        lf8 = lf6;
    elseif ctx_lfp(d) == 5
        lf8 = lf5;
    end
    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
%     [lf82_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[0.05 2]);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
%      plot(t_axis,lf8_lf,'r',t_axis,wcv_lf+3)
% hold on
% plot(t_axis,lf82_lf,'k')
%     if hpc_lfp(d) == 3
%         load ./used_data lf3 lf5
%         if ismember(d,old_data_inds)
%             lf3 = lf3 + lf5; %redefine LF3 wrt gnd
%         end
%         hpc_lf = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf]);
%         hpc_hf = get_hf_features(lf3,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
%     else
%         load ./used_data lf2
%         hpc_lf = get_lf_features(lf2,raw_Fs,Fsd,[lcf hcf]);
%         hpc_hf = get_hf_features(lf2,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
%     end
    
%     if ismember(d,new_data_inds)
%         load ./pa_hsmm_state_seq8_newmec
%         load ./pa_hsmm_state_seq_newmec
%     else
%         load ./pa_hsmm_state_seq_new2
%         load ./pa_hsmm_state_seq8_new2
%     end

%     load ./pa_hsmm_state_seq_combined_fin.mat
%     load ./pa_hsmm_state_seq7_combined_fin.mat
    load ./pa_hsmm_state_seq_combined_fin_nd.mat
    load ./pa_hsmm_state_seq7_combined_fin_nd.mat
    hsmm_bbstate_seq8 = hsmm_bbstate_seq7;
    hsmm8 = hsmm7;
    
%         if you want to flip the roles of the MP and LFP
%         temp = hsmm_bbstate_seq8;
%         hsmm_bbstate_seq8 = hsmm_bbstate_seq;
%         hsmm_bbstate_seq = temp;

    [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,hsmm_bbstate_seq);
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq8);
    fract_uds_dur(d) = dur_uds/t_axis(end);
    
    %%
%     load ./combined_lf7_period_data_fin_nd.mat
%     load ./combined_lf7_period_data_fin.mat
%     load ./combined_lf8_period_data
    load ./combined_lf7_period_data_fin_nd %to reverse roles of MP and LFP
%     load ./combined_mp_period_data_fin_nd %to reverse roles of MP and LFP
%    load ./combined_mp_period_data_fin %to reverse roles of MP and LFP

    lf8_period_vec = nan(size(wcv_lf));
    for i = 1:size(new_seg_inds,1)
        lf8_period_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = lf8_period_f{i};
    end
    
    %% compute state duration distributions
    [mp_state_durations{d}] = compute_state_durations_seg(hsmm_bbstate_seq,Fsd);
    [lfp_state_durations{d}] = compute_state_durations_seg(hsmm_bbstate_seq8,Fsd);

    [up_amps,down_amps] = get_state_amplitudes(lf8_lf,up_trans_inds8,down_trans_inds8);
    
    up_dur_hist(d,:) = histc(mp_state_durations{d}{2},lin_dur_grid);
    up_dur_hist8(d,:) = histc(lfp_state_durations{d}{2},lin_dur_grid);    
    down_dur_hist(d,:) = histc(mp_state_durations{d}{1},lin_dur_grid);
    down_dur_hist8(d,:) = histc(lfp_state_durations{d}{1},lin_dur_grid);
    up_dur_loghist(d,:) = histc(mp_state_durations{d}{2},log_dur_grid);
    up_dur_loghist8(d,:) = histc(lfp_state_durations{d}{2},log_dur_grid); 
    down_dur_loghist(d,:) = histc(mp_state_durations{d}{1},log_dur_grid);
    down_dur_loghist8(d,:) = histc(lfp_state_durations{d}{1},log_dur_grid);

    mp_updurs = mp_state_durations{d}{2};
    mp_downdurs = mp_state_durations{d}{1};
    lf8_updurs= lfp_state_durations{d}{2};
    lf8_downdurs = lfp_state_durations{d}{1};

    mp_dutycycs = mp_updurs./(mp_updurs+mp_downdurs);
    lf8_dutycycs = lf8_updurs./(lf8_updurs+lf8_downdurs);
    median_dc(d) = nanmedian(mp_dutycycs);
    median_dc8(d) = nanmedian(lf8_dutycycs);
    mean_dc(d) = nanmean(mp_dutycycs);
    mean_dc8(d) = nanmean(lf8_dutycycs);
    
    mean_updur8(d) = nanmean(lfp_state_durations{d}{2});
    mean_downdur8(d) = nanmean(lfp_state_durations{d}{1});
    median_updur8(d) = nanmedian(lfp_state_durations{d}{2});
    median_downdur8(d) = nanmedian(lfp_state_durations{d}{1});
    mean_updur(d) = nanmean(mp_state_durations{d}{2});
    mean_downdur(d) = nanmean(mp_state_durations{d}{1});
    median_updur(d) = nanmedian(mp_state_durations{d}{2});
    median_downdur(d) = nanmedian(mp_state_durations{d}{1});
    max_updur8(d) = nanmax(lfp_state_durations{d}{2});
    max_updur(d) = nanmax(mp_state_durations{d}{2});
    max_downdur8(d) = nanmax(lfp_state_durations{d}{1});
    max_downdur(d) = nanmax(mp_state_durations{d}{1});
   
    
    
    %% compute transition-triggered averages
    lf8_utrig_lf8(d,:) = get_event_trig_avg(lf8_lf,up_trans_inds8,forwardlag,backlag);
    lf8_dtrig_lf8(d,:) = get_event_trig_avg(lf8_lf,down_trans_inds8,forwardlag,backlag);
    lf8_utrig_wcv(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds8,forwardlag,backlag);
    lf8_dtrig_wcv(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds8,forwardlag,backlag);
%     lf8_utrig_hpc(d,:) = get_event_trig_avg(hpc_lf,up_trans_inds8,forwardlag,backlag);
%     lf8_dtrig_hpc(d,:) = get_event_trig_avg(hpc_lf,down_trans_inds8,forwardlag,backlag);
%     lf8_utrig_hpchf(d,:) = get_event_trig_avg(hpc_hf,up_trans_inds8,forwardlag,backlag);
%     lf8_dtrig_hpchf(d,:) = get_event_trig_avg(hpc_hf,down_trans_inds8,forwardlag,backlag);

    wcv_utrig_lf8(d,:) = get_event_trig_avg(lf8_lf,up_trans_inds,forwardlag,backlag);
    wcv_dtrig_lf8(d,:) = get_event_trig_avg(lf8_lf,down_trans_inds,forwardlag,backlag);
    wcv_utrig_wcv(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds,forwardlag,backlag);
    wcv_dtrig_wcv(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds,forwardlag,backlag);
%     wcv_utrig_hpc(d,:) = get_event_trig_avg(hpc_lf,up_trans_inds,forwardlag,backlag);
%     wcv_dtrig_hpc(d,:) = get_event_trig_avg(hpc_lf,down_trans_inds,forwardlag,backlag);
%     wcv_utrig_hpchf(d,:) = get_event_trig_avg(hpc_hf,up_trans_inds,forwardlag,backlag);
%     wcv_dtrig_hpchf(d,:) = get_event_trig_avg(hpc_hf,down_trans_inds,forwardlag,backlag);
    
    lags = (-backlag:forwardlag)/Fsd;
        
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
    
    %compute transition lags relative to ctx state durations
    mp_reldownlags{d} = nan(size(up_trans_inds));
    mp_reluplags{d} = nan(size(up_trans_inds));
    for i = 1:length(up_trans_inds)
        if ~isnan(corresp_lf8_downinds{d}(i))
            mp_reldownlags{d}(i) = mp_downlags{d}(i)/lf8_downdurs(corresp_lf8_downinds{d}(i))/Fsd;
        end
        if ~isnan(corresp_lf8_upinds{d}(i))
            mp_reluplags{d}(i) = mp_uplags{d}(i)/lf8_updurs(corresp_lf8_upinds{d}(i))/Fsd;
        end
    end

    %% compute persistence
    [mp_upskipped{d},mp_downskipped{d}] = greedy_find_skipped_ncx_states(...
        corresp_lf8_upinds{d},corresp_lf8_downinds{d},lfp_state_durations{d}{2},lfp_state_durations{d}{1},thresh_lf8_downdur,thresh_lf8_updur);
            
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
    
    mp_rel_updurs{d} = nan(size(up_trans_inds));
    mp_rel_downdurs{d} = nan(size(up_trans_inds));
    mp_rel_updurs{d}(robust_non_skipped_mp_ups{d}) = mp_state_durations{d}{2}(robust_non_skipped_mp_ups{d}) - ...
        lfp_state_durations{d}{2}(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups{d}));
    mp_rel_downdurs{d}(robust_non_skipped_mp_downs{d}) = mp_state_durations{d}{1}(robust_non_skipped_mp_downs{d}) - ...
        lfp_state_durations{d}{1}(corresp_lf8_downinds{d}(robust_non_skipped_mp_downs{d}));
    mp_rel_updurs_t2{d} = mp_rel_updurs{d}(rt2_ups{d});
    mp_rel_updurs_nt2{d} = mp_rel_updurs{d}(nrt2_ups{d});
    mp_rel_downdurs_t2{d} = mp_rel_downdurs{d}(rt2_downs{d});
    mp_rel_downdurs_nt2{d} = mp_rel_downdurs{d}(nrt2_downs{d});
    
    fract_rt1_ups(d) = nansum(mp_rel_updurs{d} > 0)/sum(~isnan(mp_rel_updurs{d}));
    fract_rt1_ups_nt2(d) = nansum(mp_rel_updurs_nt2{d} > 0)/sum(~isnan(mp_rel_updurs_nt2{d}));
    fract_rt1_downs(d) = nansum(mp_rel_downdurs{d} > 0)/sum(~isnan(mp_rel_downdurs{d}));
    
    %% compute durations in units of Ncx UDS cycles
    [mp_updurs_lfpc{d},mp_downdurs_lfpc{d}] = find_duration_ncx_uds_cycles(up_trans_inds,down_trans_inds,...
        mp_uplags{d},mp_downlags{d},lf8_period_vec);
    mp_corresp_lf8_mindur{d} = [mp_upskipped{d}.min_dur];
    mp_corresp_lf8_mindur{d}(nrt2_ups{d}) = lf8_updurs(corresp_lf8_upinds{d}(nrt2_ups{d}));
    rmp_updurs_lfpc{d} = mp_updurs_lfpc{d};
    rmp_updurs_lfpc{d}(mp_upskipped{d}.min_dur < thresh_lf8_downdur) = [];
    rmp_downdurs_lfpc{d} = mp_downdurs_lfpc{d};
    rmp_downdurs_lfpc{d}(mp_downskipped{d}.min_dur < thresh_lf8_updur) = [];
    
    %%
    median_uplag(d) = nanmedian(mp_uplags{d});
    median_uplag_t2(d) = nanmedian(mp_uplags{d}(rt2_ups{d}));
    median_uplag_nt2(d) = nanmedian(mp_uplags{d}(nrt2_ups{d}));
    median_reluplag_t2(d) = nanmedian(mp_reluplags{d}(rt2_ups{d}));
    median_reluplag_nt2(d) = nanmedian(mp_reluplags{d}(nrt2_ups{d}));
    median_downlag(d) = nanmedian(mp_downlags{d});
    median_downlag_t2(d) = nanmedian(mp_downlags{d}(rt2_ups{d}));
    median_downlag_nt2(d) = nanmedian(mp_downlags{d}(nrt2_ups{d}));
    median_reldownlag_t2(d) = nanmedian(mp_reldownlags{d}(rt2_ups{d}));
    median_reldownlag_nt2(d) = nanmedian(mp_reldownlags{d}(nrt2_ups{d}));
    mean_uplag(d) = nanmean(mp_uplags{d});
    mean_downlag(d) = nanmean(mp_downlags{d});
    median_reluplag(d) = nanmedian(mp_reluplags{d});
    median_reldownlag(d) = nanmedian(mp_reldownlags{d});
    mean_reluplag(d) = nanmean(mp_reluplags{d});
    mean_reldownlag(d) = nanmean(mp_reldownlags{d});
    
    %tests to see whether lags are different in pers vs non-pers states
    uplag_pers_test(d) = ranksum(mp_uplags{d}(rt2_ups{d}),mp_uplags{d}(nrt2_ups{d}));
    downlag_pers_test(d) = ranksum(mp_downlags{d}(rt2_ups{d}),mp_downlags{d}(nrt2_ups{d}));
    reluplag_pers_test(d) = ranksum(mp_reluplags{d}(rt2_ups{d}),mp_reluplags{d}(nrt2_ups{d}));
    reldownlag_pers_test(d) = ranksum(mp_reldownlags{d}(rt2_ups{d}),mp_reldownlags{d}(nrt2_ups{d}));
   
    %test whether up-transition lags are correlated with depth of
    %anesthesia properties
    robust_non_skipped_mp_ups{d}(end) = [];
    temp = corrcoef(mp_uplags{d}(robust_non_skipped_mp_ups{d}),lf8_updurs(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups{d})));
    uplag_updur8_corr(d) = temp(2,1);
    temp = corrcoef(mp_uplags{d}(robust_non_skipped_mp_ups{d}),lf8_dutycycs(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups{d})));
    uplag_dc8_corr(d) = temp(2,1);
    temp = corrcoef(mp_uplags{d}(robust_non_skipped_mp_ups{d}),lf8_downdurs(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups{d})));
    uplag_downdur8_corr(d) = temp(2,1);
    
    temp = corrcoef(mp_reluplags{d}(robust_non_skipped_mp_ups{d}),lf8_updurs(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups{d})));
    reluplag_updur8_corr(d) = temp(2,1);
    temp = corrcoef(mp_reluplags{d}(robust_non_skipped_mp_ups{d}),lf8_dutycycs(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups{d})));
    reluplag_dc8_corr(d) = temp(2,1);
    temp = corrcoef(mp_reluplags{d}(robust_non_skipped_mp_ups{d}),lf8_downdurs(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups{d})));
    reluplag_downdur8_corr(d) = temp(2,1);
    
end

cd C:\WC_Germany\sven_thomas_combined\
% save combined_core_analysis_lfpmprev_fin_nd
save distal_core_analysis_fin_nd_np

