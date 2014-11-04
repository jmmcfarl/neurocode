
clear all
close all

cur_dir_letter = 'C';
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
cd C:\WC_Germany\sven_thomas_combined\

load .\polar_dir


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
cnt = 1;
for d = 1:length(polar_dir)
    for j = 1:length(polar_dir(d).subdir)
        cd(polar_dir(d).dir)
        cd(polar_dir(d).subdir(j).dir)
        
        load ./used_data lf7
        lf8 = lf7;
        [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
        %     [lf82_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[0.05 2]);
        wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
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
        fract_uds_dur(d,j) = dur_uds/t_axis(end);
        
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
        [mp_state_durations{cnt}] = compute_state_durations_seg(hsmm_bbstate_seq,Fsd);
        [lfp_state_durations{cnt}] = compute_state_durations_seg(hsmm_bbstate_seq8,Fsd);
        
        [up_amps,down_amps] = get_state_amplitudes(lf8_lf,up_trans_inds8,down_trans_inds8);
        
        up_dur_hist(cnt,:) = histc(mp_state_durations{cnt}{2},lin_dur_grid);
        up_dur_hist8(cnt,:) = histc(lfp_state_durations{cnt}{2},lin_dur_grid);
        down_dur_hist(cnt,:) = histc(mp_state_durations{cnt}{1},lin_dur_grid);
        down_dur_hist8(cnt,:) = histc(lfp_state_durations{cnt}{1},lin_dur_grid);
        up_dur_loghist(cnt,:) = histc(mp_state_durations{cnt}{2},log_dur_grid);
        up_dur_loghist8(cnt,:) = histc(lfp_state_durations{cnt}{2},log_dur_grid);
        down_dur_loghist(cnt,:) = histc(mp_state_durations{cnt}{1},log_dur_grid);
        down_dur_loghist8(cnt,:) = histc(lfp_state_durations{cnt}{1},log_dur_grid);
        
        mp_updurs = mp_state_durations{cnt}{2};
        mp_downdurs = mp_state_durations{cnt}{1};
        lf8_updurs= lfp_state_durations{cnt}{2};
        lf8_downdurs = lfp_state_durations{cnt}{1};
        
        mp_dutycycs = mp_updurs./(mp_updurs+mp_downdurs);
        lf8_dutycycs = lf8_updurs./(lf8_updurs+lf8_downdurs);
        median_dc(cnt) = nanmedian(mp_dutycycs);
        median_dc8(cnt) = nanmedian(lf8_dutycycs);
        mean_dc(cnt) = nanmean(mp_dutycycs);
        mean_dc8(cnt) = nanmean(lf8_dutycycs);
        
        mean_updur8(cnt) = nanmean(lfp_state_durations{cnt}{2});
        mean_downdur8(cnt) = nanmean(lfp_state_durations{cnt}{1});
        median_updur8(cnt) = nanmedian(lfp_state_durations{cnt}{2});
        median_downdur8(cnt) = nanmedian(lfp_state_durations{cnt}{1});
        mean_updur(cnt) = nanmean(mp_state_durations{cnt}{2});
        mean_downdur(cnt) = nanmean(mp_state_durations{cnt}{1});
        median_updur(cnt) = nanmedian(mp_state_durations{cnt}{2});
        median_downdur(cnt) = nanmedian(mp_state_durations{cnt}{1});
        max_updur8(cnt) = nanmax(lfp_state_durations{cnt}{2});
        max_updur(cnt) = nanmax(mp_state_durations{cnt}{2});
        max_downdur8(cnt) = nanmax(lfp_state_durations{cnt}{1});
        max_downdur(cnt) = nanmax(mp_state_durations{cnt}{1});
        
        
        
        %% compute transition-triggered averages
        lf8_utrig_lf8(cnt,:) = get_event_trig_avg(lf8_lf,up_trans_inds8,forwardlag,backlag);
        lf8_dtrig_lf8(cnt,:) = get_event_trig_avg(lf8_lf,down_trans_inds8,forwardlag,backlag);
        lf8_utrig_wcv(cnt,:) = get_event_trig_avg(wcv_lf,up_trans_inds8,forwardlag,backlag);
        lf8_dtrig_wcv(cnt,:) = get_event_trig_avg(wcv_lf,down_trans_inds8,forwardlag,backlag);
        %     lf8_utrig_hpc(d,:) = get_event_trig_avg(hpc_lf,up_trans_inds8,forwardlag,backlag);
        %     lf8_dtrig_hpc(d,:) = get_event_trig_avg(hpc_lf,down_trans_inds8,forwardlag,backlag);
        %     lf8_utrig_hpchf(d,:) = get_event_trig_avg(hpc_hf,up_trans_inds8,forwardlag,backlag);
        %     lf8_dtrig_hpchf(d,:) = get_event_trig_avg(hpc_hf,down_trans_inds8,forwardlag,backlag);
        
        wcv_utrig_lf8(cnt,:) = get_event_trig_avg(lf8_lf,up_trans_inds,forwardlag,backlag);
        wcv_dtrig_lf8(cnt,:) = get_event_trig_avg(lf8_lf,down_trans_inds,forwardlag,backlag);
        wcv_utrig_wcv(cnt,:) = get_event_trig_avg(wcv_lf,up_trans_inds,forwardlag,backlag);
        wcv_dtrig_wcv(cnt,:) = get_event_trig_avg(wcv_lf,down_trans_inds,forwardlag,backlag);
        %     wcv_utrig_hpc(d,:) = get_event_trig_avg(hpc_lf,up_trans_inds,forwardlag,backlag);
        %     wcv_dtrig_hpc(d,:) = get_event_trig_avg(hpc_lf,down_trans_inds,forwardlag,backlag);
        %     wcv_utrig_hpchf(d,:) = get_event_trig_avg(hpc_hf,up_trans_inds,forwardlag,backlag);
        %     wcv_dtrig_hpchf(d,:) = get_event_trig_avg(hpc_hf,down_trans_inds,forwardlag,backlag);
        
        lags = (-backlag:forwardlag)/Fsd;
        
        %% compute corresponding state transitions and transition lags
        [corresp_lf8_upinds{cnt},corresp_lf8_downinds{cnt}] = greedy_find_corresponding_ncx_state_transitions_simp(...
            up_trans_inds,down_trans_inds,up_trans_inds8,down_trans_inds8);
        
        %find non-skipped states
        non_skipped_mp_up_states{cnt} = find(~isnan(corresp_lf8_upinds{cnt}));
        non_skipped_mp_down_states{cnt} = find(~isnan(corresp_lf8_downinds{cnt}));
        
        %compute transition lags for non-skipped states
        mp_uplags{cnt} = nan(size(up_trans_inds));
        mp_downlags{cnt} = nan(size(up_trans_inds));
        mp_uplags{cnt}(non_skipped_mp_up_states{cnt}) = up_trans_inds(non_skipped_mp_up_states{cnt}) - up_trans_inds8(corresp_lf8_upinds{cnt}(non_skipped_mp_up_states{cnt}));
        mp_downlags{cnt}(non_skipped_mp_down_states{cnt}) = down_trans_inds(non_skipped_mp_down_states{cnt}) - down_trans_inds8(corresp_lf8_downinds{cnt}(non_skipped_mp_down_states{cnt}));
        
        %compute transition lags relative to ctx state durations
        mp_reldownlags{cnt} = nan(size(up_trans_inds));
        mp_reluplags{cnt} = nan(size(up_trans_inds));
        for i = 1:length(up_trans_inds)
            if ~isnan(corresp_lf8_downinds{cnt}(i))
                mp_reldownlags{cnt}(i) = mp_downlags{cnt}(i)/lf8_downdurs(corresp_lf8_downinds{cnt}(i))/Fsd;
            end
            if ~isnan(corresp_lf8_upinds{cnt}(i))
                mp_reluplags{cnt}(i) = mp_uplags{cnt}(i)/lf8_updurs(corresp_lf8_upinds{cnt}(i))/Fsd;
            end
        end
        
        %% compute persistence
        [mp_upskipped{cnt},mp_downskipped{cnt}] = greedy_find_skipped_ncx_states(...
            corresp_lf8_upinds{cnt},corresp_lf8_downinds{cnt},lfp_state_durations{cnt}{2},lfp_state_durations{cnt}{1},thresh_lf8_downdur,thresh_lf8_updur);
        
        rt2_ups{cnt} = find(mp_upskipped{cnt}.rnum_skipped > 0);
        t2_ups{cnt} = find(mp_upskipped{cnt}.num_skipped > 0);
        nrt2_ups{cnt} = find(mp_upskipped{cnt}.num_skipped == 0);
        fract_rt2_ups(cnt) = length(rt2_ups{cnt})/(length(rt2_ups{cnt}) + length(nrt2_ups{cnt}));
        fract_t2_ups(cnt) = length(t2_ups{cnt})/(length(t2_ups{cnt}) + length(nrt2_ups{cnt}));
        
        rt2_downs{cnt} = find(mp_downskipped{cnt}.rnum_skipped > 0);
        nrt2_downs{cnt} = find(mp_downskipped{cnt}.num_skipped==0); %number of down states that didn't skip any LFP up states
        fract_rt2_downs(cnt) = length(rt2_downs{cnt})/(length(rt2_downs{cnt}) + length(nrt2_downs{cnt}));
        
        robust_non_skipped_mp_ups{cnt} = [rt2_ups{cnt}; nrt2_ups{cnt}];
        robust_non_skipped_mp_downs{cnt} = [rt2_downs{cnt}; nrt2_downs{cnt}];
        
        mp_rel_updurs{cnt} = nan(size(up_trans_inds));
        mp_rel_downdurs{cnt} = nan(size(up_trans_inds));
        mp_rel_updurs{cnt}(robust_non_skipped_mp_ups{cnt}) = mp_state_durations{cnt}{2}(robust_non_skipped_mp_ups{cnt}) - ...
            lfp_state_durations{cnt}{2}(corresp_lf8_upinds{cnt}(robust_non_skipped_mp_ups{cnt}));
        mp_rel_downdurs{cnt}(robust_non_skipped_mp_downs{cnt}) = mp_state_durations{cnt}{1}(robust_non_skipped_mp_downs{cnt}) - ...
            lfp_state_durations{cnt}{1}(corresp_lf8_downinds{cnt}(robust_non_skipped_mp_downs{cnt}));
        mp_rel_updurs_t2{cnt} = mp_rel_updurs{cnt}(rt2_ups{cnt});
        mp_rel_updurs_nt2{cnt} = mp_rel_updurs{cnt}(nrt2_ups{cnt});
        mp_rel_downdurs_t2{cnt} = mp_rel_downdurs{cnt}(rt2_downs{cnt});
        mp_rel_downdurs_nt2{cnt} = mp_rel_downdurs{cnt}(nrt2_downs{cnt});
        
        fract_rt1_ups(cnt) = nansum(mp_rel_updurs{cnt} > 0)/sum(~isnan(mp_rel_updurs{cnt}));
        fract_rt1_ups_nt2(cnt) = nansum(mp_rel_updurs_nt2{cnt} > 0)/sum(~isnan(mp_rel_updurs_nt2{cnt}));
        fract_rt1_downs(cnt) = nansum(mp_rel_downdurs{cnt} > 0)/sum(~isnan(mp_rel_downdurs{cnt}));
        
        %% compute durations in units of Ncx UDS cycles
        [mp_updurs_lfpc{cnt},mp_downdurs_lfpc{cnt}] = find_duration_ncx_uds_cycles(up_trans_inds,down_trans_inds,...
            mp_uplags{cnt},mp_downlags{cnt},lf8_period_vec);
        mp_corresp_lf8_mindur{cnt} = [mp_upskipped{cnt}.min_dur];
        mp_corresp_lf8_mindur{cnt}(nrt2_ups{cnt}) = lf8_updurs(corresp_lf8_upinds{cnt}(nrt2_ups{cnt}));
        rmp_updurs_lfpc{cnt} = mp_updurs_lfpc{cnt};
        rmp_updurs_lfpc{cnt}(mp_upskipped{cnt}.min_dur < thresh_lf8_downdur) = [];
        rmp_downdurs_lfpc{cnt} = mp_downdurs_lfpc{cnt};
        rmp_downdurs_lfpc{cnt}(mp_downskipped{cnt}.min_dur < thresh_lf8_updur) = [];
        
        %%
        median_uplag(cnt) = nanmedian(mp_uplags{cnt});
        median_uplag_t2(cnt) = nanmedian(mp_uplags{cnt}(rt2_ups{cnt}));
        median_uplag_nt2(cnt) = nanmedian(mp_uplags{cnt}(nrt2_ups{cnt}));
        median_reluplag_t2(cnt) = nanmedian(mp_reluplags{cnt}(rt2_ups{cnt}));
        median_reluplag_nt2(cnt) = nanmedian(mp_reluplags{cnt}(nrt2_ups{cnt}));
        median_downlag(cnt) = nanmedian(mp_downlags{cnt});
        median_downlag_t2(cnt) = nanmedian(mp_downlags{cnt}(rt2_ups{cnt}));
        median_downlag_nt2(cnt) = nanmedian(mp_downlags{cnt}(nrt2_ups{cnt}));
        median_reldownlag_t2(cnt) = nanmedian(mp_reldownlags{cnt}(rt2_ups{cnt}));
        median_reldownlag_nt2(cnt) = nanmedian(mp_reldownlags{cnt}(nrt2_ups{cnt}));
        mean_uplag(cnt) = nanmean(mp_uplags{cnt});
        mean_downlag(cnt) = nanmean(mp_downlags{cnt});
        median_reluplag(cnt) = nanmedian(mp_reluplags{cnt});
        median_reldownlag(cnt) = nanmedian(mp_reldownlags{cnt});
        mean_reluplag(cnt) = nanmean(mp_reluplags{cnt});
        mean_reldownlag(cnt) = nanmean(mp_reldownlags{cnt});
        
        %tests to see whether lags are different in pers vs non-pers states
        uplag_pers_test(cnt) = ranksum(mp_uplags{cnt}(rt2_ups{cnt}),mp_uplags{cnt}(nrt2_ups{cnt}));
        downlag_pers_test(cnt) = ranksum(mp_downlags{cnt}(rt2_ups{cnt}),mp_downlags{cnt}(nrt2_ups{cnt}));
        reluplag_pers_test(cnt) = ranksum(mp_reluplags{cnt}(rt2_ups{cnt}),mp_reluplags{cnt}(nrt2_ups{cnt}));
        reldownlag_pers_test(cnt) = ranksum(mp_reldownlags{cnt}(rt2_ups{cnt}),mp_reldownlags{cnt}(nrt2_ups{cnt}));
        
        %test whether up-transition lags are correlated with depth of
        %anesthesia properties
        robust_non_skipped_mp_ups{cnt}(end) = [];
        temp = corrcoef(mp_uplags{cnt}(robust_non_skipped_mp_ups{cnt}),lf8_updurs(corresp_lf8_upinds{cnt}(robust_non_skipped_mp_ups{cnt})));
        uplag_updur8_corr(cnt) = temp(2,1);
        temp = corrcoef(mp_uplags{cnt}(robust_non_skipped_mp_ups{cnt}),lf8_dutycycs(corresp_lf8_upinds{cnt}(robust_non_skipped_mp_ups{cnt})));
        uplag_dc8_corr(cnt) = temp(2,1);
        temp = corrcoef(mp_uplags{cnt}(robust_non_skipped_mp_ups{cnt}),lf8_downdurs(corresp_lf8_upinds{cnt}(robust_non_skipped_mp_ups{cnt})));
        uplag_downdur8_corr(cnt) = temp(2,1);
        
        temp = corrcoef(mp_reluplags{cnt}(robust_non_skipped_mp_ups{cnt}),lf8_updurs(corresp_lf8_upinds{cnt}(robust_non_skipped_mp_ups{cnt})));
        reluplag_updur8_corr(cnt) = temp(2,1);
        temp = corrcoef(mp_reluplags{cnt}(robust_non_skipped_mp_ups{cnt}),lf8_dutycycs(corresp_lf8_upinds{cnt}(robust_non_skipped_mp_ups{cnt})));
        reluplag_dc8_corr(cnt) = temp(2,1);
        temp = corrcoef(mp_reluplags{cnt}(robust_non_skipped_mp_ups{cnt}),lf8_downdurs(corresp_lf8_upinds{cnt}(robust_non_skipped_mp_ups{cnt})));
        reluplag_downdur8_corr(cnt) = temp(2,1);
        
        dvals(cnt) = d;
        jvals(cnt) = j;
        cnt = cnt + 1;
    end
end
cd C:\WC_Germany\sven_thomas_combined\
save pol_core_analysis_fin_nd_np

