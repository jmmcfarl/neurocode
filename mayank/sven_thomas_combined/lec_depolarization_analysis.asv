clear all
close all

lec_ci_dir{1} = 'C:\wc_data\2012_06_21_A\2012-6-21_14-33-23'; %before 900 = depol, after 950 = spontaneous
lec_ci_dir{2} = 'C:\wc_data\2012_06_21_A\2012-6-21_14-33-23'; %before 900 = depol, after 950 = spontaneous

lec_ci_dir{3} = 'C:\wc_data\2012_06_23_C\2012-6-23_19-39-37';
lec_ci_dir{4} = 'C:\wc_data\2012_06_23_C\2012-6-23_19-45-51';

lec_ci_dir{5} = 'C:\wc_data\2012_06_24_A\2012-6-24_14-49-46';
lec_ci_dir{6} = 'C:\wc_data\2012_06_24_A\2012-6-24_15-1-29';

lec_ci_dir{7} = 'C:\wc_data\2012_06_24_B\2012-6-24_15-46-53';
lec_ci_dir{8} = 'C:\wc_data\2012_06_24_B\2012-6-24_15-54-33';

lec_ci_dir{9} = 'C:\wc_data\2012_06_24_C\2012-6-24_16-24-8';
lec_ci_dir{10} = 'C:\wc_data\2012_06_24_C\2012-6-24_16-17-53';

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

up_fract = 0.5;
down_phase = 180;

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
numBins = 60;
lin_dur_grid = linspace(up_range(1),up_range(2),numBins+1);
log_dur_grid = logspace(log10(up_range(1)),log10(up_range(2)),numBins+1);

%%
for d = 1:10
    cd(lec_ci_dir{d})
    pwd
    
    
    load ./used_data lf5 wcv_minus_spike
    lf8 = lf5;
    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    
    load ./pa_hsmm_state_seq_combined_fin_nd.mat
    load ./pa_hsmm_state_seq7_combined_fin_nd.mat
    hsmm_bbstate_seq8 = hsmm_bbstate_seq7;
    hsmm8 = hsmm7;
    
    load ./spike_time_jmm
    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    if d < 3
       range1 = find(t_axis < 900);
       range2 = find(t_axis > 950);
       if d == 2
           lf8_lf = lf8_lf(range1);
           t_axis = t_axis(range1);
           wcv_lf = wcv_lf(range1);
           hsmm_bbstate_seq8{1} = hsmm_bbstate_seq8{1}(range1);
           hsmm_bbstate_seq{1} = hsmm_bbstate_seq{1}(range1);
           synct_d = synct_d(range1);
       else 
             lf8_lf = lf8_lf(range2);
           t_axis = t_axis(range2);
           wcv_lf = wcv_lf(range2);
           hsmm_bbstate_seq8{1} = hsmm_bbstate_seq8{1}(range2);
           hsmm_bbstate_seq{1} = hsmm_bbstate_seq{1}(range2);
           synct_d = synct_d(range2);
       end
    end
    
    mp_spike_rate = hist(synct(spkid),synct_d)*Fsd; 
    mp_spike_rate([1 end]) = 0;
    wcv_up_log = nan(size(synct_d));
    lf8_up_log = nan(size(synct_d));
    
        [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,hsmm_bbstate_seq);
for ns = 1:hmm.Nsegs
        cur_seg = new_seg_inds(ns,1):new_seg_inds(ns,2);
        wcv_up_log(cur_seg) = logical(hsmm_bbstate_seq{ns}-1);
        lf8_up_log(cur_seg) = logical(hsmm_bbstate_seq8{ns}-1);
    end

    mp_avg_rate(d) = mean(mp_spike_rate(~isnan(wcv_up_log)));
    mp_up_rate(d) = mean(mp_spike_rate(wcv_up_log == 1));
    
    %         if you want to flip the roles of the MP and LFP
    %         temp = hsmm_bbstate_seq8;
    %         hsmm_bbstate_seq8 = hsmm_bbstate_seq;
    %         hsmm_bbstate_seq = temp;
    
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq8);
    fract_uds_dur(d) = dur_uds/t_axis(end);
        
    lfp_state_seq = hsmm_bbstate_seq8;
    seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
    for ns = 1:hsmm8.Nsegs
        cur_up_trans = find(lfp_state_seq{ns}(1:end-1) == 1 & lfp_state_seq{ns}(2:end) == 2);
        cur_down_trans = find(lfp_state_seq{ns}(1:end-1) == 2 & lfp_state_seq{ns}(2:end) == 1);
        cur_up_trans(cur_up_trans > cur_down_trans(end)) = [];
        cur_down_trans(cur_down_trans < cur_up_trans(1)) = [];
        
        lf8_period_f{ns} = nan(seg_durs(ns),1);
        lf8_period_p{ns} = nan(seg_durs(ns),1);
        lf8_uds_amp{ns} = nan(length(cur_up_trans)-1,1);
        for i = 1:length(cur_up_trans)-1
            cur_up_times = cur_up_trans(i):cur_down_trans(i);
            cur_down_times = cur_down_trans(i):cur_up_trans(i+1);
            max_up = max(lf8_lf(cur_up_times));
            min_down = min(lf8_lf(cur_down_times));
            lf8_uds_amp{ns}(i) = max_up-min_down;
            
            period_samps_up = cur_down_trans(i)-cur_up_trans(i);
            period_samps_down = cur_up_trans(i+1)-cur_down_trans(i);
            
            lf8_period_f{ns}(cur_up_trans(i)+1:cur_down_trans(i)) = ...
                i-1+linspace(1,period_samps_up,period_samps_up)/period_samps_up*up_fract;
            lf8_period_f{ns}(cur_down_trans(i)+1:cur_up_trans(i+1)) = ...
                i-1+up_fract+linspace(1,period_samps_down,period_samps_down)/period_samps_down*(1-up_fract);
            
            lf8_period_p{ns}(cur_up_trans(i)+1:cur_down_trans(i)) = ...
                linspace(1,period_samps_up,period_samps_up)/period_samps_up*down_phase;
            lf8_period_p{ns}(cur_down_trans(i)+1:cur_up_trans(i+1)) = ...
                down_phase+linspace(1,period_samps_down,period_samps_down)/period_samps_down*(360-down_phase);           
        end
    end
    
    %%
    [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,hsmm_bbstate_seq);
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq8);
    fract_uds_dur(d) = dur_uds/t_axis(end);
    
    %%
%     load ./combined_lf7_period_data_fin_nd %
    %     load ./combined_mp_period_data_fin_nd %to reverse roles of MP and LFP
    
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
    
%     fract_rt1_ups(d) = nansum(mp_rel_updurs{d} > 0)/sum(~isnan(mp_rel_updurs{d}));
%     fract_rt1_ups_nt2(d) = nansum(mp_rel_updurs_nt2{d} > 0)/sum(~isnan(mp_rel_updurs_nt2{d}));
%     fract_rt1_downs(d) = nansum(mp_rel_downdurs{d} > 0)/sum(~isnan(mp_rel_downdurs{d}));
    %%
    median_uplag(d) = nanmedian(mp_uplags{d});
    median_uplag_t2(d) = nanmedian(mp_uplags{d}(rt2_ups{d}));
    median_uplag_nt2(d) = nanmedian(mp_uplags{d}(nrt2_ups{d}));
    median_downlag(d) = nanmedian(mp_downlags{d});
    median_downlag_t2(d) = nanmedian(mp_downlags{d}(rt2_ups{d}));
    median_downlag_nt2(d) = nanmedian(mp_downlags{d}(nrt2_ups{d}));
    mean_uplag(d) = nanmean(mp_uplags{d});
    mean_downlag(d) = nanmean(mp_downlags{d});
    median_reluplag(d) = nanmedian(mp_reluplags{d});
end

%%
% %%
% close all
% raw_Fs = 2016;
% dsf = 8;
% Fsd = raw_Fs/dsf;
% niqf = raw_Fs/2;
% hcf = 4;
% lcf = 0.05;
% lcf_hf = 15;
% hcf_hf = 80;
% hcf_sm = 0.025;
% load ./used_data lf5 wcv_minus_spike
% lf8 = lf5;
% [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
% wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
% hold on
% plot(t_axis,wcv_lf+3)
% hold on
% plot(t_axis,lf8_lf,'k')
% 
% load ./pa_hsmm_state_seq_combined_fin_nd.mat
% load ./pa_hsmm_state_seq7_combined_fin_nd.mat
% hsmm_bbstate_seq8 = hsmm_bbstate_seq7;
% hsmm8 = hsmm7;
% [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,hsmm_bbstate_seq);
% dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
% [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
% [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq8);
% 
% plot(t_axis(up_trans_inds8),lf8_lf(up_trans_inds8),'ro')
% plot(t_axis(down_trans_inds8),lf8_lf(down_trans_inds8),'go')
% plot(t_axis(up_trans_inds),wcv_lf(up_trans_inds)+3,'r*')
% plot(t_axis(down_trans_inds),wcv_lf(down_trans_inds)+3,'g*')
