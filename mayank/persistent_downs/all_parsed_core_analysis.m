
clear all
cd C:\WC_Germany\persistent_downs\
% load ./overall_EC_dir
load ./new_pdown_dir

addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
addpath('C:\WC_Germany\overall_EC')
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

seg_win_dur = round(Fsd*100);


%%
% for d = 1:length(sess_data)
%     cd(sess_data(d).directory)
uset = 1:length(new_pdown_dir);
uset(12) = [];
for d = uset
    cd(new_pdown_dir{d})
    pwd
    
    load ./used_data lf7 wcv_minus_spike
    [lf8_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    
    load ./all_combined_mp_uds.mat
    load ./all_combined_lf7_uds.mat
    
    if ~isempty(hmm7)
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
        fract_uds_dur(d) = dur_uds/t_axis(end);
        
        %%
        load ./allEC_ctx_period_data
        
        lf8_period_vec = nan(size(wcv_lf));
        for i = 1:size(new_seg_inds,1)
            lf8_period_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = lf8_period_f{i};
        end
        
        %% compute state duration distributions
        [mp_state_durations{d}] = compute_state_durations_seg(mp_state_seq,Fsd);
        [lfp_state_durations{d}] = compute_state_durations_seg(lfp_state_seq,Fsd);
        mp_updurs = mp_state_durations{d}{2};
        mp_downdurs = mp_state_durations{d}{1};
        lf8_updurs= lfp_state_durations{d}{2};
        lf8_downdurs = lfp_state_durations{d}{1};
        
        mp_dutycycs = mp_updurs./(mp_updurs+mp_downdurs);
        lf8_dutycycs = lf8_updurs./(lf8_updurs+lf8_downdurs);
        
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
        
        rt2_downs{d} = find(mp_downskipped{d}.rnum_skipped > 0);
        nrt2_downs{d} = find(mp_downskipped{d}.num_skipped==0); %number of down states that didn't skip any LFP up states
        robust_non_skipped_mp_ups{d} = [rt2_ups{d}; nrt2_ups{d}];
        robust_non_skipped_mp_downs{d} = [rt2_downs{d}; nrt2_downs{d}];
        
        mp_rel_updurs{d} = nan(size(up_trans_inds));
        mp_rel_downdurs{d} = nan(size(up_trans_inds));
        mp_rel_updurs{d}(robust_non_skipped_mp_ups{d}) = mp_state_durations{d}{2}(robust_non_skipped_mp_ups{d}) - ...
            lfp_state_durations{d}{2}(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups{d}));
        mp_rel_downdurs{d}(robust_non_skipped_mp_downs{d}) = mp_state_durations{d}{1}(robust_non_skipped_mp_downs{d}) - ...
            lfp_state_durations{d}{1}(corresp_lf8_downinds{d}(robust_non_skipped_mp_downs{d}));
        
        %%
        tot_dur = length(t_axis)/Fsd;
        min_states = 5;
        n_segs(d) = floor(tot_dur/seg_win_dur*Fsd);
        mp_up_times = t_axis(up_trans_inds);
        mp_down_times = t_axis(down_trans_inds);
        lfp_up_times = t_axis(up_trans_inds8);
        lfp_down_times = t_axis(down_trans_inds8);
        
        seg_mp_updur{d} = nan(n_segs(d),1);
        seg_mp_downdur{d} = nan(n_segs(d),1);
        seg_lfp_updur{d} = nan(n_segs(d),1);
        seg_lfp_downdur{d} = nan(n_segs(d),1);
        seg_mp_relupdur{d} = nan(n_segs(d),1);
        seg_mp_reldowndur{d} = nan(n_segs(d),1);
        seg_mp_uplag{d} = nan(n_segs(d),1);
        seg_mp_downlag{d} = nan(n_segs(d),1);
        seg_fract_rt2_ups{d} = nan(n_segs(d),1);
        seg_fract_rt2_downs{d} = nan(n_segs(d),1);
        for ss = 1:n_segs(d)
            time_bounds = [(ss-1)*seg_win_dur ss*seg_win_dur]/Fsd;
            cur_mp_utrans = find(mp_up_times >= time_bounds(1) & mp_up_times < time_bounds(2));
            cur_lfp_utrans = find(lfp_up_times >= time_bounds(1) & lfp_up_times < time_bounds(2));
            n_seg_mp_ups{d}(ss) = length(cur_mp_utrans);
            n_seg_lfp_ups{d}(ss) = length(cur_lfp_utrans);
            
            if n_seg_mp_ups{d}(ss) > min_states & n_seg_lfp_ups{d}(ss) > min_states
                seg_mp_updur{d}(ss) = mean(mp_updurs(cur_mp_utrans));
                seg_mp_downdur{d}(ss) =mean(mp_downdurs(cur_mp_utrans));
                seg_lfp_updur{d}(ss) = mean(lf8_updurs(cur_lfp_utrans));
                seg_lfp_downdur{d}(ss) = mean(lf8_downdurs(cur_lfp_utrans));
                
                seg_mp_relupdur{d}(ss) = nanmean(mp_rel_updurs{d}(cur_mp_utrans));
                seg_mp_reldowndur{d}(ss) = nanmean(mp_rel_downdurs{d}(cur_mp_utrans));
                
                seg_mp_uplag{d}(ss) = nanmean(mp_uplags{d}(cur_mp_utrans));
                seg_mp_downlag{d}(ss) = nanmean(mp_downlags{d}(cur_mp_utrans));
                
                seg_fract_rt2_ups{d}(ss) = sum(ismember(cur_mp_utrans,rt2_ups{d}))/length(cur_mp_utrans);
                seg_fract_rt2_downs{d}(ss) = sum(ismember(cur_mp_utrans,rt2_downs{d}))/length(cur_mp_utrans);
            end
        end
        
        bad_sets(d) = 0;
    else
        bad_sets(d) = 1;
    end
end

cd C:\WC_Germany\persistent_downs\
save allEC_segcore_analysis
% save new_down_core_analysis

%%
cd C:\WC_Germany\persistent_downs\
load ./allEC_segcore_analysis

close all
for d = 29:35
    uset(d)
    tt = (1:n_segs(d))*seg_win_dur/Fsd;
    subplot(2,3,1)
    plot(tt,seg_fract_rt2_ups{d},'o-')
    subplot(2,3,4)
    plot(tt,seg_fract_rt2_downs{d},'o-')

%         subplot(2,3,2)
%     plot(tt,seg_mp_relupdur{d},'o-')
%     subplot(2,3,5)
%     plot(tt,seg_mp_reldowndur{d},'o-')
        subplot(2,3,2)
    plot(tt,seg_mp_updur{d},'o-')
hold on
    plot(tt,seg_lfp_updur{d},'ro-')
subplot(2,3,5)
    plot(tt,seg_mp_downdur{d},'o-')
hold on
    plot(tt,seg_lfp_downdur{d},'ro-')

        subplot(2,3,3)
    plot(tt,seg_mp_uplag{d}/Fsd,'o-')
    subplot(2,3,6)
    plot(tt,seg_mp_downlag{d}/Fsd,'o-')
    pause
    clf
end