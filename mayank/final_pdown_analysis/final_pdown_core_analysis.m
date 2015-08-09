% addpath('C:/WC_Germany/parietal_cortical_2010/');
% addpath('C:/Code/general_functions/');
% addpath('C:/WC_Germany/persistent_9_27_2010/');
% addpath('C:/WC_Germany/hsmm_state_detection/');
% addpath('C:/WC_Germany/persistent_downs/');
addpath('~/James_scripts/mayank/parietal_cortical_2010/');
addpath('~/James_scripts/mayank/persistent_9_27_2010/');
addpath('~/James_scripts/mayank/hsmm_state_detection/');
addpath('~/James_scripts/mayank/persistent_downs/');
addpath('~/James_scripts/mayank/final_pdown_analysis/');

clear all
% load C:/WC_Germany/final_pdown_analysis/compiled_data.mat
load ~/Analysis/Mayank/final_pdown_analysis/compiled_data.mat

min_rec_dur = 500; %minimum recording duration (in sec)
data_ids = [data(:).id];
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur); %only analyze recs where we have at least this much data to work with
% used_dirs([data(used_dirs).id] == 76) = []; %this rec does not have usable LFP UDS
used_dirs(ismember(data_ids(used_dirs),no_cell) | ismember(data_ids(used_dirs),unclear_uds)) = [];
% used_dirs(ismember(data_ids(used_dirs),unclear_uds)) = [];

data = data(used_dirs);
data_ids = data_ids(used_dirs);

%%
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
lcf = 0.05; hcf = 10; %main low-pass filter range 
% lcf_hf = 15; hcf_hf = 80; %filter range for computing HF power
% hcf_sm = 0.025; %smoothing sigma (s) for computing HF power
% rate_sm = round(Fsd*0.05); %smoothing sigma for computing firing rates
 
% maxlag = round(4*Fsd);
% backlag = 4*Fsd;
% forwardlag = 4*Fsd;

% params.Fs = raw_Fs;
% params.fpass = [0 100];
% params.tapers = [3 5];
% win = 25;

%robust persistence requires MP skips states longer than this
thresh_lfp_updur = 0.5;
thresh_lfp_downdur = 0.5;

% %for state duration distribution calculations
% up_range = [0.1 15];
% down_range = [0.1 15];
% numBins = 60;
% log_dur_grid = logspace(log10(up_range(1)),log10(up_range(2)),numBins+1);

%%
for d = 1:length(data)
%     cd(data(d).dir)
    cd(data(d).new_dir);
pwd
    
    core_data(d).data_id = data(d).id;
    
    load ./used_data lf7 wcv_minus_spike
    
    %get low-pass filtered ctx LFP and MP signals
    [lfp_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    
    %get synced time stamps, at desired temporal resolution
    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    
    %cut data at the soonest of ep or dp (when data becomes unstable or
    %some other experimental manipulation is applied
    end_time = min(data(d).ep,data(d).dp);
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        synct_d(ep+1:end) = []; lfp_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; 
    else
        ep = length(t_axis);
    end
    
    %load the classified state seq data
    load ./pa_hsmm_state_seq_combined_fin_nd.mat
    load ./pa_hsmm_state_seq7_combined_fin_nd.mat
    
    %         if you want to flip the roles of the MP and LFP
    %         temp = hsmm_bbstate_seq7;
    %         hsmm_bbstate_seq7 = hsmm_bbstate_seq;
    %         hsmm_bbstate_seq = temp;

    lfp_hsmm = hsmm7;
    lfp_state_seq = hsmm_bbstate_seq7;
    mp_state_seq = hsmm_bbstate_seq;
    
    %resample UDS sequences and find state transition times
    [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,mp_state_seq);
    core_data(d).tot_uds_dur = sum(diff(new_seg_inds,[],2))/Fsd; %total duration of UDS
    
    %get state transition times for lfp and MP
    [up_trans_inds_mp,down_trans_inds_mp] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds_lfp,down_trans_inds_lfp] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    %eliminate state transitions after the data end point
    bad_mp_states = find(up_trans_inds_mp > ep | down_trans_inds_mp > ep);
    up_trans_inds_mp(bad_mp_states) = []; down_trans_inds_mp(bad_mp_states) = [];
    bad_lfp_states = find(up_trans_inds_lfp > ep | down_trans_inds_lfp > ep);
    up_trans_inds_lfp(bad_lfp_states) = []; down_trans_inds_lfp(bad_lfp_states) = [];
        
    %% LOAD LFP UDS period data for computing relative state durations
    load ./allEC_ctx_period_data_hsmm.mat
    lfp_period_vec = nan(size(wcv_lf));
    for i = 1:size(new_seg_inds,1)
        cur_inds = new_seg_inds(i,1):new_seg_inds(i,2);
        cur_inds_used = find(cur_inds <= ep);
        lfp_period_vec(cur_inds(cur_inds_used)) = lf8_period_f{i}(cur_inds_used);
    end
    
    %% COMPUTE STATE VECTORS
    mp_state_number = nan(length(wcv_lf),1);
    mp_state_vec = zeros(length(wcv_lf),1);
    lfp_state_number = nan(length(wcv_lf),1);
    lfp_state_vec = zeros(length(wcv_lf),1);
    for i = 1:length(up_trans_inds_mp)-1
        mp_state_vec(up_trans_inds_mp(i):down_trans_inds_mp(i)) = 1;
        mp_state_number(up_trans_inds_mp(i):down_trans_inds_mp(i)) = 2*(i-1)+1;
        mp_state_number(down_trans_inds_mp(i):up_trans_inds_mp(i+1)) = 2*(i-1)+2;
    end
    for i = 1:length(up_trans_inds_lfp)-1
        lfp_state_vec(up_trans_inds_lfp(i):down_trans_inds_lfp(i)) = 1;
        lfp_state_number(up_trans_inds_lfp(i):down_trans_inds_lfp(i)) = 2*(i-1)+1;
        lfp_state_number(down_trans_inds_lfp(i):up_trans_inds_lfp(i+1)) = 2*(i-1) + 2;
    end
    mp_state_number(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;
    mp_state_vec(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;

    %% compute state duration distributions
    [mp_state_durations] = compute_state_durations_seg(mp_state_seq,Fsd);
    [lfp_state_durations] = compute_state_durations_seg(lfp_state_seq,Fsd);
    mp_state_durations = cell2mat(mp_state_durations)'; lfp_state_durations = cell2mat(lfp_state_durations)';
    mp_state_durations(:,bad_mp_states) = []; 
    lfp_state_durations(:,bad_lfp_states) = []; 
        
    core_data(d).mp_up_durs = mp_state_durations(2,:);
    core_data(d).mp_down_durs = mp_state_durations(1,:);
    core_data(d).lfp_up_durs = lfp_state_durations(2,:);
     core_data(d).lfp_down_durs = lfp_state_durations(1,:);
           
    core_data(d).mp_dutycycs = core_data(d).mp_up_durs./(core_data(d).mp_up_durs+core_data(d).mp_down_durs);
    core_data(d).lfp_dutycycs = core_data(d).lfp_up_durs./(core_data(d).lfp_up_durs+core_data(d).lfp_down_durs);
        
    %% compute corresponding state transitions and transition lags
    [corresp_lfp_upinds,corresp_lfp_downinds] = find_corresponding_state_transitions_lookback(...
        up_trans_inds_mp,down_trans_inds_mp,up_trans_inds_lfp,down_trans_inds_lfp);
    
    %find non-skipped states
    non_skipped_mp_up_states = find(~isnan(corresp_lfp_upinds));
    non_skipped_mp_down_states = find(~isnan(corresp_lfp_downinds));
    
    %compute transition lags for non-skipped states
    mp_uplags = nan(size(up_trans_inds_mp));
    mp_downlags = nan(size(up_trans_inds_mp));
    mp_uplags(non_skipped_mp_up_states) = up_trans_inds_mp(non_skipped_mp_up_states) - up_trans_inds_lfp(corresp_lfp_upinds(non_skipped_mp_up_states));
    mp_downlags(non_skipped_mp_down_states) = down_trans_inds_mp(non_skipped_mp_down_states) - down_trans_inds_lfp(corresp_lfp_downinds(non_skipped_mp_down_states));
    
    core_data(d).mp_uplags = mp_uplags;
    core_data(d).mp_downlags = mp_downlags;
    
    %compute transition lags relative to individual ctx state durations
    mp_reldownlags = nan(size(up_trans_inds_mp));
    mp_reluplags = nan(size(up_trans_inds_mp));
    for i = 1:length(up_trans_inds_mp)
        if ~isnan(corresp_lfp_downinds(i))
            mp_reldownlags(i) = mp_downlags(i)/core_data(d).lfp_down_durs(corresp_lfp_downinds(i))/Fsd;
        end
        if ~isnan(corresp_lfp_upinds(i))
            mp_reluplags(i) = mp_uplags(i)/core_data(d).lfp_up_durs(corresp_lfp_upinds(i))/Fsd;
        end
    end
    
    core_data(d).mp_reluplags = mp_reluplags;
    core_data(d).mp_reldownlags = mp_reldownlags;
 
    
    %% compute persistence
    [mp_upskipped,mp_downskipped] = greedy_find_skipped_ncx_states_v2(...
        corresp_lfp_upinds,corresp_lfp_downinds,core_data(d).lfp_up_durs,core_data(d).lfp_down_durs,...
        up_trans_inds_lfp,down_trans_inds_lfp,mp_state_vec,thresh_lfp_downdur,thresh_lfp_updur);

    rt2_ups = find(mp_upskipped.rnum_skipped > 0); %indices of MP up states that are robust persistent
    t2_ups = find(mp_upskipped.num_skipped > 0); %indices of MP up states that are persistent
    nrt2_ups = find(mp_upskipped.num_skipped == 0); %MP up states that are non-persistent
    fract_rt2_ups = length(rt2_ups)/(length(rt2_ups) + length(nrt2_ups)); %fraction of up states that are robust persistent (not including ones that are non-robust persistent)
    fract_t2_ups = length(t2_ups)/(length(t2_ups) + length(nrt2_ups)); %fraction of up states that are pers
    
    core_data(d).rt2_ups = rt2_ups;
    core_data(d).t2_ups = t2_ups;
    core_data(d).nrt2_ups = nrt2_ups;
    core_data(d).fract_rt2_ups = fract_rt2_ups;
    core_data(d).fract_t2_ups = fract_t2_ups;
    
    t2_downs = find(mp_downskipped.num_skipped > 0); %indices of MP down states that skip at least one cortical UP state 
    rt2_downs = find(mp_downskipped.rnum_skipped > 0); %indices of MP down states that are robust persistent
    nrt2_downs = find(mp_downskipped.num_skipped==0); %number of down states that didn't skip any LFP up states
    fract_t2_downs = length(t2_downs)/(length(t2_downs) + length(nrt2_downs)); %fraction of MP down states that are pers
    fract_rt2_downs = length(rt2_downs)/(length(rt2_downs) + length(nrt2_downs)); %fraction of (possible) MP down states that are robust persistent
    
    core_data(d).rt2_downs = rt2_downs;
    core_data(d).t2_downs = t2_downs;
    core_data(d).nrt2_downs = nrt2_downs;
    core_data(d).fract_rt2_downs = fract_rt2_downs;
    core_data(d).fract_t2_downs = fract_t2_downs;
  
    robust_non_skipped_mp_ups = [rt2_ups; nrt2_ups]; %set of MP UP states that are clearly either pers or not, and are not skipped by cortical UDS
    robust_non_skipped_mp_downs = [rt2_downs; nrt2_downs]; %set of MP DOWN states that are clearly either pers or not, and not skipped by cortical UDS
        
    %compute duration of each MP state relative to the corresponding LFP
    %state
    mp_rel_updurs = nan(size(up_trans_inds_mp));
    mp_rel_downdurs = nan(size(up_trans_inds_mp));
    mp_rel_updurs(robust_non_skipped_mp_ups) = core_data(d).mp_up_durs(robust_non_skipped_mp_ups) - ...
        core_data(d).lfp_up_durs(corresp_lfp_upinds(robust_non_skipped_mp_ups));
    mp_rel_downdurs(robust_non_skipped_mp_downs) = core_data(d).mp_down_durs(robust_non_skipped_mp_downs) - ...
        core_data(d).lfp_down_durs(corresp_lfp_downinds(robust_non_skipped_mp_downs));
          
    core_data(d).mp_rel_updurs = mp_rel_updurs;
    core_data(d).mp_rel_downdurs = mp_rel_downdurs;
    
    %% compute durations in units of Ncx UDS cycles
    [mp_updurs_lfpc,mp_downdurs_lfpc] = find_duration_ncx_uds_cycles(up_trans_inds_mp,down_trans_inds_mp,...
        mp_uplags,mp_downlags,lfp_period_vec);
    core_data(d).mp_updurs_lfpc = mp_updurs_lfpc;
    core_data(d).mp_downdurs_lfpc = mp_downdurs_lfpc;
    
    %% compile list of skipped cortical states
    skipped_lfp_downs = [];
    for j = 1:length(up_trans_inds_mp)
        skipped_lfp_downs = [skipped_lfp_downs mp_upskipped.inds{j}];
    end
    non_skipped_lfp_downs = setdiff(1:length(down_trans_inds_lfp),skipped_lfp_downs);
    
    %without robust constraint on min state durs
    core_data(d).nr_skipped_lfp_downs = skipped_lfp_downs;
    dore_data(d).nr_non_skipped_lfp_downs = non_skipped_lfp_downs;
    
    %don't count LFP down states that are too short
    too_short_lfp_downs = find(core_data(d).lfp_down_durs < thresh_lfp_downdur);
    skipped_lfp_downs(ismember(skipped_lfp_downs,too_short_lfp_downs)) = [];
    non_skipped_lfp_downs(ismember(non_skipped_lfp_downs,too_short_lfp_downs)) = [];
    
    %these are skipped and non-skipped LFP down states with robustness
    %constraint
    core_data(d).skipped_lfp_downs = skipped_lfp_downs;
    core_data(d).non_skipped_lfp_downs = non_skipped_lfp_downs;
    
    %now repeat for up states
    skipped_lfp_ups = [];
    for j = 1:length(up_trans_inds_mp)
        skipped_lfp_ups = [skipped_lfp_ups mp_downskipped.inds{j}];
    end
    non_skipped_lfp_ups = setdiff(1:length(up_trans_inds_lfp),skipped_lfp_ups);
 
    %without robustness constraint
    core_data(d).nr_skipped_lfp_ups = skipped_lfp_ups;
    core_data(d).nr_non_skipped_lfp_ups = non_skipped_lfp_ups;
    
    %don't count LFP uptrans for up states that are too short
    too_short_lfp_ups = find(core_data(d).lfp_up_durs < thresh_lfp_updur);
    skipped_lfp_ups(ismember(skipped_lfp_ups,too_short_lfp_ups)) = [];
    non_skipped_lfp_ups(ismember(non_skipped_lfp_ups,too_short_lfp_ups)) = [];

    core_data(d).skipped_lfp_ups = skipped_lfp_ups;
    core_data(d).non_skipped_lfp_ups = non_skipped_lfp_ups;
    
end

% cd C:\WC_Germany\final_pdown_analysis\
% save fin_pdown_core_analysis core_data 
cd ~/Analysis/Mayank/final_pdown_analysis/
save fin_pdown_core_analysis_fin core_data 


%%
l3mec = strcmp({data.loc},'MEC');
l3lec = strcmp({data.loc},'LEC');



