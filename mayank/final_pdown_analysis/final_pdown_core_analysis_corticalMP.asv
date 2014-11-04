
clear all
load C:/WC_Germany/final_pdown_analysis/compiled_corticalMP_data.mat

addpath('C:/WC_Germany/persistent_9_27_2010/');

min_rec_dur = 500; %in sec
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);
data = data(used_dirs);

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
rate_sm = round(Fsd*0.05);

maxlag = round(4*Fsd);
backlag = 4*Fsd;
forwardlag = 4*Fsd;

params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 25;

%robust persistence requires MP skips states longer than this
thresh_lfp_updur = 0.5;
thresh_lfp_downdur = 0.5;

%for state duration distribution calculations
up_range = [0.1 15];
down_range = [0.1 15];
numBins = 60;
log_dur_grid = logspace(log10(up_range(1)),log10(up_range(2)),numBins+1);

%%
for d = 1:length(data)
    cd(data(d).dir)
    pwd
    
    core_data(d).data_id = data(d).id;
    
    load ./used_data lf7 wcv_minus_spike
    [lfp_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    lfp_hf = get_hf_features(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    
    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    
    %cut data at the soonest of ep or dp (when data becomes unstable or
    %some other experimental manipulation is applied
    end_time = min(data(d).ep,data(d).dp);
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        synct_d(ep+1:end) = []; lfp_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; lfp_hf(ep+1:end) = []; %hpc_hf(ep+1:end) = [];
    else
        ep = length(t_axis);
    end
    
    load ./pa_hsmm_state_seq_combined_fin_nd.mat
    load ./pa_hsmm_state_seq7_combined_fin_nd.mat
    
    hsmm8 = hsmm7;
    lfp_state_seq = hsmm_bbstate_seq7;
    mp_state_seq = hsmm_bbstate_seq;
    
    %resample UDS sequences and find state transition times
    [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,mp_state_seq);
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    %eliminate state transitions after the data end point
    bad_mp_states = find(up_trans_inds > ep | down_trans_inds > ep);
    up_trans_inds(bad_mp_states) = []; down_trans_inds(bad_mp_states) = [];
    bad_lfp_states = find(up_trans_inds8 > ep | down_trans_inds8 > ep);
    up_trans_inds8(bad_lfp_states) = []; down_trans_inds8(bad_lfp_states) = [];
        
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
    mp_state_number(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;
    mp_state_vec(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;

    %% compute state duration distributions
    [mp_state_durations] = compute_state_durations_seg(mp_state_seq,Fsd);
    [lfp_state_durations] = compute_state_durations_seg(lfp_state_seq,Fsd);
    mp_state_durations(:,bad_mp_states) = []; 
    lfp_state_durations(:,bad_lfp_states) = []; 
        
    core_data(d).mp_up_durs = mp_state_durations(2,:);
    core_data(d).mp_down_durs = mp_state_durations(1,:);
    core_data(d).lfp_up_durs = lfp_state_durations(2,:);
     core_data(d).lfp_down_durs = lfp_state_durations(1,:);
       
    mp_updurs = mp_state_durations(2,:);
    mp_downdurs = mp_state_durations(1,:);
    lfp_updurs = lfp_state_durations(2,:);
    lfp_downdurs = lfp_state_durations(1,:);
    
    core_data(d).mp_dutycycs = mp_updurs./(mp_updurs+mp_downdurs);
    core_data(d).lfp_dutycycs = lfp_updurs./(lfp_updurs+lfp_downdurs);
        
    %% compute corresponding state transitions and transition lags
    [corresp_lfp_upinds,corresp_lfp_downinds] = find_corresponding_state_transitions_lookback(...
        up_trans_inds,down_trans_inds,up_trans_inds8,down_trans_inds8);
    
    %find non-skipped states
    non_skipped_mp_up_states = find(~isnan(corresp_lfp_upinds));
    non_skipped_mp_down_states = find(~isnan(corresp_lfp_downinds));
    
    %compute transition lags for non-skipped states
    mp_uplags = nan(size(up_trans_inds));
    mp_downlags = nan(size(up_trans_inds));
    mp_uplags(non_skipped_mp_up_states) = up_trans_inds(non_skipped_mp_up_states) - up_trans_inds8(corresp_lfp_upinds(non_skipped_mp_up_states));
    mp_downlags(non_skipped_mp_down_states) = down_trans_inds(non_skipped_mp_down_states) - down_trans_inds8(corresp_lfp_downinds(non_skipped_mp_down_states));
    
    core_data(d).mp_uplags = mp_uplags;
    core_data(d).mp_downlags = mp_downlags;
    
    %compute transition lags relative to ctx state durations
    mp_reldownlags = nan(size(up_trans_inds));
    mp_reluplags = nan(size(up_trans_inds));
    for i = 1:length(up_trans_inds)
        if ~isnan(corresp_lfp_downinds(i))
            mp_reldownlags(i) = mp_downlags(i)/lfp_downdurs(corresp_lfp_downinds(i))/Fsd;
        end
        if ~isnan(corresp_lfp_upinds(i))
            mp_reluplags(i) = mp_uplags(i)/lfp_updurs(corresp_lfp_upinds(i))/Fsd;
        end
    end
    
    core_data(d).mp_reluplags = mp_reluplags;
    core_data(d).mp_reldownlags = mp_reldownlags;
 
    
    %% compute persistence
    [mp_upskipped,mp_downskipped] = greedy_find_skipped_ncx_states_v2(...
        corresp_lfp_upinds,corresp_lfp_downinds,lfp_state_durations(2,:),lfp_state_durations(1,:),...
        up_trans_inds8,down_trans_inds8,mp_state_vec,thresh_lfp_downdur,thresh_lfp_updur);

    rt2_ups = find(mp_upskipped.rnum_skipped > 0);
    t2_ups = find(mp_upskipped.num_skipped > 0);
    nrt2_ups = find(mp_upskipped.num_skipped == 0);
    fract_rt2_ups = length(rt2_ups)/(length(rt2_ups) + length(nrt2_ups));
    fract_t2_ups = length(t2_ups)/(length(t2_ups) + length(nrt2_ups));
    
    core_data(d).rt2_ups = rt2_ups;
    core_data(d).t2_ups = t2_ups;
    core_data(d).nrt2_ups = nrt2_ups;
    core_data(d).fract_rt2_ups = fract_rt2_ups;
    core_data(d).fract_t2_ups = fract_t2_ups;
    
    t2_downs = find(mp_downskipped.num_skipped > 0);
    rt2_downs = find(mp_downskipped.rnum_skipped > 0);
    nrt2_downs = find(mp_downskipped.num_skipped==0); %number of down states that didn't skip any LFP up states
    fract_t2_downs = length(t2_downs)/(length(t2_downs) + length(nrt2_downs));
    fract_rt2_downs = length(rt2_downs)/(length(rt2_downs) + length(nrt2_downs));
    
    core_data(d).rt2_downs = rt2_downs;
    core_data(d).t2_downs = t2_downs;
    core_data(d).nrt2_downs = nrt2_downs;
    core_data(d).fract_rt2_downs = fract_rt2_downs;
    core_data(d).fract_t2_downs = fract_t2_downs;
  
    robust_non_skipped_mp_ups = [rt2_ups; nrt2_ups];
    robust_non_skipped_mp_downs = [rt2_downs; nrt2_downs];
        
    mp_rel_updurs = nan(size(up_trans_inds));
    mp_rel_downdurs = nan(size(up_trans_inds));
    mp_rel_updurs(robust_non_skipped_mp_ups) = mp_state_durations(2,robust_non_skipped_mp_ups) - ...
        lfp_state_durations(2,corresp_lfp_upinds(robust_non_skipped_mp_ups));
    mp_rel_downdurs(robust_non_skipped_mp_downs) = mp_state_durations(1,robust_non_skipped_mp_downs) - ...
        lfp_state_durations(1,corresp_lfp_downinds(robust_non_skipped_mp_downs));
          
    core_data(d).mp_rel_updurs = mp_rel_updurs;
    core_data(d).mp_rel_downdurs = mp_rel_downdurs;
    
    %% compute durations in units of Ncx UDS cycles
    [mp_updurs_lfpc,mp_downdurs_lfpc] = find_duration_ncx_uds_cycles(up_trans_inds,down_trans_inds,...
        mp_uplags,mp_downlags,lfp_period_vec);
    core_data(d).mp_updurs_lfpc = mp_updurs_lfpc;
    core_data(d).mp_downdurs_lfpc = mp_downdurs_lfpc;
    
    %% compile list of skipped cortical states
    skipped_lfp_downs = [];
    for j = 1:length(up_trans_inds)
        skipped_lfp_downs = [skipped_lfp_downs mp_upskipped.inds{j}];
    end
    non_skipped_lfp_downs = setdiff(1:length(down_trans_inds8),skipped_lfp_downs);
    
    %don't count LFP down states that are too short
    too_short_lfp_downs = find(lfp_downdurs < thresh_lfp_downdur);
    skipped_lfp_downs(ismember(skipped_lfp_downs,too_short_lfp_downs)) = [];
    non_skipped_lfp_downs(ismember(non_skipped_lfp_downs,too_short_lfp_downs)) = [];
    
    core_data(d).skipped_lfp_downs = skipped_lfp_downs;
    core_data(d).non_skipped_lfp_downs = non_skipped_lfp_downs;
    
    skipped_lfp_ups = [];
    for j = 1:length(up_trans_inds)
        skipped_lfp_ups = [skipped_lfp_ups mp_downskipped.inds{j}];
    end
    non_skipped_lfp_ups = setdiff(1:length(up_trans_inds8),skipped_lfp_ups);
 
    %don't count LFP uptrans for up states that are too short
    too_short_lfp_ups = find(lfp_updurs < thresh_lfp_updur);
    skipped_lfp_ups(ismember(skipped_lfp_ups,too_short_lfp_ups)) = [];
    non_skipped_lfp_ups(ismember(non_skipped_lfp_ups,too_short_lfp_ups)) = [];

    core_data(d).skipped_lfp_ups = skipped_lfp_ups;
    core_data(d).non_skipped_lfp_ups = non_skipped_lfp_ups;
    
end

cd C:\WC_Germany\final_pdown_analysis\
save fin_pdown_core_analysis_corticalMP core_data 


%%
load ./compiled_corticalMP_data
min_rec_dur = 500; %in sec
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);
data = data(used_dirs);

load ./fin_pdown_core_analysis_corticalMP.mat

exclude_set = [6]; %very strong anesthesia causes clear problems with UDS detection
exclude_set = [exclude_set 12 13 14]; %interneurons
exclude = find(ismember([data(:).id],exclude_set));
data(exclude) = [];
core_data(exclude) = [];

ctx_fract_rt2_ups = [core_data(:).fract_rt2_ups];
ctx_fract_rt2_downs = [core_data(:).fract_rt2_downs];
frontal = find(strcmp({data(:).region},'frontal'));
prefrontal = find(strcmp({data(:).region},'prefrontal'));
parietal = find(strcmp({data(:).region},'parietal'));
barrel = find(strcmp({data(:).region},'barrel'));

group = nan(length(ctx_fract_rt2_ups),1);
group(prefrontal) = 3;
group(frontal) = 4;
group(parietal) = 5;

%%
load compiled_data
min_rec_dur = 500; %in sec
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);
data_ids = [data(:).id];
used_dirs(ismember(data_ids(used_dirs),unclear_uds)) = [];
data = data(used_dirs);
data_ids = [data(:).id];

load('fin_pdown_core_analysis.mat')
l3mec = find(strcmp({data.loc},'MEC'));
l3mec(~ismember(data_ids(l3mec),clear_l3)) = [];
l3lec = find(strcmp({data.loc},'LEC'));
fract_rt2_ups = [core_data(:).fract_rt2_ups];
fract_rt2_downs = [core_data(:).fract_rt2_downs];

newgroup = nan(length(fract_rt2_ups),1);
newgroup(l3mec) = 1;
newgroup(l3lec) = 2;

all_fract_ups = [ctx_fract_rt2_ups fract_rt2_ups];
all_fract_downs = [ctx_fract_rt2_downs fract_rt2_downs];
all_group = [group; newgroup];

%%
uset = find(~isnan(all_group));

close all
names = {'L3-MEC','L3-LEC','Prefrontal','Frontal','Parietal'};
h1 = figure;
% boxplot(all_fract_ups(uset)',names(all_group(uset)),'plotstyle','compact');
boxplot(all_fract_ups(uset)',names(all_group(uset)));
ylabel('Fraction persistent UP states');
yl = ylim();
% ylim([0 yl(2)]);
ylim([0 0.6])

h2=figure;
% boxplot(all_fract_downs(uset)',names(all_group(uset)),'plotstyle','compact');
boxplot(all_fract_downs(uset)',names(all_group(uset)));
ylabel('Fraction persistent DOWN states');
yl = ylim();
% ylim([0 yl(2)]);
ylim([0 0.6])

% figure
% plot(all_group(uset),all_fract_ups(uset),'o')

%%
fig_dir = 'C:\WC_Germany\final_pdown_analysis\figures\';

fig_width = 3.5;
rel_heigh = 0.8;

figufy(h1);
fname = [fig_dir 'allcell_pUp_dist2.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

figufy(h2);
fname = [fig_dir 'allcell_pDown_dist2.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
close(h2);

