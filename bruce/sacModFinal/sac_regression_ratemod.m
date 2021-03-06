% clear all
% close all

addpath('~/James_scripts/CircStat2011f/')
global Expt_name bar_ori

% Expt_name = 'M277';
% bar_ori = 0;

savename = 'sac_glm_data';
include_bursts = 0;
if include_bursts
    savename = [savename '_withbursts'];
end

%%
Expt_num = str2num(Expt_name(2:end));

if Expt_num >= 280
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end

save_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

if Expt_name(1) == 'G'
    n_probes = 96;
    rec_type = 'UA';
    good_coils = [1 0]; %which coils are usable
elseif Expt_name(1) == 'M'
    n_probes = 24;
    rec_type = 'LP';
    switch Expt_num
        case 266
            bar_ori = 80;
        case 270
            bar_ori = 60;
        case 275
            bar_ori = 135;
        case 277
            bar_ori = 70;
        case 281
            bar_ori = 140;
        case 287
            bar_ori = 90;
        case 289
            bar_ori = 160;
        case 294
            bar_ori = 40;
        case 296
            bar_ori = 45;
    end
    good_coils = [1 1]; %which coils are usable
    
else
    error('Unrecognized experiment name');
end

if strcmp(rec_type,'LP')
    if Expt_num >= 275
        rpt_seed = 1001; %M275 M277 M281
    else
        rpt_seed = 1e4; %M270 and 266
    end
else
    rpt_seed = nan;
end

ignore_blocks = [];
switch Expt_num
    case 270
        ignore_blocks = [5 19];
    case 289
        ignore_blocks = [27]; %27 is off somehow
    case 294
        ignore_blocks = [37 38 39]; %37-39 have slightly different dw used in these blocks
    case 86
        ignore_blocks = [16 17 28 30];
    case 87
        ignore_blocks = [15];
    case 93
        ignore_blocks = [28 52];
end

use_coils = [1 0];

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.01;

min_trial_dur = 0.75;
beg_buffer = 0.2;
end_buffer = 0.05;

backlag = round(0.2/dt);
forlag = round(0.4/dt);

sua_sm_sig = (0.01/dt);


if strcmp(Expt_name,'G081') || ismember(Expt_num,[232 235 239])
    trial_dur = 2;
else
    trial_dur = 4;
end

EP_bounds = 1;%eye position boundary
micro_thresh = 1; %microsaccade amplitude threshold (deg)
max_sac_dur = 0.1; %maximum saccade duration (otherwise likely a blink)
sac_burst_isi = 0.15; %minimum inter-saccade interval before classifying sac as part of a 'burst'
max_gsac_dur = 0.1;


%% LOAD EXPTS STRUCT
cd(data_dir)
if Expt_name(1) == 'G'
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
else
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
end
if ~strcmp(Expt_name,'G081') && ~ismember(Expt_num,[232 235 239])
    load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
end

%% PARSE EXPTS STRUCT
is_sim_msac_expt = false;
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa','rls.FaXwi','rls.FaXwiXimi','rls.AllSacB'};
if any(strcmp(Expt_name,{'G095','M266','M270'}))
    is_sim_msac_expt = true;
end
expt_has_simBlanks = false;
if any(strcmp(Expt_name,{'M296','M297'}))
    expt_has_simBlanks = true;
end

%later experiments have fully-interleaved trials
is_TBT_expt = false;
if Expt_num >= 275
    is_TBT_expt = true;
end

%extract block-wise stim parameters
has_data = cellfun(@(x) length(x),Expts) > 0;
expt_names(has_data) = cellfun(@(x) x.Header.expname,Expts(has_data),'UniformOutput',false);
expt_backstim(has_data) = cellfun(@(x) x.Stimvals.Bs,Expts(has_data),'UniformOutput',false);
expt_imback = strcmp(expt_backstim,'image');
expt_dd(has_data) = cellfun(@(x) x.Stimvals.dd,Expts(has_data),'UniformOutput',true);
expt_bar_ori(has_data) = cellfun(@(x) x.Stimvals.or,Expts(has_data),'UniformOutput',true);
expt_Fr(has_data) = cellfun(@(x) x.Stimvals.Fr,Expts(has_data),'UniformOutput',true);
expt_sac_dir(has_data) = mod(cellfun(@(x) x.Stimvals.Fa,Expts(has_data),'UniformOutput',true),180);
expt_sac_amp(has_data) = cellfun(@(x) x.Stimvals.Fs,Expts(has_data),'UniformOutput',true);
expt_ce(has_data) = cellfun(@(x) x.Stimvals.ce,Expts(has_data),'UniformOutput',true);
expt_ijump(has_data) = cellfun(@(x) x.Stimvals.ijump,Expts(has_data),'UniformOutput',true);

included_type(has_data) = ismember(expt_names(has_data),include_expts);

if strcmp(rec_type,'LP')
    expt_bar_ori(expt_bar_ori > 360) = bar_ori;
end

%select blocks for analysis
if strcmp(rec_type,'UA')
    cur_block_set = find(included_type & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_ce == 1);
else
    cur_block_set = find(included_type & expt_ce == 1 & expt_Fr == 1 & expt_ce == 1 & expt_bar_ori == bar_ori);
end
if strcmp(Expt_name,'G081')
    expt_has_ds = (expt_ijump==0)';
end
expt_has_ds(isnan(expt_has_ds)) = 0;
cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];
if length(unique(expt_dd(cur_block_set))) > 1
    fprintf('Warning, multiple dds detected!\n');
    main_dds = mode(expt_dd(cur_block_set));
    cur_block_set(expt_dd(cur_block_set) ~= main_dds) = [];
end

%identify sim-sac imback and grayback blocks
sim_sac_expts = find(expt_has_ds(cur_block_set) ~= 1);
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');
if is_TBT_expt
    sim_sac_expts = []; imback_gs_expts = []; grayback_gs_expts = [];
end

poss_orth_expts = find(mod(expt_bar_ori(cur_block_set) - expt_sac_dir(cur_block_set),180) == 90);
if ~isempty(poss_orth_expts)
    fprintf('Warning, possible orthoganol saccade expts detected\n');
end

if ~isempty(grayback_gs_expts) || ~isempty(imback_gs_expts)
    gsac_amp = unique(expt_sac_amp(cur_block_set([grayback_gs_expts; imback_gs_expts])));
else
    gsac_amp = unique(expt_sac_amp(cur_block_set));
end
if length(gsac_amp) > 1
    fprintf('Multiple guided sac amps detected!\n');
end
%minimum (parallel) amplitude for a guided saccade to be included in
%analysis
gsac_thresh = mean(gsac_amp)/2;

if strcmp(Expt_name,'G081')
    sim_sac_times = [0.7 1.4];
else
    sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
end

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_ttill_end = [];
all_blockvec = [];
all_trialvec = [];
all_trial_wi = [];
all_trial_back = [];
all_trial_Ff = [];
all_trial_Fs = [];
all_trial_imi = [];
all_trial_exvals = [];
all_trial_se = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_trial_blocknums = [];
all_bin_edge_pts = [];

all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:length(cur_block_set);
    cur_block = cur_block_set(ee);
    if ismember(ee,grayback_gs_expts)
        fprintf('Expt %s Block %d of %d; grayback GS, ori:%d\n',Expt_name,ee,length(cur_block_set),expt_bar_ori(cur_block));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Expt %s Block %d of %d; imback GS, ori:%d\n',Expt_name,ee,length(cur_block_set),expt_bar_ori(cur_block));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Expt %s Block %d of %d; SimSac, ori:%d\n',Expt_name,ee,length(cur_block_set),expt_bar_ori(cur_block));
    else
        fprintf('Expt %s Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_name,ee,length(cur_block_set));
    end
    
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    
    [un_ids,id_inds] = unique(trial_ids);
    if length(un_ids) < length(trial_ids)
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    %in this expt there is variable stripe width
    if strcmp(Expt_name,'G093')
        if isfield(Expts{cur_block}.Trials,'wi')
            trial_wi = [Expts{cur_block}.Trials(:).wi];
            trial_wi = trial_wi(id_inds);
        else
            trial_wi = ones(1,length(id_inds))*1.9959;
        end
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    end
    
    %for variable saccade amps
    if isfield(Expts{cur_block}.Trials(1),'Fs')
        trial_Fs = [Expts{cur_block}.Trials(:).Fs];
    else
        trial_Fs = nan(1,length(trial_durs));
    end
    all_trial_Fs = cat(1,all_trial_Fs,trial_Fs(use_trials)');
    
    if isfield(Expts{cur_block}.Trials,'imi')
        trial_imi = [Expts{cur_block}.Trials(:).imi];
        trial_imi = trial_imi(id_inds);
        all_trial_imi = cat(1,all_trial_imi,trial_imi(use_trials)');
    end

    if is_TBT_expt
        if isfield(Expts{cur_block}.Trials,'Bs')
            trial_back = strcmp('image',{Expts{cur_block}.Trials(:).Bs});
        else
            trial_back = nan(1,length(use_trials));
        end
        trial_back = trial_back(id_inds);
        
        if isfield(Expts{cur_block}.Trials,'exvals')
            exvals = reshape([Expts{cur_block}.Trials(:).exvals],length(Expts{cur_block}.Trials(1).exvals),[]);
            trial_exvals = exvals(:,id_inds)';
        else
            trial_exvals = nan(length(id_inds),3);
        end
        
        if isfield(Expts{cur_block}.Trials,'Ff')
            trial_Ff = [Expts{cur_block}.Trials(:).Ff];
            trial_Ff = trial_Ff(id_inds);
        else
            trial_Ff = nan(1,length(id_inds));
        end
        all_trial_back = cat(1,all_trial_back,trial_back(use_trials)');
        all_trial_Ff = cat(1,all_trial_Ff,trial_Ff(use_trials)');
        all_trial_exvals = cat(1,all_trial_exvals,trial_exvals(use_trials,:));
    end
    
    if isfield(Expts{cur_block}.Trials(1),'se')
        trial_se = [Expts{cur_block}.Trials(:).se];
        trial_se = trial_se(id_inds);
    else
        trial_se = nan(size(id_inds));
    end
    all_trial_se = cat(1,all_trial_se,trial_se(use_trials)');
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        cur_ttill_end = trial_end_times(use_trials(tt)) - cur_t_axis;
        
        all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
        all_tsince_start = [all_tsince_start; cur_tsince_start];
        all_ttill_end = [all_ttill_end; cur_ttill_end];
        all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    end
    if n_trials > 0
        trial_cnt = trial_cnt + n_trials;
        if strcmp(rec_type,'LP')
            trial_toffset(ee) = all_t_bin_edges(end);
            cur_toffset = trial_toffset(ee);
            cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
        end
    end
end

if is_TBT_expt
    if any(~isnan(all_trial_exvals(:,3)))
        fprintf('Using exvals to define trial-by-trial conditions\n');
        all_trial_Ff(all_trial_exvals(:,3) == 1) = 0; %these are sim sac trials
        all_trial_Ff(all_trial_exvals(:,3) > 1) = 70;
        all_trial_back(all_trial_exvals(:,3) == 2) = 0; %these are gray-back trials
        all_trial_back(ismember(all_trial_exvals(:,3),[1 3])) = 1;
    end
end

trial_start_inds = [1; find(diff(all_trialvec) > 0)+1];

%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = n_probes;
if strcmp(rec_type,'LP')
    clust_params.exclude_adjacent = true;
else
clust_params.exclude_adjacent = false;
end
[all_binned_mua,all_binned_sua,Clust_data] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;
SU_isodist = Clust_data.SU_isodists;
SU_Lratio = Clust_data.SU_Lratios;

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & all_ttill_end >= end_buffer);

%for G093 use only data where stripe width is 2 deg
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi == un_wi_vals(2));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end

%get rid of repeat trials
rpt_trials = find(all_trial_se == rpt_seed);
used_inds(ismember(all_trialvec(used_inds),rpt_trials)) = [];

%% PROCESS EYE TRACKING DATA
if isfield(Expts{cur_block_set(1)}.Header,'exptno')
    em_block_nums = cellfun(@(X) X.Header.exptno,Expts(cur_block_set),'uniformoutput',1); %block numbering for EM/LFP data sometimes isnt aligned with Expts struct
else
    em_block_nums = cur_block_set;
end

[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,em_block_nums,Expt_name,trial_toffset,good_coils);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data_v2(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

[saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params);

is_blink = detect_blinks(all_eye_ts,all_eye_vals,saccades,et_params);

[saccades,is_blink] = merge_blinks(saccades,is_blink);

sac_start_times = [saccades(:).start_time];
sac_stop_times = [saccades(:).stop_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
sac_stop_times(bad_sacs) = [];
saccades(bad_sacs) = [];
is_blink(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

%% PROCESS SACCADE STATS

%interpolate saccade start times and get rid of saccades that aren't within
%the t-axis binning
sac_start_times = [saccades(:).start_time];
sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
sac_start_inds(isnan(sac_start_inds)) = 1;
sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).stop_time]));
sac_error = abs(sac_start_times - all_t_axis(sac_start_inds)');
usac_set = find(ismember(sac_start_inds,used_inds) & ismember(sac_stop_inds,used_inds) & sac_error <= dt & sac_start_times >= all_t_axis(1) & [saccades(:).stop_time] <= all_t_axis(end));
saccades = saccades(usac_set);
used_is_blink = is_blink(usac_set);

sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).start_time]));
sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).stop_time]));
sac_peak_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).peak_time]));

%start and stop positions of each saccade
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));

%find saccades that start or end out of window
out_bounds = abs(sac_prepos(2,:)) > EP_bounds | abs(sac_postpos(2,:)) > EP_bounds;

sacburst_set = find([saccades(:).isi] < sac_burst_isi | [saccades(:).next_isi] < sac_burst_isi);
micro_set = find([saccades(:).amplitude] < micro_thresh & ~used_is_blink' & ~out_bounds);
msac_bursts = micro_set(ismember(micro_set,sacburst_set));
if ~include_bursts
micro_set(ismember(micro_set,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'
end

%saccade amplitude along parallel axis
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

%guided saccades are those whose parallel component is large enough and
%that aren't blinks
sac_durs = [saccades(:).duration];
gsac_set = find(abs(sac_deltaX) > gsac_thresh & ~used_is_blink' & ~out_bounds & sac_durs <= max_gsac_dur);

%compile indices of simulated saccades
all_sim_sacs = [];
all_sim_msacs = [];
if is_sim_msac_expt  %if there are simulated microsaccades in this dataset
    big_simsac_trials = find(all_trial_Fs > 1);
    small_simsac_trials = find(all_trial_Fs < 1);
else
    big_simsac_trials = 1:length(all_trial_start_times);
    small_simsac_trials = [];
end
if is_TBT_expt
    sim_sac_trials = find(all_trial_Ff == 0);
    sim_trial_inds = find(ismember(all_trialvec,sim_sac_trials));
    sim_sacs = cell(length(sim_sac_times),1);
    for ii = 1:length(sim_sac_times)
        sim_sacs{ii} = sim_trial_inds(all_tsince_start(sim_trial_inds(1:end-1)) < sim_sac_times(ii) & ...
            all_tsince_start(sim_trial_inds(2:end)) >= sim_sac_times(ii));
        cur_big_set = find(ismember(all_trialvec(sim_sacs{ii}),big_simsac_trials));
        cur_micro_set = find(ismember(all_trialvec(sim_sacs{ii}),small_simsac_trials));
        all_sim_sacs = [all_sim_sacs; sim_sacs{ii}(cur_big_set)];
        all_sim_msacs = [all_sim_msacs; sim_sacs{ii}(cur_micro_set)];
    end
else
    sim_expt_inds = find(ismember(all_blockvec,sim_sac_expts));
    sim_sacs = cell(length(sim_sac_times),1);
    for ii = 1:length(sim_sac_times)
        sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
            all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
        cur_big_set = find(ismember(all_trialvec(sim_sacs{ii}),big_simsac_trials));
        cur_micro_set = find(ismember(all_trialvec(sim_sacs{ii}),small_simsac_trials));
        all_sim_sacs = [all_sim_sacs; sim_sacs{ii}(cur_big_set)];
        all_sim_msacs = [all_sim_msacs; sim_sacs{ii}(cur_micro_set)];
    end
end
all_sim_sacs = sort(all_sim_sacs(ismember(all_sim_sacs,used_inds)));
all_sim_msacs = sort(all_sim_msacs(ismember(all_sim_msacs,used_inds)));

%for TBT expts, determine which saccades were part of image-back trials vs
%gray-back trials
sac_trial_inds = all_trialvec(sac_start_inds);
if is_TBT_expt
    grayback_gs_trials = find(all_trial_back == 0);
    imback_trials = find(all_trial_back  == 1);
    imback_gs_trials = imback_trials;
    imback_gs_trials(ismember(imback_gs_trials,sim_sac_trials)) = [];
    
    gback_sacs = find(ismember(sac_trial_inds,grayback_gs_trials));
    iback_sacs = find(ismember(sac_trial_inds,imback_trials));
    ss_sacs = find(ismember(sac_trial_inds,sim_sac_trials));
else
    gback_sacs = find(ismember(all_blockvec(sac_start_inds),grayback_gs_expts));
    iback_sacs = find(ismember(all_blockvec(sac_start_inds),imback_gs_expts));
    ss_sacs = find(ismember(all_blockvec(sac_start_inds),sim_sac_expts));
end

gray_msac_set = intersect(gback_sacs,micro_set);
im_msac_set = intersect(iback_sacs,micro_set);

%% SMOOTH AND NORMALIZE SPIKING DATA
all_sua_rate = nan(size(all_binned_sua));
sua_block_mean_rates = nan(length(cur_block_set),length(SU_numbers));
sua_block_n_spikes = nan(length(cur_block_set),length(SU_numbers));
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
        for ss = 1:length(SU_numbers)
            if sua_sm_sig > 0
                all_sua_rate(cur_block_inds,ss) = jmm_smooth_1d_cor(all_binned_sua(cur_block_inds,ss),sua_sm_sig);
            else
                all_sua_rate(cur_block_inds,ss) = all_binned_sua(cur_block_inds,ss);
            end
        end
        sua_block_mean_rates(ee,:) = mean(all_sua_rate(cur_block_inds,:));
        sua_block_n_spikes(ee,:) = sum(all_binned_sua(cur_block_inds,:));
    end
end

%normalized firing rates within each block
all_sua_rate_norm = nan(size(all_sua_rate));
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
        all_sua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_sua_rate(cur_block_inds,:),...
            sua_block_mean_rates(ee,:));
    end
end

if is_TBT_expt
    grayback_trial_inds = find(ismember(all_trialvec,grayback_gs_trials));
    all_sua_rate_norm(grayback_trial_inds,:) = bsxfun(@rdivide,all_sua_rate_norm(grayback_trial_inds,:),nanmean(all_sua_rate_norm(grayback_trial_inds,:)));
    
    imback_trial_inds = find(ismember(all_trialvec,imback_gs_trials));
    all_sua_rate_norm(imback_trial_inds,:) = bsxfun(@rdivide,all_sua_rate_norm(imback_trial_inds,:),nanmean(all_sua_rate_norm(imback_trial_inds,:)));
    
    simsac_trial_inds = find(ismember(all_trialvec,sim_sac_trials));
    all_sua_rate_norm(simsac_trial_inds,:) = bsxfun(@rdivide,all_sua_rate_norm(simsac_trial_inds,:),nanmean(all_sua_rate_norm(simsac_trial_inds,:)));
end


%%
n_blocks = length(unique(all_blockvec));
Xblock = zeros(length(all_t_axis),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%%
sim_sac_trials = all_trialvec(all_sim_sacs);

slags = -backlag:forlag;
n_sac_bins = length(slags);
Xsac = zeros(length(all_t_axis),length(slags));
Xmsac = zeros(length(all_t_axis),length(slags));
Xsimsac = zeros(length(all_t_axis),length(slags));
for ii = 1:n_sac_bins
    cur_sac_target = sac_start_inds(gsac_set) + slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < length(all_t_axis));
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(cur_sac_target) ~= sac_trial_inds(gsac_set(uu))) = [];
    Xsac(cur_sac_target,ii) = 1;
    
    cur_sac_target = sac_start_inds(micro_set) + slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < length(all_t_axis));
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(cur_sac_target) ~= sac_trial_inds(micro_set(uu))) = [];
    Xmsac(cur_sac_target,ii) = 1;

    cur_sac_target = all_sim_sacs + slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < length(all_t_axis));
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(cur_sac_target) ~= sim_sac_trials(uu)) = [];
    Xsimsac(cur_sac_target,ii) = 1;
end

%%
cd(save_dir)
sname = 'sac_trig_avg_data';
if include_bursts
    sname = [sname '_withburst'];
end
sname = [sname sprintf('_ori%d',bar_ori)];
load(sname);

tlags = trig_avg_params.lags;
tdt = trig_avg_params.dt;
tlags = tlags*tdt;

%% COMPUTE CROSS CORRELATION BETWEEN MSAC TIMES WITH GSACS AND SIMSACS10
binned_msacs = Xmsac(:,slags==0);
binned_gsacs = Xsac(:,slags==0);
binned_simsacs = Xsimsac(:,slags==0);

non_simsac_trials = setdiff(1:length(unique(all_trialvec)),sim_sac_trials);

maxlag = round(1/dt);
cur_uinds = find(ismember(all_trialvec(used_inds),sim_sac_trials));
% [msac_simsac_xcorr,xcorr_lags] = xcov(binned_msacs,binned_simsacs,maxlag,'coeff');
[msac_simsac_xcorr,xcorr_lags] = get_event_trig_avg_v3(binned_msacs(cur_uinds),find(binned_simsacs(cur_uinds)==1),maxlag,maxlag,all_trialvec(used_inds(cur_uinds)));

cur_uinds = find(ismember(all_trialvec(used_inds),non_simsac_trials));
% [msac_gsac_xcorr,xcorr_lags] = xcov(binned_msacs,binned_gsacs,maxlag,'coeff');
[msac_gsac_xcorr,xcorr_lags] = get_event_trig_avg_v3(binned_msacs(cur_uinds),find(binned_gsacs(cur_uinds)==1),maxlag,maxlag,all_trialvec(used_inds(cur_uinds)));

cur_uinds = 1:length(used_inds);
[msac_msac_xcorr,xcorr_lags] = get_event_trig_avg_v3(binned_msacs(cur_uinds),find(binned_msacs(cur_uinds)==1),maxlag,maxlag,all_trialvec(used_inds(cur_uinds)));
msac_msac_xcorr(xcorr_lags == 0) = 0;

xcorr_data.msac_simsac_xcorr = msac_simsac_xcorr;
xcorr_data.msac_gsac_xcorr = msac_gsac_xcorr;
xcorr_data.msac_msac_xcorr = msac_msac_xcorr;
xcorr_data.lags = xcorr_lags*dt;

xcorr_data.ov_msac_rate = sum(Xmsac(:,slags==0))/size(Xmsac,1)/dt;
%%
lambda_d2T = 50;
lambda_L2 = 5;
min_Nsacs = 100;
silent = 1;

stim_params(1) = NMMcreate_stim_params(n_blocks);
stim_params(2:4) = NMMcreate_stim_params(length(slags));
reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);

X{1} = Xblock; %block ID predictor
X{2} = Xsac;
X{3} = Xmsac;
X{4} = Xsimsac;

for ss = 1:length(SU_numbers)
    fprintf('Computing GLMs for SU %d of %d\n',ss,length(SU_numbers));
    
    cc_uinds = used_inds(~isnan(all_binned_sua(used_inds,ss)));
    cur_Robs = all_binned_sua(cc_uinds,ss);
    if ~isempty(cc_uinds)
        
        %construct model structure based on what kinds of saccades we have
        %sufficient data for
        mod_signs = [1];
        mod_Xtargs = [1];
        if sua_data(ss).N_gsacs >= min_Nsacs %if enough guided sacs
            mod_signs = [mod_signs 1];
            mod_Xtargs = [mod_Xtargs 2];
        end
        if sua_data(ss).N_msacs >= min_Nsacs %if enough microsacs
            mod_signs = [mod_signs 1];
            mod_Xtargs = [mod_Xtargs 3];
        end 
        if sua_data(ss).N_simsacs >= min_Nsacs %if enough sim sacs
            mod_signs = [mod_signs 1];
            mod_Xtargs = [mod_Xtargs 4];
        end
        
        glm = NMMinitialize_model(stim_params,mod_signs,repmat({'lin'},1,length(mod_signs)),reg_params,mod_Xtargs);
        glm.mods(1).reg_params = NMMcreate_reg_params();
        
        glm = NMMfit_filters(glm,cur_Robs,get_Xcell_tInds(X,cc_uinds),[],[],silent);
        glm = NMMfit_logexp_spkNL( glm, cur_Robs, get_Xcell_tInds(X,cc_uinds));
        
        [LL, penLL, pred_rate, G, gint, fgint] = NMMmodel_eval(glm, cur_Robs, get_Xcell_tInds(X,cc_uinds));
        avg_mod_outs = nanmean(fgint); %avg output of each subunit
        temp_kerns = [glm.mods(2:end).filtK]; %these are the saccade kernels
        
        sua_data(ss).glm_ovavg_rate = mean(cur_Robs);
        sua_data(ss).glm_lags = slags;
        sua_data(ss).glm_dt = dt;
                
        gsac_kern = find(mod_Xtargs(2:end) == 2);
        if ~isempty(gsac_kern)
            block_arate = nan(n_blocks,1); 
            block_reldurs = nan(n_blocks,1);
            block_rate_out = nan(n_blocks,length(slags));
            for bb = 1:n_blocks
                cur_set = find(all_blockvec(cc_uinds) == bb);
                cur_set(~ismember(all_trialvec(cc_uinds(cur_set)),non_simsac_trials)) = []; %exclude sim sac trials (important for TBT recs)
                block_arate(bb) = nanmean(cur_Robs(cur_set));
                block_reldurs(bb) = length(cur_set);
                
                %compute dependence of rate on gsac while holding other
                %predictors at their avg values
                other_kerns = find(mod_Xtargs ~= 2 & mod_Xtargs ~= 1);
                g_out = temp_kerns(:,gsac_kern) + sum(avg_mod_outs(other_kerns)) + glm.mods(1).filtK(bb);
                exp_g = (g_out + glm.spk_NL_params(1))*glm.spk_NL_params(2);
                block_rate_out(bb,:) = glm.spk_NL_params(3)*log(1+exp(exp_g));
            end
            block_reldurs = block_reldurs/nansum(block_reldurs);
            cur_norm_rate = bsxfun(@rdivide,block_rate_out,block_arate); %normalize within each block
            sua_data(ss).glm_gsac_rate = nansum(bsxfun(@times,cur_norm_rate,block_reldurs)); %compute weighted avg across blocks
            
%             trig_avg = get_event_trig_avg_v3(cur_Robs,find(Xsac(cc_uinds,slags==0)==1),backlag,forlag,[],all_trialvec(cc_uinds));
%             sua_data(ss).tavg_gsac_rate = jmm_smooth_1d_cor(trig_avg,0.01/dt);
            trig_avg = get_event_trig_avg_v3(all_sua_rate_norm(cc_uinds,ss),find(Xsac(cc_uinds,slags==0)==1),backlag,forlag,[],all_trialvec(cc_uinds));
            sua_data(ss).tavg_gsac_rate = trig_avg;
        else
            sua_data(ss).glm_gsac_rate = nan;
            sua_data(ss).tavg_gsac_rate = nan;
        end
        any_gsac = find(any(Xsac(cc_uinds,:) > 0,2));
        sua_data(ss).glm_gsac_avgrate = mean(cur_Robs(any_gsac));
        
        msac_kern = find(mod_Xtargs(2:end) == 3);
        if ~isempty(msac_kern)            
            block_arate = nan(n_blocks,1);
            block_reldurs = nan(n_blocks,1);
            block_rate_out = nan(n_blocks,length(slags));
            for bb = 1:n_blocks
                cur_set = find(all_blockvec(cc_uinds) == bb);
                block_arate(bb) = nanmean(cur_Robs(cur_set));
                block_reldurs(bb) = length(cur_set);
                
                other_kerns = find(mod_Xtargs ~= 3 & mod_Xtargs ~= 1);
                g_out = temp_kerns(:,msac_kern) + sum(avg_mod_outs(other_kerns)) + glm.mods(1).filtK(bb);
                exp_g = (g_out + glm.spk_NL_params(1))*glm.spk_NL_params(2);
                block_rate_out(bb,:) = glm.spk_NL_params(3)*log(1+exp(exp_g));
            end
            block_reldurs = block_reldurs/nansum(block_reldurs);
            cur_norm_rate = bsxfun(@rdivide,block_rate_out,block_arate); %normalize within each block
            sua_data(ss).glm_msac_rate = nansum(bsxfun(@times,cur_norm_rate,block_reldurs)); %compute weighted avg across blocks
            
%             trig_avg = get_event_trig_avg_v3(cur_Robs,find(Xmsac(cc_uinds,slags==0)==1),backlag,forlag,[],all_trialvec(cc_uinds));
%             sua_data(ss).tavg_msac_rate = jmm_smooth_1d_cor(trig_avg,0.01/dt);
            trig_avg = get_event_trig_avg_v3(all_sua_rate_norm(cc_uinds,ss),find(Xmsac(cc_uinds,slags==0)==1),backlag,forlag,[],all_trialvec(cc_uinds));
            sua_data(ss).tavg_msac_rate = trig_avg;
        else
            sua_data(ss).glm_msac_rate = nan;
            sua_data(ss).tavg_msac_rate = nan;
        end
        any_msac = find(any(Xmsac(cc_uinds,:) > 0,2));
        sua_data(ss).glm_msac_avgrate = mean(cur_Robs(any_msac));
        
        simsac_kern = find(mod_Xtargs(2:end) == 4);
        if ~isempty(simsac_kern)
            block_arate = nan(n_blocks,1);
            block_reldurs = nan(n_blocks,1);
            block_rate_out = nan(n_blocks,length(slags));
            for bb = 1:n_blocks
                cur_set = find(all_blockvec(cc_uinds) == bb);
                cur_set(~ismember(all_trialvec(cc_uinds(cur_set)),sim_sac_trials)) = []; %exclude guided sac trials (important for TBT recs)
                block_arate(bb) = nanmean(cur_Robs(cur_set));
                block_reldurs(bb) = length(cur_set);
                
                other_kerns = find(mod_Xtargs ~= 4 & mod_Xtargs ~= 1);
                g_out = temp_kerns(:,simsac_kern) + sum(avg_mod_outs(other_kerns)) + glm.mods(1).filtK(bb);
                exp_g = (g_out + glm.spk_NL_params(1))*glm.spk_NL_params(2);
                block_rate_out(bb,:) = glm.spk_NL_params(3)*log(1+exp(exp_g));
            end
            block_reldurs = block_reldurs/nansum(block_reldurs);
            cur_norm_rate = bsxfun(@rdivide,block_rate_out,block_arate); %normalize within each block
            sua_data(ss).glm_simsac_rate = nansum(bsxfun(@times,cur_norm_rate,block_reldurs)); %compute weighted avg across blocks
            
%             trig_avg = get_event_trig_avg_v3(cur_Robs,find(Xsimsac(cc_uinds,slags==0)==1),backlag,forlag,[],all_trialvec(cc_uinds));
%             sua_data(ss).tavg_simsac_rate = jmm_smooth_1d_cor(trig_avg,0.01/dt);
            trig_avg = get_event_trig_avg_v3(all_sua_rate_norm(cc_uinds,ss),find(Xsimsac(cc_uinds,slags==0)==1),backlag,forlag,[],all_trialvec(cc_uinds));
            sua_data(ss).tavg_simsac_rate = trig_avg;
        else
            sua_data(ss).glm_simsac_rate = nan;
            sua_data(ss).tavg_simsac_rate = nan;
        end
        any_simsac = find(any(Xsimsac(cc_uinds,:) > 0,2));
        sua_data(ss).glm_simsac_avgrate = mean(cur_Robs(any_simsac));
        
    end
end

%%
cd(save_dir)
savename = [savename sprintf('_ori%d',bar_ori)];
save(savename,'sua_data','xcorr_data','trig_avg_params');
