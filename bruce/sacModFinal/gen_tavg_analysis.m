clear all
% close all

addpath('~/James_scripts/CircStat2011f/')
global Expt_name bar_ori

Expt_name = 'G093';
bar_ori = 0;

include_bursts = 0;

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
dt = 0.005;

min_trial_dur = 0.75;
beg_buffer = 0.2;
end_buffer = 0.05;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
% mua_sm_sig = (0.005/dt);
mua_sm_sig = 0;

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
    
    %     fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    %     load(fname,'Clusters');
    %     for cc = 1:n_probes
    %         all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
    %         all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
    %         all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    %     end
    
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
trial_stop_inds = [find(diff(all_trialvec) > 0); length(all_trialvec)];

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

%if there are simulated blanks in this dataset, find their start indices
all_sim_blanks = [];
if expt_has_simBlanks
    blank_trials = find(all_trial_imi == 1);
    blank_trial_inds = find(ismember(all_trialvec,blank_trials));
    sim_blanks = cell(length(sim_sac_times),1);
    for ii = 1:length(sim_sac_times)
        sim_blanks{ii} = blank_trial_inds(all_tsince_start(blank_trial_inds(1:end-1)) < sim_sac_times(ii) & ...
            all_tsince_start(blank_trial_inds(2:end)) >= sim_sac_times(ii));
        all_sim_blanks = [all_sim_blanks; sim_blanks{ii}];
    end
end

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
    Ngback_inds = sum(ismember(all_trialvec(used_inds),grayback_gs_trials));
    Niback_inds = sum(ismember(all_trialvec(used_inds),imback_trials));
else
    gback_sacs = find(ismember(all_blockvec(sac_start_inds),grayback_gs_expts));
    iback_sacs = find(ismember(all_blockvec(sac_start_inds),imback_gs_expts));
    ss_sacs = find(ismember(all_blockvec(sac_start_inds),sim_sac_expts));
    Ngback_inds = sum(ismember(all_blockvec(used_inds),grayback_gs_expts));
    Niback_inds = sum(ismember(all_blockvec(used_inds),imback_gs_expts));
end

gray_msac_set = intersect(gback_sacs,micro_set);
im_msac_set = intersect(iback_sacs,micro_set);

%% SACCADE TIMING ANALYSIS
sac_sm = round(0.025/dt);
binned_msacs = hist(sac_start_inds(micro_set),1:length(all_t_axis));
binned_gsacs = hist(sac_start_inds(gsac_set),1:length(all_t_axis));
binned_msac_sm = jmm_smooth_1d_cor(binned_msacs,sac_sm);
binned_gsac_sm = jmm_smooth_1d_cor(binned_gsacs,sac_sm);

maxlag = round(0.5/dt);
[gen_data.msac_acorr,acorr_lags] = xcov(binned_msac_sm,maxlag,'coeff');
[gen_data.gsac_acorr,acorr_lags] = xcov(binned_gsac_sm,maxlag,'coeff');
[gen_data.msac_gsac_xcorr,acorr_lags] = xcov(binned_gsac_sm,binned_msac_sm,maxlag,'coeff');
gen_data.acorr_lags = acorr_lags*dt;

raw_sac_start_inds = round(interp1(all_eye_ts,1:length(all_eye_ts),[saccades(gsac_set).start_time]));
[gen_data.gsac_rawtavg_eyespeed,eye_tavg_lags] = get_event_trig_avg_v3(all_eye_speed,raw_sac_start_inds,round(0.15*et_params.eye_fs),round(0.15*et_params.eye_fs));
raw_sac_start_inds = round(interp1(all_eye_ts,1:length(all_eye_ts),[saccades(micro_set).start_time]));
[gen_data.msac_rawtavg_eyespeed,eye_tavg_lags] = get_event_trig_avg_v3(all_eye_speed,raw_sac_start_inds,round(0.15*et_params.eye_fs),round(0.15*et_params.eye_fs));
gen_data.raw_eye_lags = eye_tavg_lags/et_params.eye_fs;

[gen_data.gsac_tavg_eyespeed,tavg_lags] = get_event_trig_avg_v3(interp_eye_speed,sac_start_inds(gsac_set),round(0.15/dt),round(0.15/dt));
[gen_data.msac_tavg_eyespeed,tavg_lags] = get_event_trig_avg_v3(interp_eye_speed,sac_start_inds(micro_set),round(0.15/dt),round(0.15/dt));
gen_data.eye_lags = tavg_lags*dt;

gen_data.N_blinks = sum(used_is_blink);
gen_data.N_msacs = length(micro_set);
gen_data.N_gsacs = length(gsac_set);
gen_data.N_simsacs = length(all_sim_sacs);
gen_data.N_simmsacs = length(all_sim_msacs);
gen_data.N_blanks = length(all_sim_blanks);
gen_data.N_msacs_gray = length(intersect(micro_set,gback_sacs));
gen_data.N_gsacs_gray = length(intersect(gsac_set,gback_sacs));
gen_data.N_msacs_im = length(intersect(micro_set,iback_sacs));
gen_data.N_gsacs_im = length(intersect(gsac_set,iback_sacs));
gen_data.N_msacs_Par = length(msac_Par);
gen_data.N_msacs_Orth = length(msac_Orth);

gen_data.Tot_time = length(used_inds);
gen_data.Gback_time = Ngback_inds;
gen_data.Iback_time = Niback_inds;

%%
trig_avg_params = struct('dt',dt,'min_trial_dur',min_trial_dur,'beg_buffer',beg_buffer,...
    'end_buffer',end_buffer,'good_coils',good_coils,'bar_ori',bar_ori,'nboot',nboot,'micro_thresh',micro_thresh);

cd(save_dir)
sname = 'gen_trig_avg_data';
if include_bursts
    sname = [sname '_withburst'];
end
sname = [sname sprintf('_ori%d',bar_ori)];
save(sname,'gen_data','trig_avg_params');
