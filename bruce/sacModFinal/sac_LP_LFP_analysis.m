clear all
close all

global Expt_name

Expt_name = 'M275';

%%
Expt_num = str2num(Expt_name(2:end));

if Expt_num > 280 && Expt_num < 289
    data_dir = ['/media/NTlab_data2/Data/bruce/' Expt_name];
elseif Expt_num >= 289
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end

save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
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
        case 297
            if ~ismember(bar_ori,[0 90])
                error('M297 is either 0 or 90 deg bars');
            end
    end
    good_coils = [1 1]; %which coils are usable
    
else
    error('Unrecognized experiment name');
end

% TBT_expts = {'M275','M277'}; %list of expts where conditions are interleaved
is_TBT_expt = false;

if strcmp(rec_type,'LP')
    if Expt_num >= 275
        rpt_seed = 1001; %M275 M277 M281
    else
        rpt_seed = 1e4; %m270 and 266
    end
else
    rpt_seed = nan;
end

ignore_blocks = [];
switch Expt_num
    case 270
        ignore_blocks = [5 19];
    case 275
        ignore_blocks = [10];
    case 289
        ignore_blocks = [26 27]; %27 is off somehow; 26 LFPs are messed up
    case 294
        ignore_blocks = [37 38 39]; %37-39 have slightly different dw used in these blocks
    case 86
        ignore_blocks = [16 17 28 30];
    case 87
        ignore_blocks = [15];
    case 93
        ignore_blocks = [28];
end

use_coils = [1 0];

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.005;

min_trial_dur = 1;
beg_buffer = 0.25;
end_buffer = 0.05;

backlag = 0.25;
forwardlag = 0.55;

sua_sm_sig = (0.005/dt);
mua_sm_sig = (0.005/dt);

if strcmp(Expt_name,'G081')
    trial_dur = 2;
else
    trial_dur = 4;
end

EP_bounds = 1;%eye position boundary
micro_thresh = 1; %microsaccade amplitude threshold (deg)
max_sac_dur = 0.1; %maximum saccade duration (otherwise likely a blink)
sac_burst_isi = 0.15; %minimum inter-saccade interval before classifying sac as part of a 'burst'

%% LOAD EXPTS STRUCT
cd(data_dir)
if Expt_name(1) == 'G'
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
else
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
end
if ~strcmp(Expt_name,'G081')
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
    
    if strcmp(Expt_name,'G093')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    end
    if isfield(Expts{cur_block}.Trials(1),'Fs')
        trial_Fs = [Expts{cur_block}.Trials(:).Fs];
    else
        trial_Fs = nan(1,length(trial_durs));
    end
    
    if is_TBT_expt
        if isfield(Expts{cur_block}.Trials,'Bs')
            trial_back = strcmp('image',{Expts{cur_block}.Trials(:).Bs});
        elseif isfield(Expts{cur_block}.Trials,'exvals')
            exvals = reshape([Expts{cur_block}.Trials(:).exvals],3,[]);
            trial_back = exvals(3,:) ~= 2; %these are the greyback trials
        else
            fprintf('Warning no Bs or exvals for this block, skipping\n');
            use_trials = [];
        end
        trial_back = trial_back(id_inds);
        
        trial_Ff = [Expts{cur_block}.Trials(:).Ff];
        trial_Ff = trial_Ff(id_inds);
        
        all_trial_back = cat(1,all_trial_back,trial_back(use_trials)');
        all_trial_Ff = cat(1,all_trial_Ff,trial_Ff(use_trials)');
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
    trial_cnt = trial_cnt + n_trials;
    if strcmp(rec_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

trial_start_inds = [1; find(diff(all_trialvec) > 0)+1];
simsac_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),sim_sac_expts));
grayback_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),grayback_gs_expts));
imback_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),imback_gs_expts));

%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = n_probes;
% if strcmp(rec_type,'LP')
%     clust_params.exclude_adjacent = true;
% else
clust_params.exclude_adjacent = false;
% end
[all_binned_mua,all_binned_sua,Clust_data] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

%% SMOOTH AND NORMALIZE SPIKING DATA
all_mua_rate = nan(size(all_binned_mua));
mua_block_mean_rates = nan(length(cur_block_set),n_probes);
mua_block_n_spikes = nan(length(cur_block_set),n_probes);
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
        for cc = 1:n_probes
            all_mua_rate(cur_block_inds,cc) = jmm_smooth_1d_cor(all_binned_mua(cur_block_inds,cc),mua_sm_sig);
        end
        mua_block_mean_rates(ee,:) = mean(all_mua_rate(cur_block_inds,:));
        mua_block_n_spikes(ee,:) = sum(all_binned_mua(cur_block_inds,:));
    end
end

all_sua_rate = nan(size(all_binned_sua));
sua_block_mean_rates = nan(length(cur_block_set),length(SU_numbers));
sua_block_n_spikes = nan(length(cur_block_set),length(SU_numbers));
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
        for ss = 1:length(SU_numbers)
            all_sua_rate(cur_block_inds,ss) = jmm_smooth_1d_cor(all_binned_sua(cur_block_inds,ss),sua_sm_sig);
        end
        sua_block_mean_rates(ee,:) = mean(all_sua_rate(cur_block_inds,:));
        sua_block_n_spikes(ee,:) = sum(all_binned_sua(cur_block_inds,:));
    end
end

%normalized firing rates (smoothed)
all_sua_rate_norm = nan(size(all_sua_rate));
all_mua_rate_norm = nan(size(all_mua_rate));
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
        all_sua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_sua_rate(cur_block_inds,:),...
            sua_block_mean_rates(ee,:));
        all_mua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_mua_rate(cur_block_inds,:),...
            mua_block_mean_rates(ee,:));
    end
end

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

%% PARSE EYETRACKING DATA
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
micro_set(ismember(micro_set,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'

%saccade amplitude along parallel axis
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

%guided saccades are those whose parallel component is large enough and
%that aren't blinks
gsac_set = find(abs(sac_deltaX) > gsac_thresh & ~used_is_blink' & ~out_bounds);

%classify guided saccades as 'in' vs 'out' and 'pos' vs 'neg'
outsacs = gsac_set(abs(sac_prepos(1,gsac_set)) < abs(sac_postpos(1,gsac_set)));
insacs = gsac_set(abs(sac_prepos(1,gsac_set)) > abs(sac_postpos(1,gsac_set)));
negsacs = gsac_set(sac_postpos(1,gsac_set) < sac_prepos(1,gsac_set));
possacs = gsac_set(sac_postpos(1,gsac_set) > sac_prepos(1,gsac_set));
out_pos_sacs = intersect(outsacs,possacs);
out_neg_sacs = intersect(outsacs,negsacs);
in_pos_sacs = intersect(insacs,possacs);
in_neg_sacs = intersect(insacs,negsacs);

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
else
    gback_sacs = find(ismember(all_blockvec(sac_start_inds),grayback_gs_expts));
    iback_sacs = find(ismember(all_blockvec(sac_start_inds),imback_gs_expts));
    ss_sacs = find(ismember(all_blockvec(sac_start_inds),sim_sac_expts));
end

gray_msac_set = intersect(gback_sacs,micro_set);
im_msac_set = intersect(iback_sacs,micro_set);


%%
cd(data_dir);
if strcmp(rec_type,'LP')
    Fs = 1000;
    dsf = 4;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    [bb,aa] = butter(4,[0.5 niqf/dsf*0.8]/niqf);
else
    Fs = 400;
    dsf = 2;
    use_lfps = 1:2:n_probes;
end

full_lfps = [];
full_lfp_taxis = [];
cur_toffset = 0;
for ee = 1:length(cur_block_set);
    
    fprintf('Loading LFPs, Expt %d of %d\n',ee,length(cur_block_set));
    
    if strcmp(rec_type,'LP')
        fname = sprintf('lemM%dA.%d.lfp.mat',Expt_num,cur_block_set(ee));
        load(fname);
        cur_Fs = 1/LFP.Header.CRsamplerate;
        if Fs ~= cur_Fs
            error('Fs mismatch');
        end
        
        tlens = arrayfun(@(X) length(X.ftime),LFP.Trials);
        bad_trials = find(tlens == 0);
        LFP.Trials(bad_trials) = [];
        n_trials(ee) = length(LFP.Trials);
        lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
        lfp_trial_ends = [LFP.Trials(:).End]/1e4;
        expt_lfp_t_axis = [];
        expt_lfps = [];
        for tt = 1:n_trials(ee)
            %         tt
            cur_npts = size(LFP.Trials(tt).LFP,1);
            cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
            cur_t_axis = (lfp_trial_starts(tt):1/Fs:cur_t_end(tt)) + cur_toffset;
            
            if ~isempty(expt_lfp_t_axis)
                cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
            else
                cur_sp = 1;
            end
            cur_t_axis = cur_t_axis(cur_sp:end);
            if length(cur_t_axis) > 50
                cur_LFP = double([LFP.Trials(tt).LFP]);
                cur_LFP = cur_LFP(cur_sp:end,:);
                cur_LFP = filtfilt(bb,aa,cur_LFP);
                
                cur_LFP = downsample(cur_LFP,dsf);
                cur_t_axis = downsample(cur_t_axis,dsf);
                
                if size(cur_LFP,2) == n_probes
                    expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
                    expt_lfps = [expt_lfps; cur_LFP];
                end
            end
        end
    else
        lfp_fname = sprintf('Expt%d_LFP.mat',cur_block_set(ee));
        load(lfp_fname);
        Fs = lfp_params.Fsd;
        niqf = Fs/2;
        %         [bb,aa] = butter(4,[1/dsf*0.8],'low');
        [bb,aa] = butter(4,[0.5 niqf/dsf*0.8]/niqf);
        
        cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
        for ii = 1:length(use_lfps)
            cur_lfps(:,ii) = filtfilt(bb,aa,cur_lfps(:,ii));
        end
        expt_lfps = downsample(cur_lfps,dsf);
        expt_lfp_t_axis = downsample(lfp_t_ax',dsf);
    end
    
    cur_uset = find(all_blockvec == ee);
    if ~isempty(cur_uset)
        uinds = find(expt_lfp_t_axis >= all_t_axis(cur_uset(1)) & expt_lfp_t_axis <= all_t_axis(cur_uset(end)));
        full_lfps = cat(1,full_lfps,expt_lfps(uinds,:));
        full_lfp_taxis = cat(1,full_lfp_taxis,expt_lfp_t_axis(uinds));
    end
    cur_toffset = trial_toffset(ee);
end

%%
lfp_trial_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_trial_start_times));
lfp_trial_stop_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_trial_end_times));
lfp_sac_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),sac_start_times));
lfp_simsac_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_t_axis(all_sim_sacs)));
lfp_trial_inds = [lfp_trial_start_inds(:) lfp_trial_stop_inds(:)];
% lfp_trial_start_inds(isnan(lfp_trial_start_inds)) = [];
% lfp_trial_stop_inds(isnan(lfp_trial_stop_inds)) = [];

interp_lfps = interp1(full_lfp_taxis,full_lfps,all_t_axis);
%% MAKE LOOPED COHERANCE SEGC CALCULATOR
% addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))
% params.Fs = 1/dt;
% params.tapers = [3 5];
% [C,phi,S12,f] = cohmatrixc(interp_lfps(used_inds,:),params);

%% COMPUTE TRIG AVGS FOR LFPs
lforwardlag = round(forwardlag*Fsd);
lbacklag = round(backlag*Fsd);

nboot = [];
%set trial numbers to Inf so they don't get included in trig averaging
used_trialvec = ones(size(all_trialvec))*Inf;
used_trialvec(used_inds) = all_trialvec(used_inds);
used_trialvec_interp = interp1(all_t_axis,used_trialvec,full_lfp_taxis);
trialvec_interp = round(interp1(all_t_axis,all_trialvec,full_lfp_taxis));
tsince_start_interp = interp1(all_t_axis,all_tsince_start,full_lfp_taxis);

csd_params.Fs = Fsd; %sample freq
csd_params.BrainBound = 1; %first channel that is in the brain
csd_params.ChanSep = 0.05; %channel sep in mm
csd_params.diam = 2; %current disc diameter (in mm)

[lfp_data.trial_onset_csd,lags] = get_event_trig_csd(full_lfps,lfp_trial_start_inds,lbacklag,lforwardlag,csd_params);
[lfp_data.trial_offset_csd,lags] = get_event_trig_csd(full_lfps,lfp_trial_stop_inds,lbacklag,lforwardlag,csd_params);
[lfp_data.trial_onset_lfp,lags] = get_event_trig_avg(full_lfps,lfp_trial_start_inds,lbacklag,lforwardlag,csd_params);
[lfp_data.trial_offset_lfp,lags] = get_event_trig_avg(full_lfps,lfp_trial_stop_inds,lbacklag,lforwardlag,csd_params);

%general averages
[lfp_data.msac_avg,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(micro_set),lbacklag,lforwardlag,nboot,used_trialvec_interp,0);
[lfp_data.gsac_avg,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(gsac_set),lbacklag,lforwardlag,nboot,used_trialvec_interp,0);
[lfp_data.simsac_avg,lags] = get_event_trig_avg(full_lfps,lfp_simsac_start_inds,lbacklag,lforwardlag,nboot,used_trialvec_interp,0);

if strcmp(rec_type,'LP')
    [lfp_data.gsac_csd,lags] = get_event_trig_csd(full_lfps,lfp_sac_start_inds(gsac_set),lbacklag,lforwardlag,csd_params,used_trialvec_interp,0);
    [lfp_data.msac_csd,lags] = get_event_trig_csd(full_lfps,lfp_sac_start_inds(micro_set),lbacklag,lforwardlag,csd_params,used_trialvec_interp,0);
    [lfp_data.simsac_csd,lags] = get_event_trig_csd(full_lfps,lfp_simsac_start_inds,lbacklag,lforwardlag,csd_params,used_trialvec_interp,0);
end

%%
lforwardlag = round(forwardlag/dt);
lbacklag = round(backlag/dt);

trial_start_inds = [1; 1+find(diff(all_trialvec) > 0)];
[lfp_data.trial_onset_mua,mlags] = get_event_trig_avg(all_mua_rate_norm,trial_start_inds,lbacklag,lforwardlag); 

%%
%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 25;
min_freq = 2; max_freq = 80;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

ampgrams = nan(size(full_lfps,1),nwfreqs,n_probes);
% wave_trans = nan(size(full_lfps,1),nwfreqs,n_probes);
for ll = 1:n_probes
    ll
    temp = cwt(full_lfps(:,ll),scales,'cmor1-1');
    ampgrams(:,:,ll) = abs(temp)';
%     wave_trans(:,:,ll) = cwt(full_lfps(:,ll),scales,'cmor1-1')';
end

% wave_coh = nan(n_probes,n_probes,length(wfreqs));
% for ii  = 1:n_probes-1
%     ii
%     for jj = (ii+1):n_probes
%         cross_wave = squeeze(abs(mean(wave_trans(:,:,ii).*conj(wave_trans(:,:,jj)))));
%         wave_coh(ii,jj,:) = cross_wave./sqrt(mean(abs(wave_trans(:,:,ii).^2.*wave_trans(:,:,jj).^2)));
%     
%     end
% end


ampgrams = sqrt(ampgrams);
lfp_data.ampgrams_std = std(ampgrams);
lfp_data.ampgrams_mean = mean(ampgrams);
% ampgrams = zscore(ampgrams);

[lfp_data.onset_specgram,lags] = get_event_trig_avg(ampgrams,lfp_trial_start_inds,lbacklag,lforwardlag);
[lfp_data.sac_specgram,lags] = get_event_trig_avg(ampgrams,lfp_sac_start_inds(gsac_set),lbacklag,lforwardlag);
[lfp_data.msac_specgram,lags] = get_event_trig_avg(ampgrams,lfp_sac_start_inds(gsac_set),lbacklag,lforwardlag);

% onset_specgramZ = bsxfun(@rdivide,bsxfun(@minus,onset_specgram,ampgrams_mean),ampgrams_std);
% sac_specgramZ = bsxfun(@rdivide,bsxfun(@minus,sac_specgram,ampgrams_mean),ampgrams_std);
%%
cd(save_dir)
sname = 'lfp_trig_avgs';
save(sname,'lags','Fsd','lfp_data','mlags','dt','wfreqs');