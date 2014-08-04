clear all
close all

addpath('~/James_scripts/CircStat2011f/')
global Expt_name bar_ori

Expt_name = 'G086';
bar_ori = 90;

%%
Expt_num = str2num(Expt_name(2:end));

if Expt_num >= 280
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
        case 232
            bar_ori = 50;
        case 235
            bar_ori = 30;
        case 239
            bar_ori = 130;
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

min_trial_dur = 2;
beg_buffer = 0.25;
end_buffer = 0.05;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.005/dt);
mua_sm_sig = (0.005/dt);
% mua_sm_sig = 0;

if strcmp(Expt_name,'G081') || ismember(Expt_num,[232 235 239])
    trial_dur = 2;
else
    trial_dur = 4;
end

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
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa','rls.FaXwi','rls.FaXwiXimi'};
if any(strcmp(Expt_name,{'G095','M266','M270'}))
    is_sim_msac_expt = true;
    % elseif any(strcmp(Expt_name,TBT_expts))
    %     is_TBT_expt = true;
end
if Expt_num >= 275
    is_TBT_expt = true;
end

if ismember(Expt_num,[232 235 239])
    load ./random_bar_eyedata_ftime.mat bar_expts
    load ./bar_params.mat
    cur_block_set = bar_expts;
    expt_bar_ori = ones(1,length(Expts))*bar_ori;
    expt_has_ds = ones(1,length(Expts))';
    expt_imback = zeros(1,length(Expts));
    expt_sac_dir = expt_bar_ori;
    
    has_data = cellfun(@(x) length(x),Expts) > 0;
    expt_names(has_data) = cellfun(@(x) x.Header.expname,Expts(has_data),'UniformOutput',false);
    expt_sac_dir(has_data) = mod(cellfun(@(x) x.Stimvals.Fa,Expts(has_data),'UniformOutput',true),180);
    expt_sac_amp(has_data) = cellfun(@(x) x.Stimvals.Fs,Expts(has_data),'UniformOutput',true);
    fprintf('Sac amp: %d Sac dir: %d\n',expt_sac_amp(cur_block_set(1)),expt_sac_dir(cur_block_set(1)));
else
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
    
    if strcmp(rec_type,'UA')
        if strcmp(Expt_name,'G081')
            %         cur_block_set = find(included_type & expt_Fr == 1 & (expt_bar_ori == 0 | expt_bar_ori == 90 | expt_bar_ori == 45 | expt_bar_ori == 135));
%             cur_block_set = find(included_type & expt_Fr == 1 & expt_bar_ori == bar_ori);
            cur_block_set = find(included_type & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_ce == 1);
        else
            %         cur_block_set = find(included_type & expt_Fr == 1 & (expt_bar_ori == 0 | expt_bar_ori == 90));
            cur_block_set = find(included_type & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_ce == 1);
        end
    else
%         cur_block_set = find(included_type & expt_ce == 1 & expt_Fr == 1);
        cur_block_set = find(included_type & expt_ce == 1 & expt_Fr == 1 & expt_ce == 1);
    end
    if strcmp(Expt_name,'G081')
        expt_has_ds = (expt_ijump==0)';
    end
    expt_has_ds(isnan(expt_has_ds)) = 0;
end

cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

sim_sac_expts = find(expt_has_ds(cur_block_set) ~= 1);
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');
% or_expts = find(strcmp('rls.orRC',expt_names(cur_block_set)));
% expt_bar_ori(cur_block_set(or_expts)) = 0;

poss_orth_expts = find(mod(expt_bar_ori(cur_block_set) - expt_sac_dir(cur_block_set),180) == 90);
if ~isempty(poss_orth_expts)
    fprintf('Warning, possible orthoganol saccade expts detected\n');
end

gsac_amp = unique(expt_sac_amp(cur_block_set([grayback_gs_expts; imback_gs_expts])));
if length(gsac_amp) > 1
    fprintf('Multiple guided sac amps detected!\n');
end
gsac_thresh = mean(gsac_amp)/3;

if strcmp(Expt_name,'G081')
    sim_sac_times = [0.7 1.4];
else
    sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
end

if is_TBT_expt
    sim_sac_expts = []; imback_gs_expts = []; grayback_gs_expts = [];
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
    
    if strcmp(Expt_name,'G093')
        if isfield(Expts{cur_block}.Trials,'wi')
            trial_wi = [Expts{cur_block}.Trials(:).wi];
            trial_wi = trial_wi(id_inds);
        else
            trial_wi = ones(1,length(id_inds))*1.9959;
        end
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
        else
            trial_back = nan(1,length(use_trials));
            %            use_trials = [];
        end
        trial_back = trial_back(id_inds);
        
        if isfield(Expts{cur_block}.Trials,'exvals')
            exvals = reshape([Expts{cur_block}.Trials(:).exvals],3,[]);
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
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;
SU_isodist = Clust_data.SU_isodists;
SU_Lratio = Clust_data.SU_Lratios;

%% SMOOTH AND NORMALIZE SPIKING DATA
all_mua_rate = nan(size(all_binned_mua));
mua_block_mean_rates = nan(length(cur_block_set),n_probes);
mua_block_n_spikes = nan(length(cur_block_set),n_probes);
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
        for cc = 1:n_probes
            if mua_sm_sig > 0
                all_mua_rate(cur_block_inds,cc) = jmm_smooth_1d_cor(all_binned_mua(cur_block_inds,cc),mua_sm_sig);
            else
                all_mua_rate(cur_block_inds,cc) = all_binned_mua(cur_block_inds,cc);
            end
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

%% PROCESS EYE TRACKING DATA
% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,good_coils);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

% %compute corrected eye data in bar-oriented frame
% [corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
%     all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);
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

% used_sac_set = find(ismember(sac_start_inds,used_inds) & ismember(sac_stop_inds,used_inds));

sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));

%find saccades that start or end out of window
bounds = 1;
out_bounds = abs(sac_prepos(2,:)) > bounds | abs(sac_postpos(2,:)) > bounds;

micro_thresh = 1;
max_sac_dur = 0.1;
sac_burst_isi = 0.15;
sacburst_set = find([saccades(:).isi] < sac_burst_isi | [saccades(:).next_isi] < sac_burst_isi);
micro_set = find([saccades(:).amplitude] < micro_thresh & ~used_is_blink' & ~out_bounds);
msac_bursts = micro_set(ismember(micro_set,sacburst_set));
micro_set(ismember(micro_set,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'

sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);
gsac_set = find(abs(sac_deltaX) > gsac_thresh & ~used_is_blink' & ~out_bounds);

outsacs = gsac_set(abs(sac_prepos(1,gsac_set)) < abs(sac_postpos(1,gsac_set)));
insacs = gsac_set(abs(sac_prepos(1,gsac_set)) > abs(sac_postpos(1,gsac_set)));
% negsacs = outsacs(sac_postpos(1,outsacs) < 0);
% possacs = outsacs(sac_postpos(1,outsacs) > 0);
negsacs = gsac_set(sac_postpos(1,gsac_set) < sac_prepos(1,gsac_set));
possacs = gsac_set(sac_postpos(1,gsac_set) > sac_prepos(1,gsac_set));
out_pos_sacs = intersect(outsacs,possacs);
out_neg_sacs = intersect(outsacs,negsacs);
in_pos_sacs = intersect(insacs,possacs);
in_neg_sacs = intersect(insacs,negsacs);

%compile indices of simulated saccades
all_sim_sacs = [];
if is_TBT_expt
    sim_sac_trials = find(all_trial_Ff == 0);
    sim_trial_inds = find(ismember(all_trialvec,sim_sac_trials));
    sim_sacs = cell(length(sim_sac_times),1);
    for ii = 1:length(sim_sac_times)
        sim_sacs{ii} = sim_trial_inds(all_tsince_start(sim_trial_inds(1:end-1)) < sim_sac_times(ii) & ...
            all_tsince_start(sim_trial_inds(2:end)) >= sim_sac_times(ii));
        all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
    end
else
    sim_expt_inds = find(ismember(all_blockvec,sim_sac_expts));
    sim_sacs = cell(length(sim_sac_times),1);
    for ii = 1:length(sim_sac_times)
        sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
            all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
        all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
    end
end
all_sim_sacs = sort(all_sim_sacs(ismember(all_sim_sacs,used_inds)));

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
% or_sacs = find(ismember(all_blockvec(sac_start_inds),or_expts));

sac_oris = expt_bar_ori(cur_block_set(all_blockvec(sac_start_inds)));
sim_sac_oris = expt_bar_ori(cur_block_set(all_blockvec(all_sim_sacs)));
% sac_dds = expt_dd(cur_block_set(all_blockvec(sac_start_inds)));
% sac_ce = expt_ce(cur_block_set(all_blockvec(sac_start_inds)));

gray_msac_set = intersect(gback_sacs,micro_set);
% or_msac_set = intersect(micro_set,or_sacs);
% msac_set = setdiff(micro_set,or_sacs);

msac_dirs = [saccades(micro_set).direction];
large_msacs = micro_set([saccades(micro_set).amplitude] > 0.5);
small_msacs = micro_set([saccades(micro_set).amplitude] < 0.5);
rf_angle = atan2(Expts{cur_block_set(1)}.Stimvals.rf(2),Expts{cur_block_set(1)}.Stimvals.rf(1));
msac_dirs_relrf = abs(circ_dist(msac_dirs,rf_angle));
msac_towards = micro_set(msac_dirs_relrf <= pi/4);
msac_aways = micro_set(msac_dirs_relrf >= 3*pi/4);

n_msac_bins = 8;
msac_bin_edges = linspace(-pi,pi,n_msac_bins+1);
msac_bin_cents = 0.5*msac_bin_edges(1:end-1) + 0.5*msac_bin_edges(2:end);

%% FOR TBT experiments, normalize firing rates separately within each trial-type
if is_TBT_expt
    grayback_trial_inds = find(ismember(all_trialvec,grayback_gs_trials));
    all_mua_rate_norm(grayback_trial_inds,:) = bsxfun(@rdivide,all_mua_rate_norm(grayback_trial_inds,:),nanmean(all_mua_rate_norm(grayback_trial_inds,:)));
    all_sua_rate_norm(grayback_trial_inds,:) = bsxfun(@rdivide,all_sua_rate_norm(grayback_trial_inds,:),nanmean(all_sua_rate_norm(grayback_trial_inds,:)));
    
    imback_trial_inds = find(ismember(all_trialvec,imback_gs_trials));
    all_mua_rate_norm(imback_trial_inds,:) = bsxfun(@rdivide,all_mua_rate_norm(imback_trial_inds,:),nanmean(all_mua_rate_norm(imback_trial_inds,:)));
    all_sua_rate_norm(imback_trial_inds,:) = bsxfun(@rdivide,all_sua_rate_norm(imback_trial_inds,:),nanmean(all_sua_rate_norm(imback_trial_inds,:)));

    simsac_trial_inds = find(ismember(all_trialvec,sim_sac_trials));
    all_mua_rate_norm(simsac_trial_inds,:) = bsxfun(@rdivide,all_mua_rate_norm(simsac_trial_inds,:),nanmean(all_mua_rate_norm(simsac_trial_inds,:)));
    all_sua_rate_norm(simsac_trial_inds,:) = bsxfun(@rdivide,all_sua_rate_norm(simsac_trial_inds,:),nanmean(all_sua_rate_norm(simsac_trial_inds,:)));
end

%% SACCADE TIMING ANALYSIS
sac_sm = round(0.025/dt);
binned_msacs = hist(sac_start_inds(micro_set),1:length(all_t_axis));
binned_gsacs = hist(sac_start_inds(gsac_set),1:length(all_t_axis));
binned_msac_sm = jmm_smooth_1d_cor(binned_msacs,sac_sm);
binned_gsac_sm = jmm_smooth_1d_cor(binned_gsacs,sac_sm);

%trial averages (USE ALL OF TRIAL HERE)
[gen_data.msac_simsac_trial_avg,trial_lags] = get_event_trig_avg_v3(binned_msac_sm,simsac_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.msac_grayback_trial_avg,trial_lags] = get_event_trig_avg_v3(binned_msac_sm,grayback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.msac_imback_trial_avg,trial_lags] = get_event_trig_avg_v3(binned_msac_sm,imback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.gsac_simsac_trial_avg,trial_lags] = get_event_trig_avg_v3(binned_gsac_sm,simsac_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.gsac_grayback_trial_avg,trial_lags] = get_event_trig_avg_v3(binned_gsac_sm,grayback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.gsac_imback_trial_avg,trial_lags] = get_event_trig_avg_v3(binned_gsac_sm,imback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);

maxlag = round(0.5/dt);
[gen_data.msac_acorr,acorr_lags] = xcov(binned_msac_sm,maxlag,'coeff');
[gen_data.gsac_acorr,acorr_lags] = xcov(binned_gsac_sm,maxlag,'coeff');
[gen_data.msac_gsac_xcorr,acorr_lags] = xcov(binned_gsac_sm,binned_msac_sm,maxlag,'coeff');

gen_data.N_msacs = length(micro_set);
gen_data.N_gsacs = length(gsac_set);
gen_data.N_simsacs = length(all_sim_sacs);
gen_data.N_msacs_gray = length(intersect(micro_set,gback_sacs));
gen_data.N_gsacs_gray = length(intersect(gsac_set,gback_sacs));
gen_data.N_msacs_im = length(intersect(micro_set,iback_sacs));
gen_data.N_gsacs_im = length(intersect(gsac_set,iback_sacs));

%%
% cbacklag = backlag;
% cforlag = forwardlag;
% slags = -cbacklag:cforlag;
% n_sac_bins = length(slags);
% 
% saccade_trial_inds = all_trialvec(sac_start_inds);
% 
% NT = length(all_trialvec);
% Xsac = zeros(NT,length(slags));
% Xmsac = zeros(NT,length(slags));
% for ii = 1:n_sac_bins
%     cur_sac_target = sac_start_inds(gsac_set) + slags(ii);
%     uu = find(cur_sac_target > 1 & cur_sac_target < NT);
%     cur_sac_target = cur_sac_target(uu);
%     cur_sac_target(all_trialvec((cur_sac_target)) ~= saccade_trial_inds(gsac_set(uu))) = [];
%     Xsac(cur_sac_target,ii) = 1;
%     
%     cur_sac_target = sac_start_inds(micro_set) + slags(ii);
%     uu = find(cur_sac_target > 1 & cur_sac_target < NT);
%     cur_sac_target = cur_sac_target(uu);
%     cur_sac_target(all_trialvec((cur_sac_target)) ~= saccade_trial_inds(micro_set(uu))) = [];
%     Xmsac(cur_sac_target,ii) = 1;
% end
% any_sac_inds = find(any(Xsac > 0,2) | any(Xmsac > 0,2));
% 
% X{1} = Xsac(any_sac_inds,:);
% X{2} = Xmsac(any_sac_inds,:);
% clear Xsac Xmsac
% 
%%

% stim_params(1) = NMMcreate_stim_params(length(slags));
% stim_params(2) = NMMcreate_stim_params(length(slags));
% reg_params = NMMcreate_reg_params('lambda_d2T',100,'lambda_L2',5);
% gsac_filts = nan(n_probes,length(slags));
% msac_filts = nan(n_probes,length(slags));
% for cc = 1:n_probes;
%     fprintf('Fitting filter model for MUA %d of %d\n',cc,n_probes);
%     cur_mod = NMMinitialize_model(stim_params,[1 1],{'lin','lin'},reg_params,[1 2]);
%     cur_mod = NMMfit_filters(cur_mod,all_binned_mua(any_sac_inds,cc),X);
%     gsac_filts(cc,:) = cur_mod.mods(1).filtK;
%     msac_filts(cc,:) = cur_mod.mods(2).filtK;
% end
% mua_data.gsac_filts = gsac_filts;
% mua_data.msac_filts = msac_filts;

%%

% gsac_filts = nan(length(SU_probes),length(slags));
% msac_filts = nan(length(SU_probes),length(slags));
% for cc = 1:length(SU_probes)
%     fprintf('Fitting filter model for SUA %d of %d\n',cc,length(SU_probes));
%     cur_Robs = all_binned_sua(any_sac_inds,cc);
%     uset = find(~isnan(cur_Robs));
%     if ~isempty(uset)
%         cur_mod = NMMinitialize_model(stim_params,[1 1],{'lin','lin'},reg_params,[1 2]);
%         cur_mod = NMMfit_filters(cur_mod,cur_Robs(uset),get_Xcell_tInds(X,uset));
%         gsac_filts(cc,:) = cur_mod.mods(1).filtK;
%         msac_filts(cc,:) = cur_mod.mods(2).filtK;
%     end
% end
% sua_data.gsac_filts = gsac_filts;
% sua_data.msac_filts = msac_filts;

%% COMPUTE TRIG AVGS FOR MUA
nboot = [];
%set trial numbers to Inf so they don't get included in trig averaging
used_trialvec = ones(size(all_trialvec))*Inf;
used_trialvec(used_inds) = all_trialvec(used_inds);

fprintf('Computing trig avgs for MUA\n',cc,n_probes);
%trial averages (USE ALL OF TRIAL HERE)
[mua_data.simsac_trial_avg,trial_lags] = get_event_trig_avg_v3(all_mua_rate_norm,simsac_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0);
[mua_data.grayback_trial_avg,trial_lags] = get_event_trig_avg_v3(all_mua_rate_norm,grayback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0);
[mua_data.imback_trial_avg,trial_lags] = get_event_trig_avg_v3(all_mua_rate_norm,imback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0);

%general averages
[mua_data.msac_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.simsac_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,all_sim_sacs,backlag,forwardlag,nboot,used_trialvec,0);

[mua_data.msac_burst_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(msac_bursts),backlag,forwardlag,nboot,used_trialvec,0);

[mua_data.large_msac_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(large_msacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.small_msac_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(small_msacs),backlag,forwardlag,nboot,used_trialvec,0);
% for ii = 1:n_msac_bins
%     cur_msac_set = micro_set(msac_dirs >= msac_bin_edges(ii) & msac_dirs < msac_bin_edges(ii+1));
%     [mua_data.msac_dir_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(cur_msac_set),backlag,forwardlag,nboot,used_trialvec,0);
% end

[mua_data.msac_end_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_stop_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_end_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_stop_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.msac_peak_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_peak_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_peak_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_peak_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);

%background dependent
[mua_data.msac_gray_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(intersect(micro_set,gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_gray_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(intersect(gsac_set,gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.msac_im_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(intersect(micro_set,iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_im_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(intersect(gsac_set,iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);

%sac-location dependent
[mua_data.gsac_out_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(outsacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_in_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(insacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_neg_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(negsacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_pos_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(possacs),backlag,forwardlag,nboot,used_trialvec,0);

[mua_data.gsac_outpos_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(out_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_inpos_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(in_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_outneg_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(out_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_inneg_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(in_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);

[mua_data.msac_towards_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(msac_towards),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.msac_away_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(msac_aways),backlag,forwardlag,nboot,used_trialvec,0);

mua_data.avg_rates = mean(all_binned_mua(used_inds,:));
mua_data.tot_nspikes = sum(all_binned_mua(used_inds,:));

%% COMPUTE TRIG AVGS FOR SUA

nboot = 200;
%general averages
[sua_data.msac_avg,lags,sua_data.msac_std] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_avg,lags,sua_data.gsac_std] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.simsac_avg,lags,sua_data.simsac_std] = get_event_trig_avg_v3(all_sua_rate_norm,all_sim_sacs,backlag,forwardlag,nboot,used_trialvec,0);

%background dependent
[sua_data.msac_gray_avg,lags,sua_data.msac_gray_std] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(intersect(micro_set,gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_gray_avg,lags,sua_data.gsac_gray_std] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(intersect(gsac_set,gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.msac_im_avg,lags,sua_data.msac_im_std] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(intersect(micro_set,iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_im_avg,lags,sua_data.gsac_im_std] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(intersect(gsac_set,iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);

%dont compute error bars for these
nboot = [];
[sua_data.msac_burst_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(msac_bursts),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.msac_end_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_stop_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_end_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_stop_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.msac_peak_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_peak_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_peak_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_peak_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);

%sac-location dependent
[sua_data.gsac_out_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(outsacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_in_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(insacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_neg_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(negsacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_pos_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(possacs),backlag,forwardlag,nboot,used_trialvec,0);

[sua_data.gsac_outpos_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(out_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_inpos_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(in_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_outneg_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(out_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_inneg_avg,lags] = get_event_trig_avg_v3(all_sua_rate_norm,sac_start_inds(in_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);

%%
block_SU_rates = nan(length(cur_block_set),length(SU_probes));
for ii = 1:length(cur_block_set)
    block_SU_rates(ii,:) = nanmean(all_binned_sua(all_blockvec == ii,:));
end

for ss = 1:length(SU_numbers)
    unit_data(ss).SU_numbers = SU_numbers(ss);
    unit_data(ss).probe_numbers = SU_probes(ss);
    
    unit_data(ss).avg_rates = nanmean(all_binned_sua(used_inds,ss));
    unit_data(ss).tot_nspikes = nansum(all_binned_sua(used_inds,ss));
    unit_data(ss).N_used_samps = nansum(~isnan(all_binned_sua(used_inds,ss)));
    
    unit_data(ss).stability_cv = nanstd(block_SU_rates(:,ss))./nanmean(block_SU_rates(:,ss));
    
    unit_data(ss).SU_Lratio = Clust_data.SU_Lratios(ss);
    unit_data(ss).SU_isodist = Clust_data.SU_isodists(ss);
    unit_data(ss).SU_refract = Clust_data.SU_refract(ss);
    unit_data(ss).SU_dprime = Clust_data.SU_dprimes(ss);
    
    unit_data(ss).N_msacs = sum(~isnan(all_binned_sua(sac_start_inds(micro_set),ss)));
    unit_data(ss).N_gsacs = sum(~isnan(all_binned_sua(sac_start_inds(gsac_set),ss)));
    unit_data(ss).N_simsacs = sum(~isnan(all_binned_sua(all_sim_sacs,ss)));
    unit_data(ss).N_msacs_gray = sum(~isnan(all_binned_sua(sac_start_inds(intersect(micro_set,gback_sacs)),ss)));
    unit_data(ss).N_gsacs_gray = sum(~isnan(all_binned_sua(sac_start_inds(intersect(gsac_set,gback_sacs)),ss)));
    unit_data(ss).N_msacs_im = sum(~isnan(all_binned_sua(sac_start_inds(intersect(micro_set,iback_sacs)),ss)));
    unit_data(ss).N_gsacs_im = sum(~isnan(all_binned_sua(sac_start_inds(intersect(gsac_set,iback_sacs)),ss)));
end
%%

sua_msacs = mat2cell(sua_data.msac_avg,length(lags),ones(1,length(SU_numbers)));
sua_gsacs = mat2cell(sua_data.gsac_avg,length(lags),ones(1,length(SU_numbers)));
sua_simsacs = mat2cell(sua_data.simsac_avg,length(lags),ones(1,length(SU_numbers)));
sua_msacs_SD = mat2cell(sua_data.msac_std,length(lags),ones(1,length(SU_numbers)));
sua_gsacs_SD = mat2cell(sua_data.gsac_std,length(lags),ones(1,length(SU_numbers)));
sua_simsacs_SD = mat2cell(sua_data.simsac_std,length(lags),ones(1,length(SU_numbers)));

sua_msacs_gray = mat2cell(sua_data.msac_gray_avg,length(lags),ones(1,length(SU_numbers)));
sua_gsacs_gray = mat2cell(sua_data.gsac_gray_avg,length(lags),ones(1,length(SU_numbers)));
sua_msacs_gray_SD = mat2cell(sua_data.msac_gray_std,length(lags),ones(1,length(SU_numbers)));
sua_gsacs_gray_SD = mat2cell(sua_data.gsac_gray_std,length(lags),ones(1,length(SU_numbers)));

sua_msacs_im = mat2cell(sua_data.msac_im_avg,length(lags),ones(1,length(SU_numbers)));
sua_gsacs_im = mat2cell(sua_data.gsac_im_avg,length(lags),ones(1,length(SU_numbers)));
sua_msacs_im_SD = mat2cell(sua_data.msac_im_std,length(lags),ones(1,length(SU_numbers)));
sua_gsacs_im_SD = mat2cell(sua_data.gsac_im_std,length(lags),ones(1,length(SU_numbers)));

sua_trig_avgs = struct('msac',sua_msacs,'gsac',sua_gsacs,'simsac',sua_simsacs,'msac_gray',sua_msacs_gray,...
'gsac_gray',sua_gsacs_gray,'msac_im',sua_msacs_im,'gsac_im',sua_gsacs_im,...
'msac_SD',sua_msacs_SD,'gsac_SD',sua_gsacs_SD,'simsac_SD',sua_simsacs_SD,'msac_gray_SD',sua_msacs_gray_SD,...
'gsac_gray_SD',sua_gsacs_gray_SD,'msac_im_SD',sua_msacs_im_SD,'gsac_im_SD',sua_gsacs_im_SD,'unit_data',num2cell(unit_data));


%%
trig_avg_params = struct('mua_sm_sig',mua_sm_sig,'sua_sm_sig',sua_sm_sig,'dt',dt,'lags',lags,'min_trial_dur',min_trial_dur,...
    'beg_buffer',beg_buffer,'end_buffer',end_buffer,'good_coils',good_coils,'bar_ori',bar_ori);

cd(save_dir)
sname = 'sac_trig_avg_data6';
if strcmp(rec_type,'UA') && bar_ori == 90
    sname = [sname '_vbars'];
end
save(sname,'sua_data','mua_data','gen_data','trig_avg_params','sua_trig_avgs');

