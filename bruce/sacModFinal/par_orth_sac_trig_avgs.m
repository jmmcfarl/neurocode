clear all
close all

global Expt_name bar_ori

Expt_name = 'G089';

fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';
save_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod'];

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

min_trial_dur = 0.75;
beg_buffer = 0.2;
end_buffer = 0.05;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.005/dt);

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
if ~strcmp(Expt_name,'G081')
    load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
end

%% PARSE EXPTS STRUCT
is_sim_msac_expt = false;
include_expts = {'rls.seXFaRC','rls.seXFaXFrRC','rls.orXFaXsl'};
if any(strcmp(Expt_name,{'G095','M266','M270'}))
    is_sim_msac_expt = true;
    % elseif any(strcmp(Expt_name,TBT_expts))
    %     is_TBT_expt = true;
end
if Expt_num >= 275
    is_TBT_expt = true;
end

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

cur_block_set = find(included_type & expt_ntrials' >= 10);
if strcmp(Expt_name,'G081')
    expt_has_ds = (expt_ijump==0)';
end
expt_has_ds(isnan(expt_has_ds)) = 0;

cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];
if length(unique(expt_dds(cur_block_set))) > 1
    fprintf('Warning, multiple dds detected!\n');
    main_dds = mode(expt_dds(cur_block_set));
    cur_block_set(expt_dds(cur_block_set) ~= main_dds) = [];
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

all_stim_times = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_ttill_end = [];
all_blockvec = [];
all_trialvec = [];
all_trial_wi = [];
all_frame_dur = [];
all_bar_or = [];
all_sac_dir = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
cur_spkind_offset = 0;
cur_toffset = 0;
trial_toffset = zeros(length(cur_block_set),1);

all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
for ee = 1:length(cur_block_set);
    cur_block = cur_block_set(ee);
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
    end
    
    use_trials = find(trial_durs >= min_trial_dur);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    
    if isfield(Expts{cur_block}.Trials(1),'Fr')
        trial_Fr = [Expts{cur_block}.Trials(:).Fr];
        all_frame_dur = [all_frame_dur; trial_Fr(use_trials)'];
    elseif isfield(Expts{cur_block}.Trials(1),'sl')
        trial_sl = [Expts{cur_block}.Trials(:).sl];
        trial_sl(trial_sl == 0) = 1;
        all_frame_dur = [all_frame_dur; trial_sl(use_trials)'];
    else
        all_frame_dur = [all_frame_dur; ones(length(use_trials),1)*expt_Fr(cur_block)];
    end
    if isfield(Expts{cur_block}.Trials(1),'or')
        trial_or = [Expts{cur_block}.Trials(:).or];
        all_bar_or = [all_bar_or; trial_or(use_trials)'];
    else
        all_bar_or = [all_bar_or; ones(length(use_trials),1)*expt_bar_ori(cur_block)];
    end
    if isfield(Expts{cur_block}.Trials(1),'Fa')
        trial_Fa = [Expts{cur_block}.Trials(:).Fa];
        all_sac_dir = [all_sac_dir; trial_Fa(use_trials)'];
    else
        all_sac_dir = [all_sac_dir; ones(length(use_trials),1)*expt_sac_dir(cur_block)];
    end
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        if length(cur_stim_times) == 1
            cur_stim_times = trial_start_times(use_trials(tt)):1/stim_fs:trial_end_times(use_trials(tt));
        end
        cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        cur_ttill_end = trial_end_times(use_trials(tt)) - cur_t_axis;
        
        all_stim_times = [all_stim_times; cur_stim_times'];
        all_t_axis = [all_t_axis; cur_t_axis];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
        all_tsince_start = [all_tsince_start; cur_tsince_start];
        all_ttill_end = [all_ttill_end; cur_ttill_end];
        all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    end
    trial_cnt = trial_cnt + n_trials;
end
all_sac_dir = mod(all_sac_dir,180);


%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = n_probes;
clust_params.exclude_adjacent = false;
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

%normalize firing rates within trials
n_trials = max(unique(all_trialvec));
for ii = 1:n_trials
    cur_trial_inds = find(all_trialvec == ii);
    if ~isempty(cur_trial_inds)
        all_sua_rate_norm(cur_trial_inds,:) = bsxfun(@rdivide,all_sua_rate(cur_trial_inds,:),...
            mean(all_sua_rate(cur_trial_inds,:)));
        all_mua_rate_norm(cur_trial_inds,:) = bsxfun(@rdivide,all_mua_rate(cur_trial_inds,:),...
            mean(all_mua_rate(cur_trial_inds,:)));
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
micro_set(ismember(micro_set,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'

%saccade amplitude along parallel axis
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);
sac_deltaY = sac_postpos(2,:) - sac_prepos(2,:);

sac_durs = [saccades(:).duration];
gsac_set = find((abs(sac_deltaX) > gsac_thresh | abs(sac_deltaY) > gsac_thresh) & ~used_is_blink' & sac_durs <= max_gsac_dur);

p100_trials = find(all_frame_dur == 1 & all_sac_dir == all_bar_or);
p30_trials = find(all_frame_dur == 3 & all_sac_dir == all_bar_or);
o100_trials = find(all_frame_dur == 1 & all_sac_dir ~= all_bar_or);
o30_trials = find(all_frame_dur == 3 & all_sac_dir ~= all_bar_or);

gsac_start_inds = sac_start_inds(gsac_set);
msac_start_inds = sac_start_inds(micro_set);

p100_gsacs = gsac_start_inds(ismember(all_trialvec(gsac_start_inds),p100_trials));
p100_msacs = msac_start_inds(ismember(all_trialvec(msac_start_inds),p100_trials));
p30_gsacs = gsac_start_inds(ismember(all_trialvec(gsac_start_inds),p30_trials));
p30_msacs = msac_start_inds(ismember(all_trialvec(msac_start_inds),p30_trials));
o100_gsacs = gsac_start_inds(ismember(all_trialvec(gsac_start_inds),o100_trials));
o100_msacs = msac_start_inds(ismember(all_trialvec(msac_start_inds),o100_trials));
o30_gsacs = gsac_start_inds(ismember(all_trialvec(gsac_start_inds),o30_trials));
o30_msacs = msac_start_inds(ismember(all_trialvec(msac_start_inds),o30_trials));

%% COMPUTE TRIG AVGS FOR MUA
nboot = [];
%set trial numbers to Inf so they don't get included in trig averaging
used_trialvec = ones(size(all_trialvec))*Inf;
used_trialvec(used_inds) = all_trialvec(used_inds);

%general averages
[mua_data.p100_gsac_avgs,lags] = get_event_trig_avg_v3(all_mua_rate_norm,p100_gsacs,backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.p30_gsac_avgs,lags] = get_event_trig_avg_v3(all_mua_rate_norm,p30_gsacs,backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.o100_gsac_avgs,lags] = get_event_trig_avg_v3(all_mua_rate_norm,o100_gsacs,backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.o30_gsac_avgs,lags] = get_event_trig_avg_v3(all_mua_rate_norm,o30_gsacs,backlag,forwardlag,nboot,used_trialvec,0);

mua_data.p100_gsac_avgs(:,16) = nan;    
mua_data.p30_gsac_avgs(:,16) = nan;    
mua_data.o100_gsac_avgs(:,16) = nan;    
mua_data.o30_gsac_avgs(:,16) = nan;    

mu_avg_rates = mean(all_binned_mua(used_inds,:))/dt;

%%
sname = 'par_orth_trig_avgs';
cd(save_dir)
save(sname,'mua_data','lags','dt','mu_avg_rates');
%%

min_MUA_rate = 25;
use_MUs = find(mu_avg_rates >= min_MUA_rate);

xl = [-0.15 0.4];
yl = [0.7 1.4];

f = figure(); hold on
h1=shadedErrorBar(lags*dt,nanmean(mua_data.p100_gsac_avgs(:,use_MUs),2),nanstd(mua_data.p100_gsac_avgs(:,use_MUs),[],2)/sqrt(length(use_MUs)),{'color','b'});
h3=shadedErrorBar(lags*dt,nanmean(mua_data.o100_gsac_avgs(:,use_MUs),2),nanstd(mua_data.o100_gsac_avgs(:,use_MUs),[],2)/sqrt(length(use_MUs)),{'color','k'});
xlabel('Time (s)');
ylabel('Relative Rate');
% legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'P-100','P-30','O-100','O-30'});
legend([h1.mainLine h3.mainLine],{'P-100','O-100'});
xlim(xl);
ylim(yl);
line([0 0],yl,'color','k','linestyle','--');
line(xl,[1 1],'color','k','linestyle','--');


fig_width = 3.5; rel_height = 0.8;
figufy(f);
fname = [fig_dir 'Par_vs_orth_MUA.pdf'];
exportfig(f,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f);

%%
p100_avgs = mua_data.p100_gsac_avgs(:,use_MUs)';
o100_avgs = mua_data.o100_gsac_avgs(:,use_MUs)';

search_range = [0 0.2];
[p_Sfact,p_inhtime] = get_tavg_peaks(-(p100_avgs-1),lags*dt,search_range);
[o_Sfact,o_inhtime] = get_tavg_peaks(-(o100_avgs-1),lags*dt,search_range);

search_range = [0.1 0.3];
[p_Efact,p_exctime] = get_tavg_peaks(p100_avgs-1,lags*dt,search_range);
[o_Efact,o_exctime] = get_tavg_peaks(o100_avgs-1,lags*dt,search_range);

