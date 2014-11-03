% clear all

addpath('~/James_scripts/bruce/eye_tracking/');
addpath('~/James_scripts/bruce/processing/');

global Expt_name bar_ori use_LOOXV

% Expt_name = 'M296';
% bar_ori = 90; %bar orientation to use (only for UA recs)
% use_LOOXV = 1;

recompute_init_mods = 0; %use existing initial models?
use_measured_pos = 3; %1 for init with coils, 2 for init with trial-sub coils, 3 for random init, 0 for init to perfect fixation
use_sac_kerns = 1; %use sac-modulation kernels

%%

Expt_num = str2num(Expt_name(2:end));
if Expt_name(1) == 'M'
    rec_type = 'LP';
elseif Expt_name(1) == 'G'
    rec_type = 'UA';
end

if strcmp(rec_type,'LP')
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
end

if Expt_num > 280
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end

cd(data_dir);

if strcmp(rec_type,'LP')
    if Expt_num >= 275
        rpt_seed = 1001; %M275 M277 M281
    else
        rpt_seed = 1e4; %m270 and 266
    end
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
    good_coils = [1 1]; %which coils are usable
    use_coils = [1 1]; %[L R] Use info from coils?
    n_probes = 24;
    use_nPix = 32;
elseif strcmp(rec_type,'UA')
    rpt_seed = nan;
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
    good_coils = [1 0]; %which coils are usable
    use_coils = [0 0]; %[L R] Use info from coils?
    n_probes = 96;
    use_nPix = 16;
end

load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

mod_data_name = 'full_eyetrack_initmods';
old_anal_name = 'full_eyetrack';

%if using coil initialization
if use_measured_pos == 1
    mod_data_name = [mod_data_name '_Cinit'];
    old_anal_name = [old_anal_name '_Cinit'];
end
%if using coil init with trial-sub
if use_measured_pos == 2
    mod_data_name = [mod_data_name '_CPinit'];
    old_anal_name = [old_anal_name '_CPinit'];
end
if use_measured_pos == 3
    mod_data_name = [mod_data_name '_Rinit'];
    old_anal_name = [old_anal_name '_Rinit'];
end
%if using coil info
if any(use_coils > 0)
    old_anal_name = [old_anal_name '_Cprior'];
end

hr_anal_name = [old_anal_name '_hres'];
hr_mod_name = [mod_data_name 'hres'];

mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];
old_anal_name = [old_anal_name sprintf('_ori%d',bar_ori)];
hr_anal_name = [hr_anal_name sprintf('_ori%d',bar_ori)];
hr_mod_name = [hr_mod_name sprintf('_ori%d',bar_ori)];

%%
%dont fit stim models using these blocks
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
        ignore_blocks = [28];
end

%problem with M270 where vd was wrong, need this correction factor to get
%correct units
if Expt_num==270
    scale_fac = 1.72;
else
    scale_fac = 1;
end


%%
if ~isnan(rpt_seed)
    xv_type = 'rpt';
    xv_frac = nan;
else
    xv_type = 'uni';
    xv_frac = 0.2;
end

flen = 12;
old_spatial_usfac = 2;

%these recs have larger bar widths
if ismember(Expt_num,[287 289 294])
    use_nPix = 15;
    old_spatial_usfac = 4;
elseif ismember(Expt_num,[296 297])
    use_nPix = 22;
    old_spatial_usfac = 2;
end

spatial_usfac = old_spatial_usfac*2;

n_drift_inf_it = 1; %3

fix_prior_sigma = 0.15;
fix_noise_sigma = 0.1;

%this sets the posterior sigma to be ~0.004 in all cases
drift_noise_sigma = 0;
if all(use_coils)
    drift_prior_sigma = sqrt(0.004^2*3);
    drift_noise_sigma = sqrt(0.004^2*3);
elseif any(use_coils)
    drift_prior_sigma = sqrt(0.004^2*2);
    drift_noise_sigma = sqrt(0.004^2*2);
else
    drift_prior_sigma = 0.004;
end
drift_jump_sigma = 0.075; %0.05 start
drift_dsf = 2;

min_trial_dur = 0.75;

stim_fs = 100; %in Hz
dt = 0.01;
Fr = 1;

full_nPix=36;
switch Expt_num
    case 270
        full_nPix=32;
    case  287
        full_nPix = 22;
    case 289
        full_nPix = 22;
    case 294
        full_nPix = 20;
end

%exclude data at beginning and end of each trial
beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

%saccade kernel time axis
sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

stim_params = NIMcreate_stim_params([flen full_nPix],dt);

%% SELECT BLOCKS FOR ANALYSIS
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa','rls.FaXwi','rls.FaXwiXimi','rls.AllSacB'};
expt_names = cell(1,length(Expts));
expt_dds = nan(1,length(Expts));
expt_bar_ori = nan(1,length(Expts));
expt_sac_dir = nan(1,length(Expts));
expt_Fr = nan(1,length(Expts));
expt_imback = nan(1,length(Expts));
included_type = false(1,length(Expts));
for ii = 1:length(Expts)
    if ~isempty(Expts{ii})
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    end
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;

if strcmp(rec_type,'LP')
    expt_bar_ori(expt_bar_ori > 360) = bar_ori;
end

cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);
cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];
if length(unique(expt_dds(cur_block_set))) > 1
    fprintf('Warning, multiple dds detected!\n');
    main_dds = mode(expt_dds(cur_block_set));
    cur_block_set(expt_dds(cur_block_set) ~= main_dds) = [];
end

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

n_blocks = length(cur_block_set);

all_nfs = cellfun(@(x) x.Stimvals.nf,Expts(cur_block_set));
if length(unique(all_nfs)) > 1
    fprintf('Warning, multiple different nfs detected: %.4f\n',all_nfs);
end
%%
all_dws = cellfun(@(x) x.Stimvals.dw,Expts(cur_block_set));
base_sp_dx = mode(all_dws);
if length(unique(all_dws)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac/scale_fac; %model dx in deg

max_shift = round(0.8475/sp_dx);
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

max_Dshift = round(0.452/sp_dx);
Dshifts = -max_Dshift:dshift:max_Dshift;
n_Dshifts = length(Dshifts);
zero_Dframe = find(Dshifts==0);

max_Tshift = max_shift + max_Dshift;
Tshifts = -max_Tshift:dshift:max_Tshift;

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_wi = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_trial_rptframes = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
for ee = 1:n_blocks;
    if ismember(ee,grayback_gs_expts)
        fprintf('Expt %d Block %d of %d; grayback GS, ori:%d\n',Expt_num,ee,n_blocks,expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Expt %d Block %d of %d; imback GS, ori:%d\n',Expt_num,ee,n_blocks,expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Expt %d Block %d of %d; SimSac, ori:%d\n',Expt_num,ee,n_blocks,expt_bar_ori(cur_block_set(ee)));
    else
        fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,n_blocks);
    end
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
    rpt_trials = false;
    if length(un_ids) < length(trial_ids)
        rpt_trials = true;
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
    
    if isfield(Expts{cur_block}.Trials,'wi')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    else
        all_trial_wi = cat(1,all_trial_wi,nan(length(use_trials),1));
    end
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    %     fname = sprintf('%s/stims/Expt%d_oldStim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
    if buffer_pix == -1
        for ii = 1:length(left_stim_mats)
            left_stim_mats{ii} = [zeros(size(left_stim_mats{ii},1),1) left_stim_mats{ii} zeros(size(left_stim_mats{ii},1),1)];
        end
        buffer_pix = 0;
    end
    cur_use_pix = (1:full_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    cur_nrpt_frames = zeros(n_trials,1);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if isfield(Expts{cur_block}.Trials(use_trials(tt)),'rptframes')
            cur_nrpt_frames(tt) = length(Expts{cur_block}.Trials(use_trials(tt)).rptframes);
        end
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
                cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
            end
        end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times+cur_toffset < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
    all_trial_rptframes = [all_trial_rptframes; cur_nrpt_frames];
    
    %need to keep track of block time offsets for LP recordings
    if strcmp(rec_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

%%
full_nPix_us = spatial_usfac*full_nPix;
if spatial_usfac > 1
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:spatial_usfac
            all_stimmat_up(:,spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
        end
    end
elseif spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
end
stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);
all_Xmat = create_time_embedding(all_stim_mat,stim_params);
all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);

%% select submatrix with central pixels

%pick out the predictors (of the time-embedded stimulus-matrices)
%corresponding to the central pixels we're using in the models
[Xinds,~] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix)); %use this set of dimensions in the Xmatrix

%repeat for up-sampled versions of the Xmatrix
[Xinds_up,~] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
cur_use_pix = (1/spatial_usfac:1/spatial_usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

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
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);

%for G093 use only data where stripe width is AT LEAST 2 deg
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi >= un_wi_vals(2));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end
NT = length(used_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_stim_times),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% PROCESS EYE TRACKING DATA
if isfield(Expts{cur_block_set(1)}.Header,'exptno')
    em_block_nums = cellfun(@(X) X.Header.exptno,Expts(cur_block_set),'uniformoutput',1); %block numbering for EM/LFP data sometimes isnt aligned with Expts struct
else
    em_block_nums = cur_block_set;
end

% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,em_block_nums,Expt_name,trial_toffset,good_coils);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
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

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);
used_is_blink = is_blink(used_saccade_set);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro & ~used_is_blink');
micro_sacs = find(is_micro & ~used_is_blink');
sac_durs = [saccades(:).duration];

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

if use_sac_kerns
    Xsac = zeros(NT,n_sac_bins);
    Xmsac = zeros(NT,n_sac_bins);
    for ii = 1:n_sac_bins
        cur_sac_target = saccade_start_inds(big_sacs) + sac_bincents(ii);
        uu = find(cur_sac_target > 1 & cur_sac_target < NT);
        cur_sac_target = cur_sac_target(uu);
        cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(big_sacs(uu))) = [];
        Xsac(cur_sac_target,ii) = 1;
        
        cur_sac_target = saccade_start_inds(micro_sacs) + sac_bincents(ii);
        uu = find(cur_sac_target > 1 & cur_sac_target < NT);
        cur_sac_target = cur_sac_target(uu);
        cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(micro_sacs(uu))) = [];
        Xmsac(cur_sac_target,ii) = 1;
    end
end

%% DEFINE FIXATION POINTS
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink));
n_fixs = length(fix_start_inds);

%push the effects of saccades forward in time
sac_shift = round(0.05/dt);
pfix_start_inds = fix_start_inds;
pfix_stop_inds = fix_stop_inds;
for i = 1:length(fix_start_inds)
    next_trial = trial_start_inds(find(trial_start_inds >= fix_start_inds(i),1,'first'));
    if next_trial > fix_start_inds(i) + sac_shift
        pfix_start_inds(i) = fix_start_inds(i) + sac_shift;
    end
    next_trial = trial_start_inds(find(trial_start_inds >= fix_stop_inds(i),1,'first'));
    if next_trial > fix_stop_inds(i) + sac_shift
        pfix_stop_inds(i) = fix_stop_inds(i) + sac_shift;
    end
end

%for all times within the (forward-projected) saccade, use 'prior'
%state-transition model
use_prior = zeros(NT,1);
for i = 1:n_fixs-1
    use_prior((pfix_stop_inds(i)+1):pfix_start_inds(i+1)) = 1;
end
fix_ids = nan(NT,1);
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
    cur_inds = pfix_start_inds(ii):(pfix_stop_inds(ii));
    pfix_ids(cur_inds) = ii;
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

if strcmp(xv_type,'rpt')
    xv_trials = find(all_trial_Se==rpt_seed);
    n_xv_trials = length(xv_trials);
else
    n_xv_trials = round(xv_frac*nuse_trials);
    xv_trials = randperm(nuse_trials);
    xv_trials(n_xv_trials+1:end) = [];
    xv_trials = use_trials(xv_trials);
end
tr_trials = setdiff(use_trials,xv_trials);
n_tr_trials = length(tr_trials);
fprintf('Initializing models with %d training trials and %d xval trials\n',n_tr_trials,n_xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

%%
% make Robs_mat
tot_sus = size(all_binned_sua,2);
Robs_mat = nan(length(used_inds),n_probes + tot_sus);
for ss = 1:size(Robs_mat,2)
    if ss > n_probes
        Robs_mat(:,ss) = all_binned_sua(used_inds,ss-n_probes);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end


%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
%push the effects of saccades forward in time
trial_ids = nan(NT,1);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
end

%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = spdiags( ones(full_nPix_us,1), -shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    shift_mat{xx} = temp(:,use_kInds_up);
end
Dshift_mat = cell(n_Dshifts,1);
for xx = 1:n_Dshifts
    temp = spdiags( ones(full_nPix_us,1), -Dshifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    Dshift_mat{xx} = temp(:,use_kInds_up);
end

%shift matrices for images (non-time-embedded)
n_Tshifts = length(Tshifts);
Tshift_mat = cell(n_Tshifts,1);
for xx = 1:n_Tshifts
    Tshift_mat{xx} = spdiags( ones(full_nPix_us,1), -Tshifts(xx), full_nPix_us, full_nPix_us);
end

%%
cd(anal_dir)

fprintf('Loading pre-computed initial models\n');
load(mod_data_name);
load(old_anal_name,'dit_mods','et_tr_set','it_fix_post_mean','drift_post_mean','it_fix_post_std','drift_post_std');
old_best_mods = dit_mods{end};
tr_set = et_tr_set;

up_fac = spatial_usfac/old_spatial_usfac;
best_fix_cor = it_fix_post_mean(end,:)*up_fac;
best_fix_std = it_fix_post_std(end,:)*up_fac;
best_drift_cor = drift_post_mean(end,:)*up_fac;
best_drift_std = drift_post_std(end,:)*up_fac;

clear it_fix_post_mean drift_post_mean dit_mods


%%
[best_tot_corr,best_tot_std] = construct_eye_position(best_fix_cor,best_fix_std,...
    best_drift_cor,best_drift_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);

best_tot_corr = interp1(find(~isnan(fix_ids)),best_tot_corr(~isnan(fix_ids)),1:NT);
best_tot_corr(isnan(best_tot_corr)) = 0;
all_post_cor = round(best_tot_corr) + max_Tshift + 1;

%RECOMPUTE XMAT
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
end
all_Xmat_cor = create_time_embedding(all_shift_stimmat_up,stim_params_us);

%% SELECT USABLE UNITS AND make Robs_mat

n_tr_chs = length(tr_set);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = Robs_mat(:,tr_set);

%don't use separate xv set for eye-tracking
if strcmp(rpt_trials,'uni')
    tr_inds = sort([tr_inds; xv_inds]);
end

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
fin_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
fin_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
if use_sac_kerns
    fin_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
    fin_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);
end

null_stim_params = fin_stim_params(2:end);

block_L2 = 1;
silent = 0;
sac_d2t = 100;

if expt_dds(cur_block_set(1)) == 67
    base_lambda_d2XT = 60;
    base_lambda_L1 = 10;
else
    base_lambda_d2XT = 20;
    base_lambda_L1 = 5;
end
init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-9;

sac_reg_params = NMMcreate_reg_params('lambda_d2T',sac_d2t);
if use_sac_kerns
    null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);
else
    null_reg_params = NMMcreate_reg_params('lambda_L2',block_L2);
end

n_squared_filts = 2;
mod_signs = ones(1,n_squared_filts+2);
%if using more than 2 quadratic filters, take the additional ones as
%suppressive
if n_squared_filts > 2
    mod_signs(n_squared_filts+1) = -1;
end
NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
init_d2XT = [ones(n_squared_filts+1,1); 0;];
init_L2 = [zeros(n_squared_filts+1,1); block_L2];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2);
init_Xtargs = [ones(n_squared_filts+1,1); 2];

init_filts = cell(length(mod_signs),1);

if ~exist(['./' hr_mod_name '.mat'],'file') || recompute_init_mods == 1
    
    for ss = 1:n_tr_chs
        fprintf('Computing base LLs for Unit %d of %d\n',ss,n_tr_chs);
        cur_tr_inds = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
        tr_NT = length(cur_tr_inds);
        Robs = Robs_mat(cur_tr_inds,ss);
        null_mod = all_nullmod(tr_set(ss));
        if ~isempty(cur_tr_inds) && nansum(Robs) > 0
            
            tr_X{1} = all_Xmat_cor(used_inds(cur_tr_inds),use_kInds_up);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
            end
            
            %spatial up-sampling of filter estimates
            base_filts = reshape([old_best_mods(tr_set(ss)).mods(find(init_Xtargs == 1)).filtK],[flen use_nPix*old_spatial_usfac n_squared_filts+1]);
            base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
            for ii = 1:use_nPix*old_spatial_usfac
                for jj = 1:up_fac
                    base_filts_up(:,up_fac*(ii-1)+jj,:) = 1/up_fac*base_filts(:,ii,:);
                end
            end
            base_filts_up = reshape(base_filts_up,use_nPix_us*flen,n_squared_filts+1);
            
            init_filts{end} = old_best_mods(tr_set(ss)).mods(find(init_Xtargs==2)).filtK;
            for ii = 1:n_squared_filts+1
                init_filts{ii} = base_filts_up(:,ii);
            end
            gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
            vgint = var(gint)'; sgint = std(gint)';
            vgint(vgint == 0) = 1; sgint(sgint == 0) = 1;
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./vgint);
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./sgint);
            if use_sac_kerns
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,old_best_mods(tr_set(ss)).mods(end-1).filtK); %sac term
                gqm2.mods(end).reg_params = sac_reg_params;
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,old_best_mods(tr_set(ss)).mods(end).filtK); %nsac term
                gqm2.mods(end).reg_params = sac_reg_params;
            end
            all_mod_fits(tr_set(ss)) = gqm2;
            all_mod_fits(tr_set(ss)).spk_NL_params(1) = old_best_mods(tr_set(ss)).spk_NL_params(1);
            all_mod_fits(tr_set(ss)) = NMMfit_filters(all_mod_fits(tr_set(ss)),Robs,tr_X,[],[],silent);
            
            all_mod_fits_withspkNL(tr_set(ss)) = NMMfit_logexp_spkNL(all_mod_fits(tr_set(ss)),Robs,tr_X);
            
            [LL,~,pred_rate] = NMMmodel_eval(all_mod_fits_withspkNL(tr_set(ss)),Robs,tr_X);
            [null_LL(tr_set(ss)),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            all_mod_LLimp(tr_set(ss)) = (LL-null_LL(tr_set(ss)))/log(2);
        else
            all_mod_LLimp(tr_set(ss)) = nan;
        end
    end
    save(hr_mod_name,'all_mod*','null_LL');
    
else
    fprintf('Loading pre-computed initial models\n');
    load(hr_mod_name);
end

clear tr_X
%%
measured_eyepos = [corrected_eye_vals_interp(:,2) corrected_eye_vals_interp(:,4)];
max_sim_pos = max_shift*sp_dx;
measured_eyepos(measured_eyepos > max_sim_pos) = max_sim_pos;
measured_eyepos(measured_eyepos < -max_sim_pos) = -max_sim_pos;
measured_eyepos(isnan(measured_eyepos)) = 0;

%smooth out fast transients in eye signal
eye_smooth_sig = round(0.025/dt);
interp_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > eye_smooth_sig*5;
        measured_eyepos(used_inds(cur_inds),1) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),1),eye_smooth_sig,2);
        measured_eyepos(used_inds(cur_inds),2) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),2),eye_smooth_sig,2);
    end
    interp_inds = [interp_inds; cur_inds'];
end
interp_inds = unique(interp_inds);
measured_eyepos(used_inds,:) = interp1(used_inds(interp_inds),measured_eyepos(used_inds(interp_inds),:),used_inds);
measured_eyepos(isnan(measured_eyepos)) = 0;
measured_eyepos = measured_eyepos(used_inds,:);

measured_drift = nan(length(used_inds),2);
measured_fix_avg = nan(n_fixs,2);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    measured_fix_avg(ii,:) = median(measured_eyepos(cur_inds,:));
    measured_drift(cur_inds,:) = bsxfun(@minus,measured_eyepos(cur_inds,:),measured_fix_avg(ii,:));
end

back_measured_drift = nan(size(measured_drift));
%back-project measured drift
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    cur_pinds = cur_inds + sac_shift;
    if length(cur_inds) > sac_shift
        back_measured_drift(cur_inds(sac_shift+1:end),:) = measured_drift(cur_inds(1:(end-sac_shift)),:);
    end
end
measured_drift = back_measured_drift;
clear back_measured_drift

measured_drift = [zeros(1,2); diff(measured_drift)];

%if using both coils
if use_coils(1) == 1 && use_coils(2) == 1
    
    measured_fix_deltas = nan(n_fixs,1);
    measured_fix_deltas(2:end) = mean(diff(measured_fix_avg),2);
    fix_noise_sigma = fix_noise_sigma/sqrt(2); %noise variance decreases by factor of sqrt(2) with 2 independent measures
    post_drift_var = 1/(2/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    post_mean_drift = post_drift_var*2*mean(measured_drift,2)/drift_noise_sigma^2;
    
    %if using left coil only
elseif use_coils(1) == 1
    measured_fix_deltas = nan(n_fixs,1);
    measured_fix_deltas(2:end) = diff(measured_fix_avg(:,1));
    post_drift_var = 1/(1/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    post_mean_drift = post_drift_var*measured_drift(:,1)/drift_noise_sigma^2;
    
    %if using right coil only
elseif use_coils(2) == 1
    measured_fix_deltas = nan(n_fixs,1);
    measured_fix_deltas(2:end) = diff(measured_fix_avg(:,2));
    post_drift_var = 1/(1/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    post_mean_drift = post_drift_var*measured_drift(:,2)/drift_noise_sigma^2;
    
else
    post_drift_var = drift_prior_sigma^2;
    post_mean_drift = zeros(NT,1);
end

post_drift_sigma = sqrt(post_drift_var);
% post_mean_drift(isnan(post_mean_drift)) = 0;

fix_prior = -(shifts*sp_dx).^2./(2*fix_prior_sigma^2);
fix_prior = fix_prior - logsumexp(fix_prior);

drift_jump_prior = -(Dshifts*sp_dx).^2./(2*drift_jump_sigma^2);
drift_jump_prior = drift_jump_prior - logsumexp(drift_jump_prior);

cdist = squareform(pdist(Dshifts'*sp_dx));
base_lA = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
base_lA = bsxfun(@minus,base_lA,logsumexp(base_lA,2)); %normalize

%% set of units to perform LOO xval on
if use_LOOXV == 2
    loo_set = 1:length(tr_set);
elseif use_LOOXV == 1
    loo_set = find(all_mod_SU(tr_set) > 0);
else
    loo_set = [];
end


%%
%back-project saccade-times
all_fix_post_mean_cor = nan(NT,1);
all_fix_post_mean_cor(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;

all_fix_post_mean_cor = round(all_fix_post_mean_cor) + max_Tshift + 1;
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_fix_post_mean_cor(i)};
end
all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);

%% NOW INFER DRIFT CORRECTIONS
dit_mods = all_mod_fits;
dit_mods_spkNL = all_mod_fits_withspkNL;
dit_LLimp = all_mod_LLimp;
dit_R2 = all_mod_R2;
for xv = 1:length(loo_set)
    dit_mods_LOO{xv} = all_mod_fits;
    dit_mods_spkNL_LOO{xv} = all_mod_fits_withspkNL;
end

%% PREPROCESS MODEL COMPONENTS
filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
lin_kerns = nan(n_tr_chs,n_blocks);
if use_sac_kerns
    sac_kerns = nan(n_tr_chs,n_sac_bins);
    msac_kerns = nan(n_tr_chs,n_sac_bins);
end
mod_spkNL_params = nan(n_tr_chs,3);
for ss = 1:n_tr_chs
    cur_Xtargs = [dit_mods(tr_set(ss)).mods(:).Xtarget];
    cur_k = [dit_mods(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = dit_mods_spkNL(tr_set(ss)).spk_NL_params(1:3);
    lin_kerns(ss,:) = dit_mods(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
    if use_sac_kerns
        sac_kerns(ss,:) = dit_mods(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
        msac_kerns(ss,:) = dit_mods(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
    end
end
filt_bank = permute(filt_bank,[2 1 3]);

%indicator predictions
block_out = Xblock(used_inds,:)*lin_kerns';
if use_sac_kerns
    sac_out = Xsac*sac_kerns';
    msac_out = Xmsac*msac_kerns';
end

%% ESTIMATE LL for each shift in each stimulus frame

frame_LLs = nan(NT,n_Dshifts);
for xx = 1:length(Dshifts)
    fprintf('Dshift %d of %d\n',xx,n_Dshifts);
    cur_stim_shift = all_Xmat_up_fixcor*Dshift_mat{xx};
    
    %outputs of stimulus models at current X-matrix shift
    gfuns = ones(NT,n_tr_chs);
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
    gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1));
    for ff = 2:(n_squared_filts+1)
        gfuns = gfuns + mod_signs(ff)*(cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
    end
    
    %add contributions from extra lin kernels
    gfuns = gfuns + block_out;
    if use_sac_kerns
        gfuns = gfuns + sac_out + msac_out;
    end
    
    %incorporate beta
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
    
    %handle numerical overflow with log(1+exp)
    too_large = gfuns > 50;
    pred_rate = log(1+exp(gfuns));
    pred_rate(too_large) = gfuns(too_large);
    
    %incorporate alpha
    pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
    
    %enforce min predicted rate
    pred_rate(pred_rate < 1e-50) = 1e-50;
    
    frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
end

%% INFER DRIFT CORRECTIONS

lgamma = nan(NT,n_Dshifts);
for ff = 1:n_fixs
    if mod(ff,100)==0
        fprintf('Fixation %d of %d\n',ff,n_fixs);
    end
    
    tset = find(pfix_ids==ff)';
    ntset = length(tset);
    if ntset > drift_dsf
        nt_pts = ceil(ntset/drift_dsf);
        tset_inds = 1+floor((0:(ntset-1))/drift_dsf);
        talpha=zeros(nt_pts,n_Dshifts);
        tbeta = zeros(nt_pts,n_Dshifts);
        
        tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf);
        tpt_loc(end) = ntset;
        
        cur_drift_mean = post_mean_drift(tset);
        cur_LL_set = frame_LLs(tset,:);
        if mod(ntset,drift_dsf) ~= 0
            dangling_pts = nt_pts*drift_dsf-ntset;
            cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_Dshifts));
            cur_drift_mean = cat(1,cur_drift_mean,nan(dangling_pts,1));
        end
        cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_Dshifts]);
        cur_LL_set = squeeze(sum(cur_LL_set,1));
        
        cur_drift_mean = drift_dsf*nanmean(reshape(cur_drift_mean,[drift_dsf nt_pts]));
        
        if all(use_coils==0)
            cur_lA = repmat(base_lA,[1 1 nt_pts]);
        else
            cur_lA = nan(n_Dshifts,n_Dshifts,nt_pts);
            for iii = 1:nt_pts
                if ~isnan(cur_drift_mean(iii))
                    cdist = pdist2(Dshifts'*sp_dx + cur_drift_mean(iii),Dshifts'*sp_dx);
                    cur_lA(:,:,iii) = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
                else
                    cur_lA(:,:,iii) = base_lA;
                end
            end
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
        end
        
        talpha(1,:) = drift_jump_prior + cur_LL_set(1,:);
        for t = 2:nt_pts
            talpha(t,:) = logmulexp(talpha(t-1,:),cur_lA(:,:,t)) + cur_LL_set(t,:);
        end
        
        tbeta(end,:)=log(ones(1,n_Dshifts));
        for t = (nt_pts-1):-1:1
            lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
            tbeta(t,:) = logmulexp(lf1,cur_lA(:,:,t+1)');
        end
        temp_gamma = talpha + tbeta;
        
        if drift_dsf > 1
            if nt_pts > 1
                int_gamma = interp1(tpt_loc,temp_gamma,1:ntset);
                lgamma(tset,:) = int_gamma;
            else
                lgamma(tset,:) = repmat(temp_gamma,ntset,1);
            end
        else
            lgamma(tset,:) = temp_gamma;
        end
    end
end
lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
gamma = exp(lgamma);
drift_post_mean = sum(bsxfun(@times,gamma,Dshifts),2);
drift_post_std = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2) - drift_post_mean.^2);

drift_post_mean(isnan(drift_post_mean)) = 0;
drift_post_std(isnan(drift_post_std)) = 0;

drift_post_mean = interp1(find(~isnan(pfix_ids)),drift_post_mean(~isnan(pfix_ids)),1:NT);
drift_post_std = interp1(find(~isnan(pfix_ids)),drift_post_std(~isnan(pfix_ids)),1:NT);
drift_post_mean(isnan(drift_post_mean)) = 0;

if use_LOOXV > 0
    for xv = 1:length(loo_set)
        fprintf('Inferring drift corrections, XV %d of %d\n',xv,length(loo_set));
               
        %% PREPROCESS MODEL COMPONENTS
        cur_uset = setdiff(1:n_tr_chs,loo_set(xv));
        n_uset = length(cur_uset);
        filt_bank = zeros(n_uset,klen_us,n_squared_filts+1);
        lin_kerns = nan(n_uset,n_blocks);
        if use_sac_kerns
            sac_kerns = nan(n_uset,n_sac_bins);
            msac_kerns = nan(n_uset,n_sac_bins);
        end
        mod_spkNL_params = nan(n_uset,3);
        for ss = 1:n_uset
            cur_Xtargs = [dit_mods_LOO{xv}(tr_set(cur_uset(ss))).mods(:).Xtarget];
            cur_k = [dit_mods_LOO{xv}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 1).filtK];
            n_used_filts = size(cur_k,2);
            filt_bank(ss,:,1:n_used_filts) = cur_k;
            mod_spkNL_params(ss,:) = dit_mods_spkNL_LOO{xv}(tr_set(cur_uset(ss))).spk_NL_params(1:3);
            lin_kerns(ss,:) = dit_mods_LOO{xv}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 2).filtK;
            if use_sac_kerns
                sac_kerns(ss,:) = dit_mods_LOO{xv}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 3).filtK;
                msac_kerns(ss,:) = dit_mods_LOO{xv}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 4).filtK;
            end
        end
        filt_bank = permute(filt_bank,[2 1 3]);
        
        %indicator predictions
        block_out = Xblock(used_inds,:)*lin_kerns';
        if use_sac_kerns
            sac_out = Xsac*sac_kerns';
            msac_out = Xmsac*msac_kerns';
        end
        %% ESTIMATE LL for each shift in each stimulus frame
        frame_LLs = nan(NT,n_Dshifts);
        for xx = 1:length(Dshifts)
            fprintf('Dshift %d of %d\n',xx,n_Dshifts);
            cur_stim_shift = all_Xmat_up_fixcor*Dshift_mat{xx};
            
            %outputs of stimulus models at current X-matrix shift
            gfuns = ones(NT,n_uset);
            gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
            gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1));
            for ff = 2:(n_squared_filts+1)
                gfuns = gfuns + mod_signs(ff)*(cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
            end
            
            %add contributions from extra lin kernels
            gfuns = gfuns + block_out;
            if use_sac_kerns
                gfuns = gfuns + sac_out + msac_out;
            end
            
            %incorporate beta
            gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
            
            %handle numerical overflow with log(1+exp)
            too_large = gfuns > 50;
            pred_rate = log(1+exp(gfuns));
            pred_rate(too_large) = gfuns(too_large);
            
            %incorporate alpha
            pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
            
            %enforce min predicted rate
            pred_rate(pred_rate < 1e-50) = 1e-50;
            
            frame_LLs(:,xx) = squeeze(nansum(Robs_mat(:,cur_uset).*log(pred_rate) - pred_rate,2));
        end
        
        %% INFER MICRO-SAC SEQUENCE
        lgamma = nan(NT,n_Dshifts);
        for ff = 1:n_fixs
            if mod(ff,100)==0
                fprintf('Fixation %d of %d\n',ff,n_fixs);
            end
            
            %                 tset = find(fix_ids==ff)';
            tset = find(pfix_ids==ff)';
            ntset = length(tset);
            if ntset > drift_dsf
                nt_pts = ceil(ntset/drift_dsf);
                tset_inds = 1+floor((0:(ntset-1))/drift_dsf);
                talpha=zeros(nt_pts,n_Dshifts);
                tbeta = zeros(nt_pts,n_Dshifts);
                
                tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf);
                tpt_loc(end) = ntset;
                
                cur_drift_mean = post_mean_drift(tset);
                cur_LL_set = frame_LLs(tset,:);
                if mod(ntset,drift_dsf) ~= 0
                    dangling_pts = nt_pts*drift_dsf-ntset;
                    cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_Dshifts));
                    cur_drift_mean = cat(1,cur_drift_mean,nan(dangling_pts,1));
                end
                cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_Dshifts]);
                cur_LL_set = squeeze(sum(cur_LL_set,1));
                
                cur_drift_mean = drift_dsf*nanmean(reshape(cur_drift_mean,[drift_dsf nt_pts]));
                
                if all(use_coils==0)
                    cur_lA = repmat(base_lA,[1 1 nt_pts]);
                else
                    cur_lA = nan(n_Dshifts,n_Dshifts,nt_pts);
                    for iii = 1:nt_pts
                        if ~isnan(cur_drift_mean(iii))
                            cdist = pdist2(Dshifts'*sp_dx + cur_drift_mean(iii),Dshifts'*sp_dx);
                            cur_lA(:,:,iii) = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
                        else
                            cur_lA(:,:,iii) = base_lA;
                        end
                    end
                    cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
                end
                
                talpha(1,:) = drift_jump_prior + cur_LL_set(1,:);
                for t = 2:nt_pts
                    talpha(t,:) = logmulexp(talpha(t-1,:),cur_lA(:,:,t)) + cur_LL_set(t,:);
                end
                
                tbeta(end,:)=log(ones(1,n_Dshifts));
                for t = (nt_pts-1):-1:1
                    lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
                    tbeta(t,:) = logmulexp(lf1,cur_lA(:,:,t+1)');
                end
                temp_gamma = talpha + tbeta;
                
                if drift_dsf > 1
                    if nt_pts > 1
                        int_gamma = interp1(tpt_loc,temp_gamma,1:ntset);
                        lgamma(tset,:) = int_gamma;
                    else
                        lgamma(tset,:) = repmat(temp_gamma,ntset,1);
                    end
                else
                    lgamma(tset,:) = temp_gamma;
                end
            end
        end
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        
        gamma = exp(lgamma);
        drift_post_mean_LOO(xv,:) = sum(bsxfun(@times,gamma,Dshifts),2);
        drift_post_std_LOO(xv,:) = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2)) - drift_post_mean_LOO(xv,:)'.^2;
        
        drift_post_mean_LOO(xv,isnan(drift_post_mean_LOO(xv,:))) = 0;
        drift_post_std_LOO(xv,isnan(drift_post_std_LOO(xv,:))) = 0;
        
        drift_post_mean_LOO(xv,:) = interp1(find(~isnan(pfix_ids)),squeeze(drift_post_mean_LOO(xv,~isnan(pfix_ids))),1:NT);
        drift_post_std_LOO(xv,:) = interp1(find(~isnan(pfix_ids)),squeeze(drift_post_std_LOO(xv,~isnan(pfix_ids))),1:NT);
    end
end


%% SAVE EYE-TRACKING RESULTS
et_params = struct('beg_buffer',beg_buffer,'end_buffer',end_buffer,'min_trial_dur',min_trial_dur,'bar_ori',bar_ori,'good_coils',good_coils,...
    'use_nPix',use_nPix,'flen',flen,'dt',dt,'drift_jump_sigma',drift_jump_sigma,'drift_prior_sigma',drift_prior_sigma,...
    'fix_prior_sigma',fix_prior_sigma,'fix_noise_sigma',fix_noise_sigma,'drift_noise_sigma',drift_noise_sigma,...
    'drift_dsf',drift_dsf,'use_sac_kerns',use_sac_kerns,'shifts',shifts,...
    'use_measured_pos',use_measured_pos,'sac_bincents',sac_bincents,'spatial_usfac',spatial_usfac,'old_spatial_usfac',old_spatial_usfac,'sac_shift',sac_shift,'use_coils',use_coils,'sp_dx',sp_dx);

et_used_inds = used_inds;
et_tr_set = tr_set;
et_saccades = saccades;
et_is_blink = is_blink;
et_clust_data = Clust_data;
cd(anal_dir);
save(hr_anal_name,'best_fix_*','drift_post_*','fix_ids','dit_*','et_used_inds','et_tr_set','et_clust_data','et_saccades','et_is_blink','et_params');


%%
% close all
% n_trials = length(unique(all_trialvec));
% for tt = 1:n_trials
%     % for tt = [96 137 154 179 376 409]
%     uu = find(all_trialvec(used_inds) == tt);
%     if ~isempty(uu)
%         bt = all_t_axis(used_inds(uu(1)));
%         et = all_t_axis(used_inds(uu(end)));
%         dur = et-bt;
%         if dur > 3.5
%             hold on
% %             h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_fix_corr(uu),fin_fix_std(uu),{'color','m'});
%             h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','k'});
%             h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'r','linewidth',2);
%             h4=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),4),'k','linewidth',2);
%             %                 h4=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),4)-median(corrected_eye_vals_interp(used_inds(uu),4)),'color',[0.2 0.8 0.2],'linewidth',2);
%             %             plot(all_t_axis(used_inds(uu))-bt,nanmean(Robs_mat(uu,:),2)/5,'k');
%
%             %             legend([h1.mainLine h2.mainLine h3 h4],{'Fixation corrections','Drift corrections','Left eye','Right eye'})
%             xlim([0 dur]);
%             ylim([-0.5 0.5]);
%             xlabel('Time (s)','fontsize',10);
%             ylabel('Orthoganol position (deg)','fontsize',10);
%             title(sprintf('Trial %d',tt));
%             set(gca,'fontsize',8,'fontname','arial');
%             fillPage(gcf,'papersize',[8 5]);
%             pause
%             clf
%         end
%     end
% end

% %%
% close all
% f1 = figure();
% f2 = figure();
% for ss = 1:length(tr_set)
% % sbeg = find(all_mod_SU(tr_set) > 0,1);
% % for ss = sbeg:length(tr_set)
%     ss
%     init_mod = all_mod_fits(tr_set(ss));
%     xtargs = [init_mod.mods(:).Xtarget];
%     kmat = [init_mod.mods(xtargs == 1).filtK];
%     figure(f1); clf
%     subplot(2,2,1)
%     imagesc(reshape(kmat(:,1),flen,use_nPix_us));
%     ca = max(abs(kmat(:,1))); caxis([-ca ca]);
%     for ii = 1:(size(kmat,2)-1)
%         subplot(2,2,2+ii)
%         imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us));
%         ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
%     end
%     colormap(gray)
%
%     fin_mod = dit_mods{end}(tr_set(ss));
%     xtargs = [fin_mod.mods(:).Xtarget];
%     kmat = [fin_mod.mods(xtargs == 1).filtK];
%     figure(f2); clf
%     subplot(2,2,1)
%     imagesc(reshape(kmat(:,1),flen,use_nPix_us));
%     ca = max(abs(kmat(:,1))); caxis([-ca ca]);
%     for ii = 1:(size(kmat,2)-1)
%         subplot(2,2,2+ii)
%         imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us));
%         ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
%     end
%     colormap(gray)
%
%     fprintf('Cell %d of %d\n',ss,length(tr_set));
%     fprintf('Original: %.4f  Fin: %.4f\n',all_mod_LLimp(tr_set(ss)),dit_LLimp(end,tr_set(ss)));
%     pause
% end
%
