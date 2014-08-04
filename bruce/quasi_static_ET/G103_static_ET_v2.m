clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

%run on 86 [0, 90];
Expt_num = 103;
Expt_name = sprintf('G%.3d',Expt_num);

use_LOOXV = 0; %[0 no LOOXV; 1 SU LOOXV; 2 all LOOXV]
bar_ori = 0; %bar orientation to use (only for UA recs)

recompute_init_mods = 0; %use existing initial models?
use_measured_pos = 0; %1 for init with coils, 2 for init with trial-sub coils, 0 for init to perfect fixation
use_sac_kerns = 1; %use sac-modulation kernels


data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

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
    end
end

if Expt_num >= 280
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
anal_name = 'full_eyetrack';

if bar_ori == 90 & strcmp(rec_type,'UA')
    mod_data_name = [mod_data_name '_vbars'];
    anal_name = [anal_name '_vbars'];
end

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

%if using coil info
if any(use_coils > 0)
    anal_name = [anal_name '_Cprior'];
end
%if using coil initialization
if use_measured_pos == 1
    mod_data_name = [mod_data_name '_Cinit'];
    anal_name = [anal_name '_Cinit'];
end
%if using coil init with trial-sub
if use_measured_pos == 2
    mod_data_name = [mod_data_name '_CPinit'];
    anal_name = [anal_name '_CPinit'];
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
spatial_usfac = 2;

%these recs have larger bar widths
if ismember(Expt_num,[287 289 294])
    use_nPix = 15;
    spatial_usfac = 4;
end

n_fix_inf_it = 3; %3
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
drift_dsf = 3;

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

%%
include_expts = {'rls.sl'};
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
    included_type(ii) = any(strcmp(expt_names{ii},include_expts));
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);


cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

n_blocks = length(cur_block_set);

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
% all_Xmat = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_sl = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:n_blocks;
    fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,n_blocks);
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
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    trial_sl = [Expts{cur_block}.Trials(:).sl];
    trial_sl = trial_sl(id_inds);
    all_trial_sl = cat(1,all_trial_sl,trial_sl(use_trials)');
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
    cur_use_pix = (1:full_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
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
            %             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            %             all_Xmat = [all_Xmat; bar_Xmat];
            all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
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
all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);

%% select submatrix with central pixels

buffer_pix = floor((full_nPix - use_nPix)/2);

%repeat for up-sampled versions of the Xmatrix
[Xinds_up,Tinds_up] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
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
used_trials = find(all_trial_sl == 50);
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer & ...
    ismember(all_trialvec,used_trials));
NT = length(used_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_stim_times),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% PROCESS EYE TRACKING DATA
% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,good_coils);
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
%%
fprintf('Loading pre-computed initial models\n');
cd(anal_dir)
load(mod_data_name);

%if we used eye-trackers for initial model fitting, reconstruct uncorrected
%stimulus
if use_measured_pos > 0
    all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);
end

load(anal_name);
%% SELECT USABLE UNITS AND make Robs_mat
% usable_units = find(all_mod_xvLLimp > 0);
usable_units = find(~isnan(all_mod_xvLLimp));

n_used_sus = sum(all_mod_SU(usable_units) ~= 0);
n_used_mus = sum(all_mod_SU(usable_units) == 0);
fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
tr_set = usable_units;
n_tr_chs = length(tr_set);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    if all_mod_SU(tr_set(ss)) > 0
        su_probe_ind = find(Clust_data.SU_numbers == all_mod_SUnum(tr_set(ss)));
        Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,tr_set(ss));
    end
end

%don't use separate xv set for eye-tracking
if strcmp(rpt_trials,'uni')
    tr_inds = sort([tr_inds; xv_inds]);
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
%repeat for up-sampled versions of the Xmatrix
buffer_pix = floor((full_nPix - use_nPix)/2);
[Xinds_up,Tinds_up] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1);
cur_use_pix = (1/spatial_usfac:1/spatial_usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

all_static_stim = all_stimmat_up(used_inds,use_kInds_up);

%%
dstim = diff(all_stimmat_up);
stim_flips = find(any(dstim ~= 0,2))+1;
stim_flips = [1; stim_flips];
stim_flips = find(ismember(used_inds,stim_flips));
stim_inds = [stim_flips(1:end-1) stim_flips(2:end)-1];

fix_start_inds = sort([stim_inds(:,1); saccade_stop_inds]);
fix_stop_inds = sort([stim_inds(:,2); saccade_start_inds]);

fix_durs = fix_stop_inds-fix_start_inds;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink));
n_fixs = length(fix_start_inds);

t_since_fix_start = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    t_since_fix_start(cur_inds) = 0:(length(cur_inds)-1);
end

nposs_lags = 50;
XT = zeros(NT,nposs_lags);
for ii = 1:nposs_lags
    curset = find(t_since_fix_start == ii);
    XT(curset,ii) = 1;
end


%%
filt_bank = zeros(n_tr_chs,use_nPix_us,3);
for cc = 1:size(Robs_mat,2);
    
    fprintf('Fitting static mod %d/%d\n',cc,size(Robs_mat,2));
    
    cur_mod = dit_mods{end}(cc);
    xtargs = [cur_mod.mods(:).Xtarget];
    cur_filts = [cur_mod.mods(xtargs==1).filtK];
    re_filts = reshape(cur_filts,[flen use_nPix_us 3]);
    tempkern = squeeze(mean(std(re_filts,[],2),3));
    [~,bestlag] = max(tempkern);
    
    filt_slices = squeeze(re_filts(bestlag,:,:));
    filt_bank(cc,:,:) = filt_slices;
    cur_Robs = Robs_mat(:,cc);
    filt_outs = all_static_stim*filt_slices;
    filt_outs(:,2:end) = filt_outs(:,2:end).^2;
    
    gout = sum(filt_outs,2);
    
    X{1} = XT;
    X{2} = bsxfun(@times,XT,gout);
    X{3} = Xblock(used_inds,:);
    
    nim_stim_params(1:2) = NMMcreate_stim_params(50);
    nim_stim_params(3) = NMMcreate_stim_params(n_blocks);
    reg_params = NMMcreate_reg_params('lambda_d2T',100);
    
    nim1 = NMMinitialize_model(nim_stim_params,[1 1 1],{'lin','lin','lin'},reg_params,[1 2 3]);
    nim1 = NMMfit_filters(nim1,cur_Robs,X,[],[],1);
    
    all_static_mod(cc) = nim1;
    
    nullMod = NMMinitialize_model(nim_stim_params,[1 1],{'lin','lin'},reg_params,[1 3]);
    nullMod = NMMfit_filters(nullMod,cur_Robs,X,[],[],1);
    
    all_null_mod(cc) = nullMod;
    
    all_static_LLimp(cc) = all_static_mod(cc).LL_seq(end) - all_null_mod(cc).LL_seq(end);
end

filt_bank = permute(filt_bank,[2 1 3]);

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

fix_prior = -(shifts*sp_dx).^2./(2*fix_prior_sigma^2);
fix_prior = fix_prior - logsumexp(fix_prior);

%%
n_loo = 10;
loo_set = randperm(n_tr_chs);
loo_set = sort(loo_set(1:n_loo));

%%
it_mods{1} = all_static_mod;
it_LLimp(1,:) = all_static_LLimp;
clear it_fix_post*
it_LLimp_LOO(length(loo_set),1,:) = all_static_LLimp;
for xv = 1:length(loo_set)
    it_mods_LOO{xv,1} = all_static_mod;
end
for nn = 1:n_fix_inf_it
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,n_fix_inf_it);
    
    %% PREPROCESS MODEL COMPONENTS
    Tfilt_bank = zeros(n_tr_chs,nposs_lags,2);
    lin_kerns = nan(n_tr_chs,n_blocks);
    
    mod_spkNL_params = nan(n_tr_chs,3);
    for ss = 1:n_tr_chs
        cur_Xtargs = [it_mods{nn}(ss).mods(:).Xtarget];
        cur_k = [it_mods{nn}(ss).mods(cur_Xtargs == 1 | cur_Xtargs == 2).filtK];
        cur_bk = [it_mods{nn}(ss).mods(cur_Xtargs == 3).filtK];
        Tfilt_bank(ss,:,:) = cur_k;
        mod_spkNL_params(ss,:) = it_mods{nn}(ss).spk_NL_params(1:3);
        lin_kerns(ss,:) = cur_bk';
    end
    
    %indicator predictions
    block_out = Xblock(used_inds,:)*lin_kerns';
    toff_out = XT*squeeze(Tfilt_bank(:,:,1))';
    
    %% ESTIMATE LL for each shift in each stimulus frame
    cur_Xmat = all_stimmat_up(used_inds,:);
    %precompute LL at all shifts for all units
    frame_LLs = nan(NT,n_shifts);
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',xx,n_shifts);
        %     cur_stim_shift = cur_Xmat*shift_mat{xx};
        
        cur_stim_shift = shift_matrix_Nd(cur_Xmat,-shifts(xx),2);
        cur_stim_shift = cur_stim_shift(:,use_kInds_up);
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,n_tr_chs);
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
        
        gouts = cur_stim_shift*squeeze(filt_bank(:,:,1));
        for ff = 2:3
            gouts = gouts + (cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
        end
        
        cell_gouts = nan(NT,n_tr_chs);
        for cc = 1:n_tr_chs
            curX = bsxfun(@times,XT,gouts(:,cc));
            cell_gouts(:,cc) = curX*squeeze(Tfilt_bank(cc,:,2))';
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + cell_gouts + block_out + toff_out;
        
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
%         frame_LLs(:,xx,:) = Robs_mat.*log(pred_rate) - pred_rate;
    end
    
    %% INFER MICRO-SAC SEQUENCE
    fix_LLs = nan(n_fixs,n_shifts);
    for ii = 1:n_fixs
        cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
    end
    
    lgamma = bsxfun(@plus,fix_LLs,fix_prior);
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    gamma = exp(lgamma);
    
    it_fix_post_mean(nn,:) = sum(bsxfun(@times,gamma,shifts),2);
    cur_diff = bsxfun(@minus,it_fix_post_mean(nn,:)',shifts).^2;
    it_fix_post_std(nn,:) = sqrt(sum(cur_diff.*gamma,2));
    
    %back-project saccade-times
    all_fix_post_mean_cor = nan(NT,1);
    all_fix_post_mean_cor(~isnan(fix_ids)) = it_fix_post_mean(nn,fix_ids(~isnan(fix_ids)));
    all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
    all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;

    %%
    all_shift_stimmat_up = all_stimmat_up(used_inds,:);
    for i=1:NT
        all_shift_stimmat_up(i,:) = shift_matrix_Nd(all_shift_stimmat_up(i,:),-round(all_fix_post_mean_cor(i)),2);
    end
    all_shift_stimmat_up = all_shift_stimmat_up(:,use_kInds_up);
    
    %%
    for cc = 1:size(Robs_mat,2);
        fprintf('Fitting static mod %d/%d\n',cc,size(Robs_mat,2));
        
        cur_mod = dit_mods{end}(cc);
        xtargs = [cur_mod.mods(:).Xtarget];
        cur_filts = [cur_mod.mods(xtargs==1).filtK];
        re_filts = reshape(cur_filts,[flen use_nPix_us 3]);
        tempkern = squeeze(mean(std(re_filts,[],2),3));
        [~,bestlag] = max(tempkern);
        
        filt_slices = squeeze(re_filts(bestlag,:,:));
        cur_Robs = Robs_mat(:,cc);
        filt_outs = all_shift_stimmat_up*filt_slices;
        filt_outs(:,2:end) = filt_outs(:,2:end).^2;
        gout = sum(filt_outs,2);
        X{2} = bsxfun(@times,XT,gout);
        
        nim1 = it_mods{nn}(cc);
        nim1 = NMMfit_filters(nim1,cur_Robs,X,[],[],1);
        
        it_mods{nn+1}(cc) = nim1;
        it_LLimp(nn+1,cc) = (nim1.LL_seq(end) - all_null_mod(cc).LL_seq(end));
        if nn > 1
            fprintf('Original: %.4f  Prev: %.4f  New: %.4f\n',all_static_LLimp(cc),it_LLimp(nn,cc),it_LLimp(nn+1,cc));
        else
            fprintf('Original: %.4f  New: %.4f\n',all_static_LLimp(cc),it_LLimp(nn+1,cc));
        end
    end
    
    for xv = 1:length(loo_set)
        fprintf('Fixation iter %d, LOO %d of %d\n',nn,xv,length(loo_set));
        cur_uset = setdiff(1:n_tr_chs,loo_set(xv));
        n_uset = length(cur_uset);
        Tfilt_bank = zeros(n_uset,nposs_lags,2);
        lin_kerns = nan(n_uset,n_blocks);
        mod_spkNL_params = nan(n_uset,3);
        for ss = 1:length(cur_uset);
            cur_Xtargs = [it_mods_LOO{xv,nn}(cur_uset(ss)).mods(:).Xtarget];
            cur_k = [it_mods_LOO{xv,nn}(cur_uset(ss)).mods(cur_Xtargs == 1| cur_Xtargs == 2).filtK];
            cur_bk = [it_mods_LOO{xv,nn}(cur_uset(ss)).mods(cur_Xtargs == 3).filtK];
            Tfilt_bank(ss,:,:) = cur_k;
            mod_spkNL_params(ss,:) = it_mods_LOO{xv,nn}(cur_uset(ss)).spk_NL_params(1:3);
            lin_kerns(ss,:) = cur_bk';
        end
        cur_filt_bank = filt_bank(:,cur_uset,:);
        
        %indicator predictions
        block_out = Xblock(used_inds,:)*lin_kerns';
        toff_out = XT*squeeze(Tfilt_bank(:,:,1))';
        
        %% ESTIMATE LL for each shift in each stimulus frame
        for xx = 1:length(shifts)
            fprintf('Shift %d of %d\n',xx,n_shifts);
            cur_stim_shift = shift_matrix_Nd(cur_Xmat,-shifts(xx),2);
            cur_stim_shift = cur_stim_shift(:,use_kInds_up);
            
            %outputs of stimulus models at current X-matrix shift
            gfuns = ones(NT,n_uset);
            gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
            
            gouts = cur_stim_shift*squeeze(cur_filt_bank(:,:,1));
            for ff = 2:3
                gouts = gouts + (cur_stim_shift*squeeze(cur_filt_bank(:,:,ff))).^2;
            end
            
            cell_gouts = nan(NT,n_uset);
            for cc = 1:n_uset
                curX = bsxfun(@times,XT,gouts(:,cc));
                cell_gouts(:,cc) = curX*squeeze(Tfilt_bank(cc,:,2))';
            end
            
            %add contributions from extra lin kernels
            gfuns = gfuns + cell_gouts + block_out + toff_out;
            
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
        fix_LLs = nan(n_fixs,n_shifts);
        for ii = 1:n_fixs
            cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
            fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
        end
        
        lgamma = bsxfun(@plus,fix_LLs,fix_prior);
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
                
        gamma = exp(lgamma);
        
        it_fix_post_mean_LOO(xv,nn,:) = sum(bsxfun(@times,gamma,shifts),2);
        cur_diff = bsxfun(@minus,squeeze(it_fix_post_mean_LOO(xv,nn,:)),shifts).^2;
        it_fix_post_std_LOO(xv,nn,:) = sqrt(sum(cur_diff.*gamma,2));
        
        %back-project saccade-times
        all_fix_post_mean_cor = nan(NT,1);
        all_fix_post_mean_cor(~isnan(fix_ids)) = squeeze(it_fix_post_mean_LOO(xv,nn,fix_ids(~isnan(fix_ids))));
        all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
        all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;
        
        %% RECOMPUTE XMAT
        all_shift_stimmat_up = all_stimmat_up(used_inds,:);
        for i=1:NT
            all_shift_stimmat_up(i,:) = shift_matrix_Nd(all_shift_stimmat_up(i,:),-round(all_fix_post_mean_cor(i)),2);
        end
        all_shift_stimmat_up = all_shift_stimmat_up(:,use_kInds_up);
        
        %% REFIT XV CELLS
%         cur_X{1} = all_Xmat_up_fixcor(:,use_kInds_up);
%         
%         silent = 1;
%         for ss = 1:length(tr_set)
%             cur_cell = tr_set(ss);
%             fprintf('LOO %d/%d, Refitting model %d of %d\n',xv,length(loo_set),ss,length(tr_set));
%             cur_unit_ind = find(tr_set == cur_cell);
%             cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
%             
%             tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
%             
%             it_mods_LOO{xv,nn+1}(cur_cell) = it_mods{nn+1}(cur_cell);
%             it_mods_LOO{xv,nn+1}(cur_cell) = NMMfit_filters(it_mods_LOO{xv,nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
%                 tr_X,[],[],silent); %fit stimulus filters
%             
%             %refit spk NL
%             it_mods_spkNL_LOO{xv,nn+1}(cur_cell) = NMMfit_logexp_spkNL(it_mods_LOO{xv,nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
%             
%             [newLL,~,new_prate] = NMMmodel_eval(it_mods_spkNL_LOO{xv,nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
%             [~,~,null_prate] = NMMmodel_eval(all_nullmod(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X(2:end));
%             [it_R2_LOO(xv,nn+1,cur_cell),it_dev_LOO(xv,nn+1,cur_cell)] = pseudo_r2(Robs_mat(cur_tr_uset,cur_unit_ind),new_prate,null_prate);
%             
%             it_LLimp_LOO(xv,nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
%         end
        
        for cc = 1:size(Robs_mat,2);
            fprintf('Fitting static mod %d/%d\n',cc,size(Robs_mat,2));
            
            cur_mod = dit_mods{end}(cc);
            xtargs = [cur_mod.mods(:).Xtarget];
            cur_filts = [cur_mod.mods(xtargs==1).filtK];
            re_filts = reshape(cur_filts,[flen use_nPix_us 3]);
            tempkern = squeeze(mean(std(re_filts,[],2),3));
            [~,bestlag] = max(tempkern);
            
            filt_slices = squeeze(re_filts(bestlag,:,:));
            cur_Robs = Robs_mat(:,cc);
            filt_outs = all_shift_stimmat_up*filt_slices;
            filt_outs(:,2:end) = filt_outs(:,2:end).^2;
            gout = sum(filt_outs,2);
            X{2} = bsxfun(@times,XT,gout);
            
            nim1 = it_mods{nn+1}(cc);
            nim1 = NMMfit_filters(nim1,cur_Robs,X,[],[],1);
            
            it_mods_LOO{xv,nn+1}(cc) = nim1;
            it_LLimp_LOO(xv,nn+1,cc) = (nim1.LL_seq(end) - all_null_mod(cc).LL_seq(end));
        end
        
    end
end

