% clear all

global Expt_name bar_ori monk_name rec_type

% Expt_name = 'M266';
% monk_name = 'lem';
% bar_ori = 80; %bar orientation to use (only for UA recs)
rec_number = 1;

Expt_num = str2num(Expt_name(2:end));

data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
load(Edata_file);

%is this a laminar probe or utah array rec?
if strcmp(Expts{1}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{1}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
fprintf('Loading %s\n',data_name);
load(data_name);

%%
recompute_init_mods = 0; %use existing initial models?
use_measured_pos = 3; %1 for init with coils, 2 for init with trial-sub coils, 3 for random init, 0 for init to perfect fixation
use_sac_kerns = 1; %use sac-modulation kernels

%%
if strcmp(rec_type,'LP')
    switch Expt_num
        case 266
            cor_ori = 80;
        case 270
            cor_ori = 60;
        case 275
            cor_ori = 135;
        case 277
            cor_ori = 70;
        case 281
            cor_ori = 140;
        case 287
            cor_ori = 90;
        case 289
            cor_ori = 160;
        case 294
            cor_ori = 40;
        case 296
            cor_ori = 45;
        case 297
            cor_ori = [0 90];
        case 309
            cor_ori = 120;
        case 5
            cor_ori = 50;
        case 9
            cor_ori = 0;
        case 10
            cor_ori = 60;
        case 11
            cor_ori = 160;
        case 12
            cor_ori = 0;
    end
else
    cor_ori = [0 90];
end
if ~ismember(bar_ori,cor_ori)
    error('this bar ori is not allowed');
end

cd(data_dir);

if strcmp(rec_type,'LP')
    n_probes = 24;
elseif strcmp(rec_type,'UA')
    n_probes = 96;
end

load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

if rec_number > 1
    cluster_dir = [cluster_dir sprintf('/rec%d',rec_number)];
end


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
if any(params.use_coils > 0)
    old_anal_name = [old_anal_name '_Cprior'];
end

hr_anal_name = [old_anal_name '_hres'];
hr_mod_name = [mod_data_name 'hres'];

mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];
old_anal_name = [old_anal_name sprintf('_ori%d',bar_ori)];
hr_anal_name = [hr_anal_name sprintf('_ori%d',bar_ori)];
hr_mod_name = [hr_mod_name sprintf('_ori%d',bar_ori)];

load([anal_dir '/' old_anal_name],'et_params');
%%

use_nPix = et_params.use_nPix;
flen = et_params.flen;

if ~isnan(params.rpt_seeds)
    xv_type = 'rpt';
    xv_frac = nan;
else
    xv_type = 'uni';
    xv_frac = 0.2;
end
old_spatial_usfac = et_params.spatial_usfac;
spatial_usfac = old_spatial_usfac*2;

n_drift_inf_it = 1; %3

fix_prior_sigma = 0.15;
fix_noise_sigma = 0.1;

%this sets the posterior sigma to be ~0.004 in all cases
drift_noise_sigma = 0;
if all(params.use_coils)
    drift_prior_sigma = sqrt(0.004^2*3);
    drift_noise_sigma = sqrt(0.004^2*3);
elseif any(params.use_coils)
    drift_prior_sigma = sqrt(0.004^2*2);
    drift_noise_sigma = sqrt(0.004^2*2);
else
    drift_prior_sigma = 0.004;
end
drift_jump_sigma = 0.075; %0.05 start
drift_dsf = 2;

stim_fs = 100; %in Hz
dt = 0.01;
Fr = 1;

full_nPix = params.full_nPix;

%exclude data at beginning and end of each trial
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
base_sp_dx = mode(expt_data.expt_dw);
if length(unique(expt_data.expt_dw)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac/params.scale_fac; %model dx in deg

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

%%
cd(anal_dir)

fprintf('Loading pre-computed initial models\n');
load(mod_data_name);
load(old_anal_name,'dit_mods','et_tr_set','it_fix_post_mean','drift_post_mean','it_fix_post_std','drift_post_std','et_saccades');
old_best_mods = dit_mods{end};
tr_set = et_tr_set;

up_fac = spatial_usfac/old_spatial_usfac;
best_fix_cor = it_fix_post_mean(end,:)*up_fac;
best_fix_std = it_fix_post_std(end,:)*up_fac;
best_drift_cor = drift_post_mean(end,:)*up_fac;
best_drift_std = drift_post_std(end,:)*up_fac;

clear it_fix_post_mean drift_post_mean dit_mods

%%
%for some recs the saccades for ET_data were detected using a slightly
%different algo. Eliminate the saccades causing the difference
if length(ET_data.saccades) ~= length(et_saccades)
    sac_start_times = [ET_data.saccades(:).start_time];
    old_sac_start_times = [et_saccades(:).start_time];
    if length(sac_start_times) > length(old_sac_start_times)
        extra_sacs = find(~ismember(sac_start_times,old_sac_start_times));
        ET_data.saccades(extra_sacs) = [];
        ET_data.is_blink(extra_sacs) = [];
        fprintf('Difference in saccade detection from eye-tracking data, eliminating %d/%d saccades\n',length(extra_sacs),length(sac_start_times));
    else
        error('Fewer saccades than in ETdata');
    end
end

%%
all_stim_mat = decompressTernNoise(stimComp);

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
all_binned_mua = spikes_int82double(spike_data.binned_mua);
all_binned_sua = spikes_int82double(spike_data.binned_sua);
Clust_data = spike_data.Clust_data;
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

%%
NT = length(used_inds);
fullNT = size(all_binned_mua,1);
n_trials = length(time_data.trial_flip_ids);
n_blocks = length(expt_data.used_blocks);

all_t_axis = time_data.t_axis;
trial_start_inds = [1; 1+time_data.trial_flip_inds(2:end)];
trial_end_inds = [time_data.trial_flip_inds(2:end); fullNT];
all_trialvec = nan(fullNT,1);
for ii = 1:n_trials
    all_trialvec(trial_start_inds(ii):trial_end_inds(ii)) = time_data.trial_flip_ids(ii);
end

block_start_inds = [1; 1+time_data.block_flip_inds(2:end)];
block_end_inds = [time_data.block_flip_inds(2:end); fullNT];
all_blockvec = nan(fullNT,1);
for ii = 1:n_blocks
    all_blockvec(block_start_inds(ii):block_end_inds(ii)) = time_data.block_flip_ids(ii);
end
Xblock = zeros(fullNT,n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% CREATE SACCADE PREDICTOR MATS
corrected_eye_vals_interp = ET_data.interp_eye_pos;
sac_start_times = [ET_data.saccades(:).start_time];
sac_stop_times = [ET_data.saccades(:).stop_time];

interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));

saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);
used_is_blink = ET_data.is_blink(used_saccade_set);

sac_amps = [ET_data.saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro & ~used_is_blink');
micro_sacs = find(is_micro & ~used_is_blink');
sac_durs = [ET_data.saccades(:).duration];

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
    xv_trials = find(ismember([trial_data(:).se],params.rpt_seeds));
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
if strcmp(xv_type,'uni')
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
silent = 1;
sac_d2t = 100;

if unique(expt_data.expt_dds) == 67
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

fprintf('Loading pre-computed initial models\n');
load(hr_mod_name);

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
if params.use_coils(1) == 1 && params.use_coils(2) == 1
    
    measured_fix_deltas = nan(n_fixs,1);
    measured_fix_deltas(2:end) = mean(diff(measured_fix_avg),2);
    fix_noise_sigma = fix_noise_sigma/sqrt(2); %noise variance decreases by factor of sqrt(2) with 2 independent measures
    post_drift_var = 1/(2/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    post_mean_drift = post_drift_var*2*mean(measured_drift,2)/drift_noise_sigma^2;
    
    %if using left coil only
elseif params.use_coils(1) == 1
    measured_fix_deltas = nan(n_fixs,1);
    measured_fix_deltas(2:end) = diff(measured_fix_avg(:,1));
    post_drift_var = 1/(1/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    post_mean_drift = post_drift_var*measured_drift(:,1)/drift_noise_sigma^2;
    
    %if using right coil only
elseif params.use_coils(2) == 1
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

%%
it_mods = all_mod_fits;
it_mods_spkNL = all_mod_fits_withspkNL;
it_LLimp = all_mod_LLimp;
it_R2 = all_mod_R2;

%% PREPROCESS MODEL COMPONENTS
loo_set = (n_probes+1):n_tr_chs; %leave out all SUs

cur_uset = setdiff(tr_set,loo_set);
n_uset = length(cur_uset);
filt_bank = zeros(n_uset,klen_us,n_squared_filts+1);
lin_kerns = nan(n_uset,n_blocks);
if use_sac_kerns
    sac_kerns = nan(n_uset,n_sac_bins);
    msac_kerns = nan(n_uset,n_sac_bins);
end
mod_spkNL_params = nan(n_uset,3);
for ss = 1:n_uset
    cur_Xtargs = [it_mods(tr_set(cur_uset(ss))).mods(:).Xtarget];
    cur_k = [it_mods(tr_set(cur_uset(ss))).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = it_mods_spkNL(tr_set(cur_uset(ss))).spk_NL_params(1:3);
    lin_kerns(ss,:) = it_mods(tr_set(cur_uset(ss))).mods(cur_Xtargs == 2).filtK;
    if use_sac_kerns
        sac_kerns(ss,:) = it_mods(tr_set(cur_uset(ss))).mods(cur_Xtargs == 3).filtK;
        msac_kerns(ss,:) = it_mods(tr_set(cur_uset(ss))).mods(cur_Xtargs == 4).filtK;
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
cur_Xmat = all_Xmat_us(used_inds,:);
%precompute LL at all shifts for all units
frame_LLs = nan(NT,n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_stim_shift = cur_Xmat*shift_mat{xx};
    
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
fix_LLs = nan(n_fixs,n_shifts);
for ii = 1:n_fixs
    cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
    fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
end

if all(params.use_coils==0) %if not using coil info
    lgamma = bsxfun(@plus,fix_LLs,fix_prior);
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
else
    lalpha=zeros(n_fixs,n_shifts);
    lbeta = zeros(n_fixs,n_shifts);
    %compute forward messages
    lalpha(1,:) = fix_prior + fix_LLs(1,:);
    for t=2:n_fixs
        if ~fix_post_blink(t) %if not following a blink
            cdist = pdist2(shifts'*sp_dx + measured_fix_deltas(t),shifts'*sp_dx);
            cur_lA = -cdist.^2/(2*fix_noise_sigma^2);
            cur_lA = bsxfun(@plus,cur_lA,fix_prior);
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + fix_LLs(t,:);
        else
            lalpha(t,:) = fix_prior + fix_LLs(t,:);
        end
    end
    
    %compute backward messages
    lbeta(n_fixs,:)=zeros(1,n_shifts);
    for t=n_fixs-1:-1:1
        if ~fix_post_blink(t+1)
            cdist = pdist2(shifts'*sp_dx + measured_fix_deltas(t+1),shifts'*sp_dx);
            cur_lA = -cdist.^2/(2*fix_noise_sigma^2);
            cur_lA = bsxfun(@plus,cur_lA,fix_prior);
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
            
            lf1 = lbeta(t+1,:) + fix_LLs(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA');
        else
            lbeta(t,:) = zeros(1,n_shifts);
        end
    end
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
end

gamma = exp(lgamma);

it_fix_post_mean = sum(bsxfun(@times,gamma,shifts),2);
cur_diff = bsxfun(@minus,it_fix_post_mean,shifts).^2;
it_fix_post_std = sqrt(sum(cur_diff.*gamma,2));

%back-project saccade-times
all_fix_post_mean_cor = nan(NT,1);
all_fix_post_mean_cor(~isnan(fix_ids)) = it_fix_post_mean(fix_ids(~isnan(fix_ids)));
all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;

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
        
        if all(params.use_coils==0)
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
drift_post_std = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2)) - drift_post_mean.^2;

drift_post_mean(isnan(drift_post_mean)) = 0;
drift_post_std(isnan(drift_post_std)) = 0;
drift_post_mean = interp1(find(~isnan(pfix_ids)),squeeze(drift_post_mean(~isnan(pfix_ids))),1:NT);
drift_post_std = interp1(find(~isnan(pfix_ids)),squeeze(drift_post_std(~isnan(pfix_ids))),1:NT);

%% SAVE EYE-TRACKING RESULTS
et_params = struct('beg_buffer',params.beg_buffer,'end_buffer',params.end_buffer,'min_trial_dur',params.min_trial_dur,'bar_ori',bar_ori,'good_coils',params.good_coils,...
    'use_nPix',use_nPix,'flen',flen,'dt',dt,'drift_jump_sigma',drift_jump_sigma,'drift_prior_sigma',drift_prior_sigma,...
    'fix_prior_sigma',fix_prior_sigma,'fix_noise_sigma',fix_noise_sigma,'drift_noise_sigma',drift_noise_sigma,...
    'drift_dsf',drift_dsf,'use_sac_kerns',use_sac_kerns,'shifts',shifts,...
    'use_measured_pos',use_measured_pos,'sac_bincents',sac_bincents,'spatial_usfac',spatial_usfac,'old_spatial_usfac',old_spatial_usfac,'sac_shift',sac_shift,'use_coils',params.use_coils,'sp_dx',sp_dx);

et_used_inds = used_inds;
et_tr_set = tr_set;
et_saccades = ET_data.saccades;
et_is_blink = ET_data.is_blink;
et_clust_data = Clust_data;
cd(anal_dir);

hr_anal_name = strcat(hr_anal_name,'_fullLOO');

save(hr_anal_name,'it_fix*','drift_post_*','fix_ids','et_used_inds','et_tr_set','et_clust_data','et_saccades','et_is_blink','et_params');

%%
% fin_fix_corr = nan(NT,1);
% fin_fix_std = nan(NT,1);
% fin_fix_corr(~isnan(fix_ids)) = best_fix_cor(end,fix_ids(~isnan(fix_ids)));
% fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
% fin_fix_std(~isnan(fix_ids)) = best_fix_std(end,fix_ids(~isnan(fix_ids)));
% fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
% 
% fin_fix_corr = fin_fix_corr*sp_dx;
% fin_fix_std = fin_fix_std*sp_dx;
% 
% fin_drift_corr = drift_post_mean(end,:)*sp_dx;
% fin_drift_std = drift_post_std(end,:)*sp_dx;
% 
% for ii = 1:length(trial_start_inds)
%     cur_inds = trial_start_inds(ii):trial_end_inds(ii);
%     fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
%     fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
% end
% fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
% fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);
% 
% 
% fin_tot_corr = fin_fix_corr + fin_drift_corr;
% fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);

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
%             %             h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_fix_corr(uu),fin_fix_std(uu),{'color','m'});
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

