clear all

addpath('~/James_scripts/bruce/eye_tracking/');
addpath('~/James_scripts/bruce/processing/');

global Expt_name bar_ori use_LOOXV monk_name rec_type rec_number

Expt_name = 'M297';
monk_name = 'lem';
bar_ori = 0; %bar orientation to use (only for UA recs)
rec_number = 1;

use_LOOXV = 1; %[0 no LOOXV; 1 SU LOOXV; 2 all LOOXV]

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

fused = find(cellfun(@(x) length(x),Expts) > 0,1,'first');
%is this a laminar probe or utah array rec?
if strcmp(Expts{fused}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{fused}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
fprintf('Loading %s\n',data_name);
load(data_name);

mod_data_name = 'full_eyetrack_initmods_FIN2';
anal_name = 'full_eyetrack_FIN2';

base_sname = 'rpt_xval_test';

%%
recompute_init_mods = 1; %use existing initial models?
use_measured_pos = 3; %1 for init with coils, 2 for init with trial-sub coils, 3 for random init,
use_sac_kerns = 1; %use sac-modulation kernels
model_pop_avg = 1; %include an explicit model of population avg (and use it as a predictor in unit models)
pop_avg_sigma = 0.05; %gaussian smoothing sigma for pop avg rates

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
        case 13
            cor_ori = 100;
        case 14
            cor_ori = 40;
        case 320
            cor_ori = 100;
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
if use_measured_pos == 3
    mod_data_name = [mod_data_name '_Rinit'];
    anal_name = [anal_name '_Rinit'];
end
%if using coil info
if any(params.use_coils > 0)
    anal_name = [anal_name '_Cprior'];
end

mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];
anal_name = [anal_name sprintf('_ori%d',bar_ori)];
sname = [base_sname sprintf('_ori%d',bar_ori)];
if rec_number > 1
    mod_data_name = strcat(mod_data_name,sprintf('r%d',rec_number));
    anal_name = strcat(anal_name,sprintf('r%d',rec_number));
    sname = strcat(sname,sprintf('r%d',rec_number));
end
%%
xv_frac = 0.2; %fraction of trials to use for cross-validation
if ~isempty(params.rpt_seeds)
    xv_type = 'rpt';
else
    xv_type = 'uni';
end

flen = 12; %number of stimulus time lags to use

%number of spatial pixels to use for models
if strcmp(rec_type,'LP')
    use_nPix = 32;
elseif strcmp(rec_type,'UA')
    use_nPix = 16;
end
%for manually specifying number of used pixs
if ismember(Expt_num,[287 289 294])
    use_nPix = 15;
elseif ismember(Expt_num,[296 297 9 10 11 13 320])
    use_nPix = 22;
elseif ismember(Expt_num,[5 309 14])
    use_nPix = 26;
elseif Expt_num == 12
    if rec_number == 1
        use_nPix = 22;
    elseif rec_number == 2
        use_nPix = 28;
    end
end

%for bar widths bigger than this do additional spatial upsampling
if mode(expt_data.expt_dw)/params.scale_fac > 1.5*.0565
    add_usfac = 2;
else
    add_usfac = 1;
end
base_spatial_usfac = 2; %baseline spatial up-sampling
spatial_usfac = add_usfac*base_spatial_usfac;

n_fix_inf_it = 3; %3 number of iterations to estimate fixation-based EP
n_drift_inf_it = 1; %1 number of iterations to estimate within-fixation drift

fix_prior_sigma = 0.15; %prior sigma on fixation-based EP
fix_noise_sigma = 0.1; %sigma on noise of eye-trackers

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
drift_jump_sigma = 0.075; %prior sigma on change
drift_dsf = 3; %temporal down-sampling for estimating drift

min_trial_dur = 0.75; %minimal trial dur (s)

stim_fs = 100; %in Hz
dt = 0.01; %in sec
Fr = 1;

full_nPix=36; %total number of stimulus pixels to keep track of
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

if full_nPix ~= params.full_nPix
    fprintf('Using full_nPix in params struct\n');
    full_nPix = params.full_nPix;
end
if use_nPix > full_nPix
    fprintf('Using npix == full_nPix\n');
    use_nPix = full_nPix;
end

%exclude data at beginning and end of each trial (numbers in sec)
beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

%saccade kernel time axis
sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

stim_params = NIMcreate_stim_params([flen full_nPix],dt);

%%
base_sp_dx = mode(expt_data.expt_dw); %size of bars in pixels
if length(unique(expt_data.expt_dw)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac/params.scale_fac; %model dx in deg
init_sp_dx = base_sp_dx/base_spatial_usfac/params.scale_fac; %dx for initial modeling

max_shift = round(0.8475/sp_dx); %max eye pos deviation (deg)
dshift = 1;
shifts = -max_shift:dshift:max_shift; %shift-range
n_shifts = length(shifts); %number of shifts
zero_frame = find(shifts==0);

%shifts for inferring drift
max_Dshift = round(0.452/sp_dx);
Dshifts = -max_Dshift:dshift:max_Dshift;
n_Dshifts = length(Dshifts);
zero_Dframe = find(Dshifts==0);

%total shifts
max_Tshift = max_shift + max_Dshift;
Tshifts = -max_Tshift:dshift:max_Tshift;

%%
all_stim_mat = decompressTernNoise(stimComp);

%% get stim-matrix for initial model estimation
full_nPix_us = spatial_usfac*full_nPix;
use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;
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

%% get trial and block data
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

%indices of saccade start and end times
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));

saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds)); %indices of saccade starts relative to used_inds
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds)); %set of saccades occuring during used_inds
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);
used_is_blink = ET_data.is_blink(used_saccade_set); %which of used saccades are actually blinks

sac_amps = [ET_data.saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1; %identify potential microsacs
big_sacs = find(~is_micro & ~used_is_blink'); %saccades
micro_sacs = find(is_micro & ~used_is_blink'); %microsaccades
sac_durs = [ET_data.saccades(:).duration];

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

%create saccade/micro Xmatrices
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
%index values of trial starts and stops
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

%index values of fixation starts and stops
fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds; %duration of each fixations
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink));
n_fixs = length(fix_start_inds);

%push the effects of saccades forward in time for inferring eye positions
sac_shift = round(0.05/dt); %how much to shift jump points forward in time
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
%compute fixation ids both raw and forward-projected
fix_ids = nan(NT,1);
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
    cur_inds = pfix_start_inds(ii):(pfix_stop_inds(ii));
    pfix_ids(cur_inds) = ii;
end

all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);

%% Create set of TR and XV trials
if strcmp(xv_type,'rpt')
    rpt_trials = find(ismember([trial_data(:).se],params.rpt_seeds));
else
    rpt_trials = [];
end

use_trials = unique(all_trialvec(used_inds));
use_trials(ismember(use_trials,rpt_trials)) = []; %don't use rpt trials for model fitting
nuse_trials = length(use_trials);

%set aside a fraction of trials for xval
n_xv_trials = round(xv_frac*nuse_trials);
xv_trials = randperm(nuse_trials);
xv_trials(n_xv_trials+1:end) = [];
xv_trials = use_trials(xv_trials);
tr_trials = setdiff(use_trials,xv_trials);
n_tr_trials = length(tr_trials);
fprintf('Initializing models with %d training trials and %d xval trials\n',n_tr_trials,n_xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

%%
% make Robs_mat
tot_sus = size(all_binned_sua,2);
tot_nUnits = length(su_probes) + n_probes;
Robs_mat = nan(length(used_inds),n_probes + tot_sus);
for ss = 1:size(Robs_mat,2)
    if ss > n_probes
        Robs_mat(:,ss) = all_binned_sua(used_inds,ss-n_probes);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end

if model_pop_avg
    norm_rates = Robs_mat(:,1:n_probes); %define pop rate over MUA only
    for tt = 1:length(trial_data) %loop over trials
        cur_trial_inds = find(all_trialvec(used_inds) == tt);
        if ~isempty(cur_trial_inds)
            for ss = 1:n_probes %smooth the spk cnt data with a Gaussian kernel
                norm_rates(cur_trial_inds,ss) = jmm_smooth_1d_cor(norm_rates(cur_trial_inds,ss),round(pop_avg_sigma/dt));
            end
        end
    end
    
    %zscore each MU
    norm_rates = bsxfun(@minus,norm_rates,nanmean(norm_rates));
    norm_rates = bsxfun(@rdivide,norm_rates,nanstd(norm_rates));
    pop_rate = nanmean(norm_rates,2); %then avg
    
    n_block_filts = n_blocks + 1;
else
    n_block_filts = n_blocks;
end

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
cd(anal_dir);
mod_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us/add_usfac],dt);
mod_stim_params(2) = NMMcreate_stim_params([n_block_filts 1],dt);
if use_sac_kerns
    mod_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
    mod_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);
end

null_stim_params = mod_stim_params(2:end);

silent = 1;
sac_d2t = 100; %temporal smoothness regularization on saccade kernels
block_L2 = 1; %small L2 penalty on block filter kernels

%baseline regularization strengths
base_lambda_d2XT = 100;
base_lambda_L1 = 5;
init_lambda_d2XT = 100;
%if using sparser stimuli, reduce base reg strenghts to approximately
%account for difference in stimulus marginal variance
if unique(expt_data.expt_dds) == 12
    base_lambda_d2XT = base_lambda_d2XT/5;
    base_lambda_L1 = base_lambda_L1/2;
    init_lambda_d2XT = init_lambda_d2XT/5;
end

%initial optimization tolerances
init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-9;

sac_reg_params = NMMcreate_reg_params('lambda_d2T',sac_d2t); %saccade kernel reg params
%make reg params for the null model
if use_sac_kerns
    null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);
else
    null_reg_params = NMMcreate_reg_params('lambda_L2',block_L2);
end

n_squared_filts = 2; %number of squared filters
mod_signs = ones(1,n_squared_filts+2); %total number of subunits
%if using more than 2 quadratic filters, take the additional ones as
%suppressive
if n_squared_filts > 2
    mod_signs(n_squared_filts+1) = -1;
end
NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}]; %NL types for full model
Xtargets = [ones(n_squared_filts+1,1); 2]; %xtargets for full model

%initial regularization parameters
init_d2XT = [init_lambda_d2XT*ones(n_squared_filts+1,1); 0;];
init_L2 = [zeros(n_squared_filts+1,1); block_L2];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2);

%range of smoothness regularization 'rescaling' to try
poss_lambda_scales = logspace(-2,2,10);
% poss_lambda_scales = 1;

%create X matrix
if add_usfac > 1
    %if doing additional spatial up-sampling use tent basis functions
    X{1} = tb_proc_stim(all_Xmat_us(used_inds,use_kInds_up),add_usfac,flen);
else
    X{1} = all_Xmat_us(used_inds,use_kInds_up);
end
if model_pop_avg
    X{2} = [Xblock(used_inds,:) pop_rate];
else
    X{2} = [Xblock(used_inds,:)];
end
if use_sac_kerns
    X{3} = Xsac;
    X{4} = Xmsac;
end

%%
fprintf('Loading pre-computed initial models\n');
load(mod_data_name);

new_xv_inds = find(ismember(all_trialvec(used_inds),rpt_trials));

clear xvLL LL null_* rpt_*LL
for ss = 1:tot_nUnits
    fprintf('Computing base LLs for Unit %d of %d\n',ss,tot_nUnits);
    cur_tr_inds = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
    cur_xv_inds = xv_inds(~isnan(Robs_mat(xv_inds,ss)));
    all_modfit_inds = union(cur_tr_inds,cur_xv_inds);
    tr_NT = length(cur_tr_inds);
    xv_NT = length(cur_xv_inds);
    
    if ~isempty(cur_tr_inds)
        Robs = Robs_mat(:,ss);
        cur_mod = all_mod_fits(ss);
        cur_nullmod = all_nullmod(ss);
        if ~model_pop_avg
           cur_mod.mods(4).filtK(end) = [];
           cur_mod.stim_params(2).stim_dims(1) = cur_mod.stim_params(2).stim_dims(1)-1;
           cur_nullmod.mods(1).filtK(end) = [];
           cur_nullmod.stim_params(2).stim_dims(1) = cur_nullmod.stim_params(2).stim_dims(1)-1;
        end
        cur_new_xvinds = new_xv_inds(~isnan(Robs(new_xv_inds)));
        
        
        new_mod_fits(ss) = NMMfit_filters(cur_mod,Robs,X,[],cur_tr_inds);
%         new_mod_fits(ss) = NMMfit_logexp_spkNL(new_mod_fits(ss),Robs,X,[],cur_tr_inds);
        null_mod_fits(ss) = NMMfit_filters(cur_nullmod,Robs,X,[],cur_tr_inds);
        [xvLL(ss),null_xvLL(ss)] = NMMeval_model(new_mod_fits(ss),Robs,X,[],cur_xv_inds);
        [LL(ss),null_LL(ss)] = NMMeval_model(new_mod_fits(ss),Robs,X,[],cur_tr_inds);
        [rpt_xvLL(ss),rpt_null_xvLL(ss)] = NMMeval_model(new_mod_fits(ss),Robs,X,[],cur_new_xvinds);
    end
end

%% DIAGNOSTICS
full_prates = nan(NT,tot_nUnits);
full_nullrates = nan(NT,tot_nUnits);
for cc = 1:tot_nUnits
    if ~isempty(new_mod_fits(cc).mods)
        [~,~,full_prates(:,cc)] = NMMeval_model(new_mod_fits(cc),Robs_mat(:,cc),X);
        [~,~,full_nullrates(:,cc)] = NMMeval_model(null_mod_fits(cc),Robs_mat(:,cc),X);
    end
end
full_modLL = Robs_mat.*log(full_prates) - full_prates;
full_nullLL = Robs_mat.*log(full_nullrates) - full_nullrates;
full_LLimp = full_modLL-full_nullLL;

use_trial_set = unique(all_trialvec(used_inds));
n_use_trials = length(use_trial_set);
trial_blocknums = nan(n_use_trials,1);
trial_LLimp = nan(n_use_trials,tot_nUnits);
trial_avg_rates = nan(n_use_trials,tot_nUnits);
for tt = 1:n_use_trials
    cur_used_inds = find(all_trialvec(used_inds) == use_trial_set(tt));
    trial_LLimp(tt,:) = nanmean(full_LLimp(cur_used_inds,:));
    trial_avg_rates(tt,:) = nanmean(Robs_mat(cur_used_inds,:));
end
trial_LLimp_zscore = nanzscore(trial_LLimp);

%%
save(sname,'new_mod_fits','xvLL','null_xvLL','LL','null_LL','rpt_xvLL','rpt_null_xvLL');
