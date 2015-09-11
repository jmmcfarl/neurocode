clear all
addpath('~/James_scripts/NIMclass/')

Expt_name = 'M012';
monk_name = 'jbe';
bar_ori = 0; %bar orientation to use (only for UA recs)
rec_number = 1;

use_LOOXV = 0; %[0 no LOOXV; 1 SU LOOXV; 2 all LOOXV]

Expt_num = str2num(Expt_name(2:end));

%get data directory
data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
load(Edata_file); %load Expts file

fused = find(cellfun(@(x) length(x),Expts) > 0,1,'first'); %first block with data
%is this a laminar probe or utah array rec?
if strcmp(Expts{fused}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{fused}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

%load packaged data
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
fprintf('Loading %s\n',data_name);
load(data_name);

%data name!

if strcmp(rec_type,'LP')
    n_probes = 24;
elseif strcmp(rec_type,'UA')
    n_probes = 96;
end

% load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

if rec_number > 1
    cluster_dir = [cluster_dir sprintf('/rec%d',rec_number)];
end

% mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];
% anal_name = [anal_name sprintf('_ori%d',bar_ori)];
% if rec_number > 1
%     mod_data_name = strcat(mod_data_name,sprintf('r%d',rec_number));
%     anal_name = strcat(anal_name,sprintf('r%d',rec_number));
% end

%%
recompute_init_mods = 0; %use existing initial models?
use_sac_kerns = 1; %use sac-modulation kernels
model_pop_avg = 1; %include an explicit model of population avg (and use it as a predictor in unit models)
pop_avg_sigma = 0.05; %gaussian smoothing sigma for pop avg rates
init_jitter_SD = 0.075; %SD for initial random EP distribution
max_fix_shift = 0.85; %maximum shift to consider when estimating fixation-corrections (in deg)
max_drift_shift = 0.45; %max shift for drift corrections (in deg)

xv_frac = 0.2; %fraction of trials to use for cross-validation
flen = 12; %number of stimulus time lags to use
if ~isempty(params.rpt_seeds)
    has_rpts = true;
else
    has_rpts = false;
end

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

%for bar widths bigger than this do additional spatial upsampling (using
%interpolation, no filter up-sampling)
if mode(expt_data.expt_dw)/params.scale_fac > 1.5*.0565
    add_usfac = 2;
else
    add_usfac = 1;
end
base_spatial_usfac = 2; %baseline spatial up-sampling (up-sampling of filters)
spatial_usfac = add_usfac*base_spatial_usfac; %total spatial up-sampling

n_fix_inf_it = 3; %3 number of iterations to estimate fixation-based EP
n_drift_inf_it = 1; %1 number of iterations to estimate within-fixation drift

% set priors on EP (all values in deg)
fix_prior_sigma = 0.15; %prior sigma on fixation-based EP
fix_delta_noise_sigma = 0.1; %sigma on noise of eye-trackers in measuring change in EP between fixations

%this sets the posterior sigma to be ~0.004 in all cases (just a heuristic)
drift_noise_sigma = 0; %gaussian noise sigma on coil-measured drift
drift_prior_sigma = 0.004; %sigma of gaussian prior on diff of eye-pos
if all(params.use_coils) %if using both coils
    drift_prior_sigma = sqrt(0.004^2*3);
    drift_noise_sigma = sqrt(0.004^2*3);
elseif any(params.use_coils) %if just using one coil
    drift_prior_sigma = sqrt(0.004^2*2);
    drift_noise_sigma = sqrt(0.004^2*2);
end
drift_jump_sigma = 0.075; %gaussian prior sigma on change in EP between fixations during drift-corrections
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

full_nPix_us = spatial_usfac*full_nPix;
use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

%exclude data at beginning and end of each trial (numbers in sec)
beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

%saccade kernel time axis
sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_lag_set = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_lag_set);

%%
base_sp_dx = mode(expt_data.expt_dw); %size of bars in pixels
if length(unique(expt_data.expt_dw)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac/params.scale_fac; %model dx in deg

%%
all_stim_mat = decompressTernNoise(stimComp);

%% get stim-matrix for initial model estimation
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
stim_params = NIM.create_stim_params([flen full_nPix_us],'stim_dt',dt);

%pick out the predictors (of the time-embedded stimulus-matrices)
%corresponding to the central pixels we're using in the models
use_kInds = get_used_kInds([flen full_nPix_us],use_nPix_us);

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
n_blocks = length(time_data.block_flip_ids);

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

%make an indicator matrix of predictors for the block number
Xblock = zeros(fullNT,n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end
Xblock = Xblock(used_inds,:);

%% get saccade data and create saccade predictor matrices
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
saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds); %if no corresponding value was found, it was after the end of the data

used_is_blink = ET_data.is_blink(used_saccade_set); %which of used saccades are actually blinks

%identify micro and macro saccades
sac_amps = [ET_data.saccades(used_saccade_set).amplitude];
is_micro = sac_amps < 1; %identify potential microsacs (smaller than 1 degree)
big_sacs = find(~is_micro & ~used_is_blink'); %saccades
micro_sacs = find(is_micro & ~used_is_blink'); %microsaccades
all_sacs = find(~used_is_blink); %pool all saccades

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds)); %trial index of each used saccade

%create saccade/micro Xmatrices, indexing the start time of the saccade
if use_sac_kerns %could do separate predictors for big sacs and micros sacs but they have similar effects on firing rates, so easiest to pool them
    Xsac = create_event_Xmat(NT,saccade_start_inds(all_sacs),sac_lag_set,all_trialvec(used_inds));
end

%% DEFINE FIXATIONS
%index values of trial starts and stops
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

%index values of fixation starts and stops (pooling trial boundaries,
%saccades, blinks)
fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds; %duration of each fixations
assert(~any(fix_durs < 0),'alignment issues');
fix_start_inds(fix_durs==0) = []; %get rid of any cases where fixation breaks are localized with trial boundaries
fix_stop_inds(fix_durs==0) = [];
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink)); %is the given fixation following a blink?
n_fixs = length(fix_start_inds);

%compute fixation id vector
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

%% INITIALIZE STIMULUS FOR INITIAL MODEL FITTING
rand_fixpos = randn(n_fixs,1)*init_jitter_SD;
init_eyepos = zeros(NT,1);
all_stimmat_shift = all_stimmat_up;
if init_jitter_SD > 0
    for nn = 1:n_fixs %draw a random position for each fixation
        init_eyepos(used_inds(fix_ids == nn)) = rand_fixpos(nn);
    end
    %maximum initial corrections
    max_sim_pos = full_nPix*sp_dx/2; %in deg
    init_eyepos(init_eyepos > max_sim_pos) = max_sim_pos;
    init_eyepos(init_eyepos < - max_sim_pos) = -max_sim_pos;
    
    %incorporate initial eye position
    init_eyepos_rnd = round(init_eyepos/sp_dx); %in pixels
    for ii = 1:NT
        all_stimmat_shift(used_inds(ii),:) = shift_matrix_Nd(all_stimmat_up(used_inds(ii),:),-init_eyepos_rnd(used_inds(ii)),2);
    end
end
all_Xmat = NIM.create_time_embedding(all_stimmat_shift,stim_params);

%% Create set of TR and XV trials
if has_rpts
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
Robs_mat = nan(length(used_inds),tot_nUnits);
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
    
    n_block_filts = n_blocks + 1; %we add this predictor into our 'block index' Xmatrix
else
    n_block_filts = n_blocks;
end

%% create X matrices
if add_usfac > 1
    %if doing additional spatial up-sampling use tent basis functions
    X{1} = tb_proc_stim(all_Xmat(used_inds,use_kInds),add_usfac,flen);
else
    X{1} = all_Xmat(used_inds,use_kInds);
end
mod_stim_params(1) = NIM.create_stim_params([flen use_nPix_us/add_usfac],'stim_dt',dt,'stim_dx',sp_dx*add_usfac);
if model_pop_avg
    X{2} = [Xblock pop_rate]; %if using a pop-rate predictor, add it to the block-index predictors
else
    X{2} = [Xblock];
end
mod_stim_params(2) = NIM.create_stim_params([n_block_filts 1]);
if use_sac_kerns
    X{3} = Xsac;
end
mod_stim_params(3) = NIM.create_stim_params(n_sac_bins,'boundary_conds',[0 0 0]);

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
silent = 1;
sac_d2t_lambda = 100; %temporal smoothness regularization on saccade kernels
block_L2_lambda = 1; %slight penalty on block coefs to ensure convergence

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

%stimulus subunits
n_Esquared_filts = 2; %number of excitatory squared filters
n_Isquared_filts = 0; %number of inhibitory squared filters
mod_signs = [1 ones(1,n_Esquared_filts) -1*ones(1,n_Isquared_filts)];
NL_types = [{'lin'} repmat({'quad'},1,n_Esquared_filts+n_Isquared_filts)]; %NL types for full model
Xtargs = ones(1, n_Esquared_filts + n_Isquared_filts + 1);

%add block-predictor subunit
mod_signs = [mod_signs 1];
NL_types = [NL_types {'lin'}];
Xtargs = [Xtargs 2];

%initial regularization parameters (only smoothness reg on filters at
%first, used to estimate filter scales and adjust relative reg strengths)
init_d2XT = [init_lambda_d2XT*ones(n_Esquared_filts + n_Isquared_filts + 1,1); 0;];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT);

%range of smoothness regularization 'rescaling' to try
% poss_lambda_scales = logspace(-2,2,10);
poss_lambda_scales = 1;

% if ~exist(['./' mod_data_name '.mat'],'file') || recompute_init_mods == 1
all_mod_SU = zeros(tot_nUnits,1);
all_mod_SUnum = zeros(tot_nUnits,1);
[all_mod_xvLLimp,all_mod_LLimp,null_LL,null_xvLL] = deal(nan(tot_nUnits,1));
for ss = 1:tot_nUnits
    fprintf('Computing base LLs for Unit %d of %d\n',ss,tot_nUnits);
    cur_tr_inds = tr_inds(~isnan(Robs_mat(tr_inds,ss))); %set of indices where this unit was isolated
    cur_xv_inds = xv_inds(~isnan(Robs_mat(xv_inds,ss)));
    all_modfit_inds = union(cur_tr_inds,cur_xv_inds);
    tr_NT = length(cur_tr_inds);
    xv_NT = length(cur_xv_inds);
    
    if ss > n_probes && ss <= length(su_probes) + n_probes
        all_mod_SU(ss) = su_probes(ss-n_probes);
        all_mod_SUnum(ss) = SU_numbers(ss-n_probes);
    end
    if ~isempty(cur_tr_inds)
        Robs = Robs_mat(:,ss);
        
        %estimate null model
        if use_sac_kerns
            null_mod = NIM(mod_stim_params,repmat({'lin'},1,2),[1 1],'Xtargets',[2 3],...
                'd2t',[0 sac_d2t_lambda],'l2',[block_L2_lambda 0]);
        else
            null_mod = NIM(mod_stim_params,'lin',1,'Xtargets',2); %otherwise null model just has block predictors
        end
        null_mod = null_mod.fit_filters(Robs,X,cur_tr_inds,'silent',silent);
        null_xvLL(ss) = null_mod.eval_model(Robs,X,cur_xv_inds);
        
        %% estimate initial models and adjust base reg strength
        init_gqm = NIM(mod_stim_params,NL_types,mod_signs,'Xtargets',Xtargs,...
            'd2xt',init_lambda_d2XT*(Xtargs==1),'l2',block_L2_lambda*(Xtargs == 2));
        init_gqm = init_gqm.fit_filters(Robs,X,all_modfit_inds,'silent',silent);
        
        %adjust regularization stregnths of each stimulus filter relative to the scale of filter outputs
        [~,~,mod_internals] = init_gqm.eval_model(Robs,X,all_modfit_inds);
        new_lambda_d2xt = base_lambda_d2XT./var(mod_internals.gint(:,Xtargs == 1));
        new_lambda_l1 = base_lambda_L1./std(mod_internals.gint(:,Xtargs == 1));
        init_gqm = init_gqm.set_reg_params('sub_inds',find(Xtargs == 1),'d2xt',new_lambda_d2xt,'l1',new_lambda_l1);
        
        %% scan a range of smoothness regularization strengths and pick the best using xval
        cur_lambda_d2XT = init_gqm.get_reg_lambdas('sub_inds',find(Xtargs == 1),'d2xt');
        new_gqm = init_gqm;
        new_gqm_xvLL = nan(length(poss_lambda_scales),1);
        clear all_gqms
        for ll = 1:length(poss_lambda_scales) %loop over possible range of smoothness lambdas
            fprintf('%d ',ll);
            new_lambda_d2XT = cur_lambda_d2XT*poss_lambda_scales(ll);
            new_gqm = new_gqm.set_reg_params('sub_inds',find(Xtargs == 1),'d2xt',new_lambda_d2XT);
            new_gqm = new_gqm.fit_filters(Robs,X,cur_tr_inds,'silent',silent);
            new_gqm_xvLL(ll) = new_gqm.eval_model(Robs,X,cur_xv_inds);
            all_gqms(ll) = new_gqm;
        end
        [best_xvLL,best_scale] = max(new_gqm_xvLL);
        fprintf('. Best at %d\n',best_scale);
        new_gqm = all_gqms(best_scale);
        up_gqm = new_gqm;
        
        
        %% add saccade kernels and fit the model
        if use_sac_kerns
            %                 up_gqm = up_gqm.add_subunits('lin',[1 1],'xtargs',[3 4],'d2t',sac_d2t_lambda);
            up_gqm = up_gqm.add_subunits('lin',[1],'xtargs',[3],'d2t',sac_d2t_lambda);
            up_gqm = up_gqm.fit_filters(Robs,X,cur_tr_inds,'silent',silent);
        end
        xvLL = up_gqm.eval_model(Robs,X,cur_xv_inds);
        all_mod_xvLLimp(ss) = (xvLL - null_xvLL(ss))/log(2);
        
        %% refit models using all data
        up_gqm = up_gqm.fit_filters(Robs,X,all_modfit_inds,'silent',silent);
        null_mod = null_mod.fit_filters(Robs,X,all_modfit_inds,'silent',silent);
        
        %% store results
        all_nullmod(ss) = null_mod;
        all_mod_fits(ss) = up_gqm;
        all_mod_fits_withspkNL(ss) = up_gqm.fit_spkNL(Robs,X,all_modfit_inds,'silent',silent); %fit spkNL params
        null_LL(ss) = null_mod.eval_model(Robs,X,all_modfit_inds);
        LL = all_mod_fits_withspkNL(ss).eval_model(Robs,X,all_modfit_inds);
        all_mod_LLimp(ss) = (LL - null_LL(ss))/log(2); %LL improvement in bits per spike
        
    end
end

%%
%     save(mod_data_name,'all_mod*','all_nullmod','Clust_data','null_xvLL','null_LL','*_trials');
% else
%     fprintf('Loading pre-computed initial models\n');
%     load(mod_data_name);

%%
% end

%%
%if using jittered initialization, reconstruct unjittered stimulus
if init_jitter_SD > 0
    all_Xmat = NIM.create_time_embedding(all_stimmat_up,stim_params);
end
all_Xmat = all_Xmat(used_inds,:);

%% SELECT USABLE UNITS AND make Robs_mat
% usable_units = find(all_mod_xvLLimp > 0);
tr_units = find(~isnan(all_mod_xvLLimp)); %use all units we have data for

n_used_sus = sum(all_mod_SU(tr_units) ~= 0);
n_used_mus = sum(all_mod_SU(tr_units) == 0);
fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
n_tr_chs = length(tr_units);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = Robs_mat(:,tr_units);

all_modfit_inds = union(tr_inds,xv_inds);

%% set of units to perform LOO xval on
if use_LOOXV == 2
    loo_set = 1:length(tr_set);
elseif use_LOOXV == 1
    loo_set = find(all_mod_SU(tr_set) > 0);
else
    loo_set = [];
end

%% DIAGNOSTICS
full_prates = nan(NT,n_tr_chs);
full_nullrates = nan(NT,n_tr_chs);
for cc = 1:n_tr_chs
    [~,full_prates(:,cc)] = all_mod_fits_withspkNL(tr_units(cc)).eval_model(Robs_mat(:,cc),X);
    [~,full_nullrates(:,cc)] = all_nullmod(tr_units(cc)).eval_model(Robs_mat(:,cc),X);
end
full_modLL = Robs_mat.*log(full_prates) - full_prates;
full_nullLL = Robs_mat.*log(full_nullrates) - full_nullrates;
full_LLimp = full_modLL-full_nullLL; %pointwise LL improvement over null model

block_LLimp = nan(n_blocks,n_tr_chs);
for tt = 1:n_blocks
    cur_used_inds = find(all_blockvec(used_inds) == tt);
    block_LLimp(tt,:) = nanmean(full_LLimp(cur_used_inds,:)); %compute within-block avg LL imps
end

block_LLimp_zscore = nanzscore(block_LLimp);
avg_LLimp_zscore = nanmean(block_LLimp_zscore,2);
questionable_blocks = find(avg_LLimp_zscore <= -1);
if ~isempty(questionable_blocks)
    fprintf('Warning, found %d questionable blocks\n',length(questionable_blocks));
    fprintf('Avg z-score LL: %.3f\n',avg_LLimp_zscore(questionable_blocks));
end

%% COMPUTE FORWARD-PROJECTED FIXATION INDICES
%push the effects of saccades forward in time for inferring eye positions
sac_shift = round(0.05/dt); %how much to shift jump points forward in time
pfix_start_inds = fix_start_inds;
pfix_stop_inds = fix_stop_inds;
for i = 1:length(fix_start_inds)
    next_trial = trial_start_inds(find(trial_start_inds >= fix_start_inds(i),1,'first'));
    if next_trial > fix_start_inds(i) + sac_shift %if the forward projection does not push you into another trial
        pfix_start_inds(i) = fix_start_inds(i) + sac_shift;
    end
    next_trial = trial_start_inds(find(trial_start_inds >= fix_stop_inds(i),1,'first'));
    if next_trial > fix_stop_inds(i) + sac_shift
        pfix_stop_inds(i) = fix_stop_inds(i) + sac_shift;
    end
end

%compute forward-projected fixation ids
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
    pfix_ids(cur_inds) = ii;
end

%% generate shift matrices. Must be applied to the stimulus (not the filters)
%shifts for inferring fixation corrections
max_fix_shift_pix = round(max_fix_shift/sp_dx); %max eye pos deviation (deg)
fix_shifts = -max_fix_shift_pix:max_fix_shift_pix; %shift-range
n_fix_shifts = length(fix_shifts); %number of shifts
fix_zero_shift = find(fix_shifts==0);

%shifts for inferring drift
max_drift_shift_pix = round(max_drift_shift/sp_dx);
drift_shifts = -max_drift_shift_pix:max_drift_shift_pix;
n_drift_shifts = length(drift_shifts);
drift_zero_shift = find(drift_shifts==0);

It = speye(flen);

fix_shift_mats = cell(n_fix_shifts,1);
for xx = 1:n_fix_shifts
    temp = spdiags( ones(full_nPix_us,1), -fix_shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    fix_shift_mats{xx} = temp(:,use_kInds);
end

drift_shift_mats = cell(n_drift_shifts,1);
for xx = 1:n_drift_shifts
    temp = spdiags( ones(full_nPix_us,1), -drift_shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    drift_shift_mats{xx} = temp(:,use_kInds);
end

%% if using eye-tracking data get relevant signals
if any(params.use_coils)
    measured_eyepos = [corrected_eye_vals_interp(:,2) corrected_eye_vals_interp(:,4)]; %orthogonal eye position from two coils
    
    %threshold measured eye positions to avoid extreme values
    max_measured_pos = 1;
    measured_eyepos(measured_eyepos > max_measured_pos) = max_measured_pos;
    measured_eyepos(measured_eyepos < -max_measured_pos) = -max_measured_pos;
    measured_eyepos(isnan(measured_eyepos)) = 0;
    
    %smooth out fast transients in eye signal during fixations
    eye_smooth_sig = round(0.025/dt);
    interp_inds = [];
    measured_drift = nan(length(used_inds),2);
    measured_fix_med = nan(n_fixs,2);
    for ii = 1:n_fixs
        cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
        if length(cur_inds) > eye_smooth_sig*5; %if theres enough data in this fixation to smooth
            measured_eyepos(used_inds(cur_inds),1) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),1),eye_smooth_sig,2);
            measured_eyepos(used_inds(cur_inds),2) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),2),eye_smooth_sig,2);
        end
        measured_fix_med(ii,:) = median(measured_eyepos(cur_inds,:));
        measured_drift(cur_inds,:) = bsxfun(@minus,measured_eyepos(cur_inds,:),measured_fix_med(ii,:));
        interp_inds = [interp_inds; cur_inds'];
    end
    interp_inds = unique(interp_inds);
    
    %interpolate eye position for data between fixations (saccades/blinks)
    measured_eyepos(used_inds,:) = interp1(used_inds(interp_inds),measured_eyepos(used_inds(interp_inds),:),used_inds);
    measured_eyepos(isnan(measured_eyepos)) = 0;
    measured_eyepos = measured_eyepos(used_inds,:);
    
    forward_measured_drift = nan(size(measured_drift));
    %forward-project measured drift so it aligns with neural data
    for ii = 1:n_fixs
        cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
        cur_pinds = cur_inds + sac_shift;
        if length(cur_inds) > sac_shift
            forward_measured_drift(cur_inds(sac_shift+1:end),:) = measured_drift(cur_inds(1:(end-sac_shift)),:);
        end
    end
    measured_drift = [zeros(1,2); diff(forward_measured_drift)];
    clear forward_measured_drift
end

%% define updated posteriors (incorporating measured coil data)
measured_fix_deltas = nan(n_fixs,1); %measured change in eye position between adjacent fixations

if all(params.use_coils) == 0 %if not using any coil signals
    post_drift_var = drift_prior_sigma^2; %posterior var after observing coil data is still the prior
    post_drift_mean = zeros(NT,1);
    n_coil_samps = 0;
else
    if all(params.use_coils == 1) %if using both coils
        measured_fix_deltas(2:end) = mean(diff(measured_fix_med),2); %avg change in eye position between fixation
        measured_drift = mean(measured_drift,2); %avg of measured drifts
        n_coil_samps = 2;
    else %if only using one coil
        measured_fix_deltas(2:end) = diff(measured_fix_med(:,params.use_coils == 1)); %just use the good coil
        measured_drift = measured_drift(:,params.use_coils == 1);
    end
    fix_delta_sigma = fix_delta_noise_sigma/sqrt(n_coil_samps); %if we are avging multiple measures of delta_EP, scale noise variance of the measure (assuming independent measurements)
    post_drift_var = 1/(n_coil_samps/drift_noise_sigma^2 + 1/drift_prior_sigma^2); %variance of posterior gaussian, given N_coil_samps observations of drift, and gaussian prior
    post_mean_drift = n_coil_samps*measured_drift*post_drift_var/drift_noise_sigma^2; %mean of the gaussian posterior
end
post_drift_sigma = sqrt(post_drift_var);

%prior on fixations
fix_Lprior = -(fix_shifts*sp_dx).^2./(2*fix_prior_sigma^2); %log-prior on within-fixation eye-pos
fix_Lprior = fix_Lprior - logsumexp(fix_Lprior); %normalize

%prior on drift jumps
drift_jump_Lprior = -(drift_shifts*sp_dx).^2./(2*drift_jump_sigma^2); %log-prior on initial position (after a jump) during drift-inference
drift_jump_Lprior = drift_jump_Lprior - logsumexp(drift_jump_Lprior);

%prior matrix for drift (if not using eye position)
cdist = squareform(pdist(drift_shifts'*sp_dx));
base_LA = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2); %baseline log-prior matrix on drift position changes
base_LA = bsxfun(@minus,base_LA,logsumexp(base_LA,2)); %normalize

%%
for ss = 1:length(all_mod_fits)
    all_mod_fits(ss).fit_props.null_LL = null_LL(ss);
end
%%
fix_boundaries = [fix_start_inds fix_stop_inds];
pfix_boundaries = [pfix_start_inds pfix_stop_inds];
trial_boundaries = [trial_start_inds trial_end_inds];
fix_data = struct('fix_boundaries',fix_boundaries,'pfix_boundaries',pfix_boundaries,'fix_post_blink',fix_post_blink,...
    'trial_boundaries',trial_boundaries);

fit_data = struct('tr_units',tr_units,'stim_params',stim_params,'used_inds',used_inds,'use_kInds',use_kInds,'modfit_inds',all_modfit_inds);
HMM_data = struct('add_usfac',add_usfac,'fix_delta_noise_sigma',fix_delta_noise_sigma,'fix_shifts',fix_shifts,...
    'n_fix_iters',n_fix_inf_it,'fix_Lprior',fix_Lprior,'use_coils',params.use_coils,'stim_dx',sp_dx,...
    'drift_dsf',drift_dsf,'drift_jump_Lprior',drift_jump_Lprior,'drift_shifts',drift_shifts,...
    'post_drift_sigma',post_drift_sigma,'n_drift_iters',n_drift_inf_it,'base_LA',base_LA);
HMM_data.fix_shift_mats = fix_shift_mats;
HMM_data.drift_shift_mats = drift_shift_mats;

ET_meas = struct('fix_deltas',measured_fix_deltas,'post_drift_mean',post_drift_mean);
[mod_fits,post_mean_EP,post_std_EP] = eyetracking_EM(all_mod_fits(tr_units),Robs_mat(:,tr_units),all_stimmat_up,X(2:end),...
    fit_data,HMM_data,ET_meas,fix_data);

%% ITERATE FIXATION-BASED CORRECTIONS
Xtargs = [all_mod_fits(1).subunits(:).Xtarg];
stim_mod_signs = [all_mod_fits(1).subunits(find(Xtargs == 1)).weight];
stim_NL_types = all_mod_fits(1).get_NLtypes(find(Xtargs == 1));
n_stim_subs = sum(Xtargs == 1);

it_mods{1} = all_mod_fits;
it_mods_spkNL{1} = all_mod_fits_withspkNL;
it_LLimp(1,:) = all_mod_LLimp;

it_fix_post_mean = nan(n_fix_inf_it,n_fixs);
it_fix_post_std = nan(n_fix_inf_it,n_fixs);

% if use_LOOXV > 0
%     it_LLimp_LOO(length(loo_set),1,:) = all_mod_LLimp;
%     for xv = 1:length(loo_set)
%         it_mods_LOO{xv,1} = all_mod_fits;
%         it_mods_spkNL_LOO{xv,1} = all_mod_fits_withspkNL;
%     end
% end
for nn = 1:n_fix_inf_it
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,n_fix_inf_it);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(klen_us/add_usfac,sum(Xtargs==1),n_tr_chs); %compile all stimulus filters
    lin_kerns = nan(n_block_filts,n_tr_chs); %compile all linear filters (except sac kernels)
    if use_sac_kerns %for storing all saccade kernels
        sac_kerns = nan(n_sac_bins,n_tr_chs);
    end
    mod_spkNL_params = nan(3,n_tr_chs); %assuming softplus spkNL (with three params including offset)
    for ss = 1:n_tr_chs
        %         filt_bank(ss,:,:) = cell2mat(it_mods{nn}(tr_units(ss)).get_filtKs(find(Xtargs == 1))');
        filt_bank(:,:,ss) = cell2mat(it_mods{nn}(tr_units(ss)).get_filtKs(find(Xtargs == 1))');
        mod_spkNL_params(:,ss) = [it_mods_spkNL{nn}(tr_units(ss)).spkNL.params it_mods_spkNL{nn}(tr_units(ss)).spkNL.theta];
        lin_kerns(:,ss) = cell2mat(it_mods{nn}(tr_units(ss)).get_filtKs(find(Xtargs == 2))');
        if use_sac_kerns
            sac_kerns(:,ss) = cell2mat(it_mods{nn}(tr_units(ss)).get_filtKs(find(Xtargs == 3))');
        end
    end
    
    %get net outputs of all non-stimulus predictors
    if model_pop_avg
        block_out = [Xblock pop_rate]*lin_kerns;
    else
        block_out = Xblock*lin_kerns;
    end
    if use_sac_kerns
        sac_out = Xsac*sac_kerns;
    end
    
    %% ESTIMATE LL for each shift in each stimulus frame
    cur_Xmat = all_Xmat;
    %precompute LL at all shifts for all units
    frame_LLs = nan(NT,n_fix_shifts);
    reverseStr = '';
    for xx = 1:length(fix_shifts)
        msg = sprintf('Calculating LLs for fixation shift %d of %d\n',xx,n_fix_shifts);
        fprintf([reverseStr msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        cur_stim_shift = cur_Xmat*fix_shift_mats{xx};
        
        %process with spatial TB if doing additional up-sampling
        if add_usfac > 1
            cur_stim_shift = tb_proc_stim(cur_stim_shift,add_usfac,flen);
        end
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,n_tr_chs);
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(3,:)); %add in constant offset terms
        for ff = 1:n_stim_subs
            cur_filt_outs = cur_stim_shift*squeeze(filt_bank(:,ff,:));
            if ~strcmp(stim_NL_types{ff},'lin')
                cur_filt_outs = all_mod_fits(1).subunits(ff).apply_NL(cur_filt_outs);
            end
            gfuns = gfuns + stim_mod_signs(ff)*cur_filt_outs;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + block_out;
        if use_sac_kerns
            gfuns = gfuns + sac_out;
        end
        
        %incorporate beta
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(1,:));
        
        %handle numerical overflow with log(1+exp)
        too_large = gfuns > 50;
        pred_rate = log(1+exp(gfuns));
        pred_rate(too_large) = gfuns(too_large);
        
        %incorporate alpha
        pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(2,:));
        
        %enforce min predicted rate
        pred_rate(pred_rate < 1e-50) = 1e-50;
        
        frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
    end
    
    %% INFER MICRO-SAC SEQUENCE
    fix_LLs = nan(n_fixs,n_fix_shifts);
    for ii = 1:n_fixs %compute total LL of each eye position for all data within each fixation
        cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii); %use forward-projected fixation indices to account for neural delay
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:)); %sum LLs for data within the fixation
    end
    
    if all(params.use_coils==0) %if not using coil info we just compute the posterior at each fixation independently
        Lgamma = bsxfun(@plus,fix_LLs,fix_Lprior); %add log-likelihood and log-prior
        Lgamma = bsxfun(@minus,Lgamma,logsumexp(Lgamma,2)); %normalize
    else
        Lalpha=zeros(n_fixs,n_fix_shifts); %log of forward messages
        Lbeta = zeros(n_fixs,n_fix_shifts); %log of backward messages
        Lalpha(1,:) = fix_Lprior + fix_LLs(1,:); %initialize using prior
        for t=2:n_fixs %compute forward messages
            if ~fix_post_blink(t) %if this fixation is NOT following a blink
                %make quadratic transition matrix centered on coil-measured deltaX, with variance
                %given by the delta_noise_sigma^2
                cdist = pdist2(fix_shifts'*sp_dx + measured_fix_deltas(t),fix_shifts'*sp_dx);
                cur_LA = -cdist.^2/(2*fix_delta_noise_sigma^2);
                cur_LA = bsxfun(@plus,cur_LA,fix_Lprior); %multiply this by the prior over fixation positions
                cur_LA = bsxfun(@minus,cur_LA,logsumexp(cur_LA,2)); %normalize log transition-matrix
                Lalpha(t,:) = logmulexp(Lalpha(t-1,:),cur_LA) + fix_LLs(t,:); %compute forward message
            else %if it is following a blink, ignore the eye-trackers
                Lalpha(t,:) = fix_Lprior + fix_LLs(t,:);
            end
        end
        
        Lbeta(n_fixs,:)=zeros(1,n_fix_shifts); %initialize log-backward message
        for t=n_fixs-1:-1:1%compute backward messages
            if ~fix_post_blink(t+1) %if the fixation isn't followed by a blink
                %construct A contingent on measured deltX between current and NEXT fixation
                cdist = pdist2(fix_shifts'*sp_dx + measured_fix_deltas(t+1),fix_shifts'*sp_dx);
                cur_LA = -cdist.^2/(2*fix_delta_noise_sigma^2);
                cur_LA = bsxfun(@plus,cur_LA,fix_Lprior); %multiply by log-prior
                cur_LA = bsxfun(@minus,cur_LA,logsumexp(cur_LA,2)); %normalize
                Lbeta(t,:) = logmulexp(Lbeta(t+1,:) + fix_LLs(t+1,:),cur_LA'); %compute backward message
            else %otherwise, dont incorporate coil-measured data
                Lbeta(t,:) = zeros(1,n_fix_shifts);
            end
        end
        Lgamma = Lalpha + Lbeta;%combine forward and backward messages
        Lgamma = bsxfun(@minus,Lgamma,logsumexp(Lgamma,2)); %normalize
    end
    gamma = exp(Lgamma); %posterior on fixation positions
    
    %compute posterior mean and SD
    it_fix_post_mean(nn,:) = sum(bsxfun(@times,gamma,fix_shifts),2);
    cur_diff = bsxfun(@minus,it_fix_post_mean(nn,:)',fix_shifts).^2;
    it_fix_post_std(nn,:) = sqrt(sum(cur_diff.*gamma,2));
    
    %back-project saccade-times and interpolate over saccades
    all_fix_post_mean_cor = nan(NT,1);
    all_fix_post_mean_cor(~isnan(fix_ids)) = it_fix_post_mean(nn,fix_ids(~isnan(fix_ids)));
    all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT); %interpolate position between fixations
    all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;
    
    %% adjust stimulus matrix for fixation eye positions
    cur_fix_shifts = round(all_fix_post_mean_cor);
    all_shift_stimmat_up = all_stimmat_up;
    for ii=1:NT
        all_shift_stimmat_up(used_inds(ii),:) = shift_matrix_Nd(all_stimmat_up(used_inds(ii),:),-cur_fix_shifts(ii),2);
    end
    X{1} = NIM.create_time_embedding(all_shift_stimmat_up,stim_params);
    
    %% REFIT ALL CELLS
    if add_usfac %project onto spatial TBs if needed
        X{1} = tb_proc_stim( X{1}(used_inds,use_kInds),add_usfac,flen);
    else
        X{1} =  X{1}(used_inds,use_kInds);
    end
    
    silent = 1;
    reverseStr = '';
    for ss = 1:length(tr_units)
        cur_cell = tr_units(ss);
        cur_unit_ind = find(tr_units == cur_cell);
        cur_fit_inds = all_modfit_inds(~isnan(Robs_mat(all_modfit_inds,cur_unit_ind))); %all non-rpt indices when current unit was isolated
        
        it_mods{nn+1}(cur_cell) = it_mods{nn}(cur_cell);
        it_mods{nn+1}(cur_cell) = it_mods{nn+1}(cur_cell).fit_filters(Robs_mat(:,cur_unit_ind),...
            X,cur_fit_inds,'silent',silent); %fit stimulus filters
        
        %refit spk NL
        it_mods_spkNL{nn+1}(cur_cell) = it_mods{nn+1}(cur_cell).fit_spkNL(Robs_mat(:,cur_unit_ind),X,cur_fit_inds,'silent',silent);
        newLL = it_mods_spkNL{nn+1}(cur_cell).eval_model(Robs_mat(:,cur_unit_ind),X,cur_fit_inds);
        
        it_LLimp(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
        if nn > 1
            fprintf('Cell %d.  Original: %.4f  Prev: %.4f  New: %.4f\n',ss,all_mod_LLimp(cur_cell),it_LLimp(nn,cur_cell),it_LLimp(nn+1,cur_cell));
        else
            fprintf('Cell %d. Original: %.4f  New: %.4f\n',ss,all_mod_LLimp(cur_cell),it_LLimp(nn+1,cur_cell));
        end
    end
end

%% NOW INFER DRIFT CORRECTIONS
dit_mods{1} = it_mods{n_fix_inf_it+1};
dit_mods_spkNL{1} = it_mods_spkNL{n_fix_inf_it+1};
dit_LLimp(1,:) = it_LLimp(n_fix_inf_it+1,:);
if use_LOOXV > 0
    for xv = 1:length(loo_set)
        dit_mods_LOO{xv,1} = it_mods_LOO{xv,n_fix_inf_it+1};
        dit_mods_spkNL_LOO{xv,1} = it_mods_spkNL_LOO{xv,n_fix_inf_it+1};
    end
    dit_LLimp_LOO(:,1,:) = it_LLimp_LOO(:,n_fix_inf_it+1,:);
end
for nn = 1:n_drift_inf_it
    fprintf('Inferring drift corrections, iter %d of %d\n',nn,n_drift_inf_it);
    if use_LOOXV > 0
        %re-create all-cell fixation-corrected Xmat
        all_fix_post_mean_cor(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
        all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
        all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;
        all_fix_post_mean_cor = round(all_fix_post_mean_cor) + max_Tshift + 1;
        for i=1:NT
            all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_fix_post_mean_cor(i)};
        end
        all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
        all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
    elseif nn == 1
        all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
    end
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_tr_chs,klen_us/add_usfac,n_squared_filts+1);
    lin_kerns = nan(n_tr_chs,n_block_filts);
    if use_sac_kerns
        sac_kerns = nan(n_tr_chs,n_sac_bins);
        msac_kerns = nan(n_tr_chs,n_sac_bins);
    end
    mod_spkNL_params = nan(n_tr_chs,3);
    for ss = 1:n_tr_chs
        cur_Xtargs = [dit_mods{nn}(tr_set(ss)).mods(:).Xtarget];
        cur_k = [dit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = dit_mods_spkNL{nn}(tr_set(ss)).spk_NL_params(1:3);
        cur_lin_kerns = dit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
        lin_kerns(ss,1:length(cur_lin_kerns)) = cur_lin_kerns;
        %         if model_pop_avg && tr_set(ss) == tot_nUnits
        %            lin_kerns(ss,end) = 0; %not using pop-avg predictor for the pop-avg itself
        %         end
        if use_sac_kerns
            sac_kerns(ss,:) = dit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = dit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
        end
    end
    filt_bank = permute(filt_bank,[2 1 3]);
    
    %indicator predictions
    if model_pop_avg
        block_out = [Xblock(used_inds,:) pop_rate]*lin_kerns';
    else
        block_out = Xblock(used_inds,:)*lin_kerns';
    end
    if use_sac_kerns
        sac_out = Xsac*sac_kerns';
        msac_out = Xmsac*msac_kerns';
    end
    
    %% ESTIMATE LL for each shift in each stimulus frame
    
    frame_LLs = nan(NT,n_Dshifts);
    for xx = 1:length(Dshifts)
        fprintf('Dshift %d of %d\n',xx,n_Dshifts);
        cur_stim_shift = all_Xmat_up_fixcor*Dshift_mat{xx};
        
        if add_usfac > 1
            cur_stim_shift = tb_proc_stim(cur_stim_shift,add_usfac,flen);
        end
        
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
    
    Lgamma = nan(NT,n_Dshifts);
    for ff = 1:n_fixs
        if mod(ff,100)==0
            fprintf('Fixation %d of %d\n',ff,n_fixs);
        end
        
        tset = find(pfix_ids==ff)'; %indices in this projected fixation
        ntset = length(tset);
        if ntset > drift_dsf
            nt_pts = ceil(ntset/drift_dsf); %number of points to infer drift at
            tset_inds = 1+floor((0:(ntset-1))/drift_dsf);
            talpha=zeros(nt_pts,n_Dshifts);
            tbeta = zeros(nt_pts,n_Dshifts);
            
            tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf); %drift inference locations
            tpt_loc(end) = ntset;
            
            cur_drift_mean = post_mean_drift(tset);
            cur_LL_set = frame_LLs(tset,:);
            if mod(ntset,drift_dsf) ~= 0 %deal with any 'dangling points'
                dangling_pts = nt_pts*drift_dsf-ntset;
                cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_Dshifts));
                cur_drift_mean = cat(1,cur_drift_mean,nan(dangling_pts,1));
            end
            %sum LL over points within a downsampled bin
            cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_Dshifts]);
            cur_LL_set = squeeze(sum(cur_LL_set,1));
            
            cur_drift_mean = drift_dsf*nanmean(reshape(cur_drift_mean,[drift_dsf nt_pts]));
            
            %get transition matrix
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
                    Lgamma(tset,:) = int_gamma;
                else
                    Lgamma(tset,:) = repmat(temp_gamma,ntset,1);
                end
            else
                Lgamma(tset,:) = temp_gamma;
            end
        end
    end
    Lgamma = bsxfun(@minus,Lgamma,logsumexp(Lgamma,2));
    gamma = exp(Lgamma);
    drift_post_mean(nn,:) = sum(bsxfun(@times,gamma,Dshifts),2);
    drift_post_std(nn,:) = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2) - drift_post_mean(nn,:)'.^2);
    
    drift_post_mean(nn,isnan(drift_post_mean(nn,:))) = 0;
    drift_post_std(nn,isnan(drift_post_std(nn,:))) = 0;
    
    %interpolate over saccades
    drift_post_mean(nn,:) = interp1(find(~isnan(pfix_ids)),drift_post_mean(nn,~isnan(pfix_ids)),1:NT);
    drift_post_std(nn,:) = interp1(find(~isnan(pfix_ids)),drift_post_std(nn,~isnan(pfix_ids)),1:NT);
    drift_post_mean(nn,isnan(drift_post_mean(nn,:))) = 0;
    
    %% construct drift-corrected X-mat
    fix_post_cor = nan(NT,1);
    fix_post_cor(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
    fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
    fix_post_cor(isnan(fix_post_cor)) = 0;
    drift_post_cor = squeeze(drift_post_mean(nn,:));
    
    for ii = 1:length(trial_start_inds)
        cur_inds = trial_start_inds(ii):trial_end_inds(ii);
        drift_post_cor(cur_inds(1:end-sac_shift)) = drift_post_cor(cur_inds(sac_shift+1:end));
    end
    drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids)),1:NT);
    drift_post_cor(isnan(drift_post_cor)) = 0;
    all_post_cor = round((fix_post_cor+drift_post_cor)) + max_Tshift + 1;
    
    %RECOMPUTE XMAT
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
    end
    X{1} = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    if add_usfac > 1
        X{1} = tb_proc_stim(X{1}(used_inds,use_kInds_up),add_usfac,flen);
    else
        X{1} = X{1}(used_inds,use_kInds_up);
    end
    %% REFIT ALL CELLS
    silent = 1;
    for ss = 1:length(tr_set)
        cur_cell = tr_set(ss);
        fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
        cur_unit_ind = find(tr_set == cur_cell);
        cur_fit_inds = all_modfit_inds(~isnan(Robs_mat(all_modfit_inds,cur_unit_ind)));
        
        %         if tr_set(ss) == tot_nUnits && model_pop_avg
        %             X{2} = Xblock(used_inds,:); %for modeling pop-avg dont use the pop-avg as a predictor
        %         end
        
        dit_mods{nn+1}(cur_cell) = dit_mods{nn}(cur_cell);
        dit_mods{nn+1}(cur_cell) = NMMfit_filters(dit_mods{nn+1}(cur_cell),Robs_mat(:,cur_unit_ind),...
            X,[],cur_fit_inds,silent); %fit stimulus filters
        
        dit_mods_spkNL{nn+1}(cur_cell) = NMMfit_logexp_spkNL(dit_mods{nn+1}(cur_cell),Robs_mat(:,cur_unit_ind),X,[],cur_fit_inds);
        
        newLL = NMMeval_model(dit_mods_spkNL{nn+1}(cur_cell),Robs_mat(:,cur_unit_ind),X,[],cur_fit_inds);
        dit_LLimp(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
        
        fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(cur_cell),dit_LLimp(nn,cur_cell),dit_LLimp(nn+1,cur_cell));
    end
    
    %     if model_pop_avg
    %           X{2} = [Xblock(used_inds,:) pop_rate]; %for modeling pop-avg dont use the pop-avg as a predictor
    %     end
    %%
    if use_LOOXV > 0
        for xv = 1:length(loo_set)
            fprintf('Inferring drift corrections, XV %d of %d\n',xv,length(loo_set));
            
            %re-create LOOXV fixation-corrected Xmat
            all_fix_post_mean_cor(~isnan(fix_ids)) = squeeze(it_fix_post_mean_LOO(xv,end,fix_ids(~isnan(fix_ids))));
            all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
            all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;
            all_fix_post_mean_cor = round(all_fix_post_mean_cor) + max_Tshift + 1;
            for i=1:NT
                all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_fix_post_mean_cor(i)};
            end
            all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
            all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
            
            %% PREPROCESS MODEL COMPONENTS
            cur_uset = setdiff(1:n_tr_chs,loo_set(xv));
            n_uset = length(cur_uset);
            filt_bank = zeros(n_uset,klen_us/add_usfac,n_squared_filts+1);
            lin_kerns = nan(n_uset,n_block_filts);
            if use_sac_kerns
                sac_kerns = nan(n_uset,n_sac_bins);
                msac_kerns = nan(n_uset,n_sac_bins);
            end
            mod_spkNL_params = nan(n_uset,3);
            for ss = 1:n_uset
                cur_Xtargs = [dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(:).Xtarget];
                cur_k = [dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 1).filtK];
                n_used_filts = size(cur_k,2);
                filt_bank(ss,:,1:n_used_filts) = cur_k;
                mod_spkNL_params(ss,:) = dit_mods_spkNL_LOO{xv,nn}(tr_set(cur_uset(ss))).spk_NL_params(1:3);
                cur_lin_kerns = dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 2).filtK;
                lin_kerns(ss,1:length(cur_lin_kerns)) = cur_lin_kerns;
                %                 if model_pop_avg && tr_set(ss) == tot_nUnits
                %                     lin_kerns(ss,end) = 0; %not using pop-avg predictor for the pop-avg itself
                %                 end
                if use_sac_kerns
                    sac_kerns(ss,:) = dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 3).filtK;
                    msac_kerns(ss,:) = dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 4).filtK;
                end
            end
            filt_bank = permute(filt_bank,[2 1 3]);
            
            %indicator predictions
            if model_pop_avg
                block_out = [Xblock(used_inds,:) pop_rate]*lin_kerns';
            else
                block_out = Xblock(used_inds,:)*lin_kerns';
            end
            if use_sac_kerns
                sac_out = Xsac*sac_kerns';
                msac_out = Xmsac*msac_kerns';
            end
            %% ESTIMATE LL for each shift in each stimulus frame
            frame_LLs = nan(NT,n_Dshifts);
            for xx = 1:length(Dshifts)
                fprintf('Dshift %d of %d\n',xx,n_Dshifts);
                cur_stim_shift = all_Xmat_up_fixcor*Dshift_mat{xx};
                
                if add_usfac > 1
                    cur_stim_shift = tb_proc_stim(cur_stim_shift,add_usfac,flen);
                end
                
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
            
            %% INFER DRIFT SEQUENCE
            Lgamma = nan(NT,n_Dshifts);
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
                            Lgamma(tset,:) = int_gamma;
                        else
                            Lgamma(tset,:) = repmat(temp_gamma,ntset,1);
                        end
                    else
                        Lgamma(tset,:) = temp_gamma;
                    end
                end
            end
            Lgamma = bsxfun(@minus,Lgamma,logsumexp(Lgamma,2));
            
            gamma = exp(Lgamma);
            drift_post_mean_LOO(xv,nn,:) = sum(bsxfun(@times,gamma,Dshifts),2);
            drift_post_std_LOO(xv,nn,:) = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2) - squeeze(drift_post_mean_LOO(xv,nn,:)).^2);
            
            drift_post_mean_LOO(xv,nn,isnan(drift_post_mean_LOO(xv,nn,:))) = 0;
            drift_post_std_LOO(xv,nn,isnan(drift_post_std_LOO(xv,nn,:))) = 0;
            
            drift_post_mean_LOO(xv,nn,:) = interp1(find(~isnan(pfix_ids)),squeeze(drift_post_mean_LOO(xv,nn,~isnan(pfix_ids))),1:NT);
            drift_post_std_LOO(xv,nn,:) = interp1(find(~isnan(pfix_ids)),squeeze(drift_post_std_LOO(xv,nn,~isnan(pfix_ids))),1:NT);
            
            %% construct drift-corrected X-mat
            fix_post_cor = nan(NT,1);
            fix_post_cor(~isnan(fix_ids)) = squeeze(it_fix_post_mean_LOO(xv,end,fix_ids(~isnan(fix_ids))));
            fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
            fix_post_cor(isnan(fix_post_cor)) = 0;
            drift_post_cor = squeeze(drift_post_mean_LOO(xv,nn,:));
            
            for ii = 1:length(trial_start_inds)
                cur_inds = trial_start_inds(ii):trial_end_inds(ii);
                drift_post_cor(cur_inds(1:end-sac_shift)) = drift_post_cor(cur_inds(sac_shift+1:end));
            end
            drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids)),1:NT);
            drift_post_cor(isnan(drift_post_cor)) = 0;
            all_post_cor = round(fix_post_cor+drift_post_cor) + max_Tshift + 1;
            
            %RECOMPUTE XMAT
            for i=1:NT
                all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
            end
            X{1} = create_time_embedding(all_shift_stimmat_up,stim_params_us);
            
            if add_usfac > 1
                X{1} = tb_proc_stim(X{1}(used_inds,use_kInds_up),add_usfac,flen);
            else
                X{1} = X{1}(used_inds,use_kInds_up);
            end
            %% REFIT XV CELLS
            silent = 1;
            for ss = 1:length(tr_set)
                cur_cell = tr_set(ss);
                fprintf('Drift LOO %d/%d, Refitting model %d of %d\n',xv,length(loo_set),ss,length(tr_set));
                cur_unit_ind = find(tr_set == cur_cell);
                cur_fit_inds = all_modfit_inds(~isnan(Robs_mat(all_modfit_inds,cur_unit_ind)));
                
                %                 if tr_set(ss) == tot_nUnits && model_pop_avg
                %                     X{2} = Xblock(used_inds,:); %for modeling pop-avg dont use the pop-avg as a predictor
                %                 end
                
                dit_mods_LOO{xv,nn+1}(cur_cell) = dit_mods{nn+1}(cur_cell);
                dit_mods_LOO{xv,nn+1}(cur_cell) = NMMfit_filters(dit_mods_LOO{xv,nn+1}(cur_cell),Robs_mat(:,cur_unit_ind),...
                    X,[],cur_fit_inds,silent); %fit stimulus filters
                
                dit_mods_spkNL_LOO{xv,nn+1}(cur_cell) = NMMfit_logexp_spkNL(dit_mods_LOO{xv,nn+1}(cur_cell),Robs_mat(:,cur_unit_ind),X,[],cur_fit_inds);
                
                newLL = NMMeval_model(dit_mods_spkNL_LOO{xv,nn+1}(cur_cell),Robs_mat(:,cur_unit_ind),X,[],cur_fit_inds);
                dit_LLimp_LOO(xv,nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
                %                 fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(cur_cell),dit_LLimp(nn,cur_cell),dit_LLimp(nn+1,cur_cell));
            end
            
            %             if model_pop_avg
            %                X{2} = [Xblock(used_inds,:) pop_rate];
            %             end
        end
    end
end

%% SAVE EYE-TRACKING RESULTS
et_params = struct('beg_buffer',params.beg_buffer,'end_buffer',params.end_buffer,'min_trial_dur',params.min_trial_dur,'bar_ori',bar_ori,'good_coils',params.good_coils,...
    'use_nPix',use_nPix,'flen',flen,'dt',dt,'drift_jump_sigma',drift_jump_sigma,'drift_prior_sigma',drift_prior_sigma,...
    'fix_prior_sigma',fix_prior_sigma,'fix_noise_sigma',fix_noise_sigma,'drift_noise_sigma',drift_noise_sigma,...
    'drift_dsf',drift_dsf,'n_fix_inf_it',n_fix_inf_it,'n_drift_inf_it',n_drift_inf_it,'use_sac_kerns',use_sac_kerns,'shifts',shifts,...
    'use_measured_pos',use_measured_pos,'sac_bincents',sac_bincents,'spatial_usfac',spatial_usfac,'add_usfac',add_usfac,...
    'sac_shift',sac_shift,'use_coils',params.use_coils,'sp_dx',sp_dx,'model_pop_avg',model_pop_avg,'poss_lambda_scales',poss_lambda_scales);

et_rand_fixpos = rand_fixpos;
et_used_inds = used_inds;
et_tr_set = tr_set;
et_tr_trials = tr_trials;
et_xv_trials = xv_trials;
et_saccades = ET_data.saccades;
et_is_blink = ET_data.is_blink;
et_clust_data = Clust_data;
% cd(anal_dir);
% save(anal_name,'it_*','drift_post_*','fix_ids','dit_*','et_used_inds','et_tr_set','et_clust_data','et_saccades','et_is_blink','et_params','et_tr_trials','et_xv_trials','et_rand_fixpos');

%%
% fin_fix_corr = nan(NT,1);
% fin_fix_std = nan(NT,1);
% fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
% fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
% fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
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
%

%%
% % close all
% n_trials = length(unique(all_trialvec));
% for tt = n_trials:-1:1
%     fprintf('Trial %d/%d, Se: %d\n',tt,n_trials,trial_data(tt).se);
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
%
%%
% close all
% f1 = figure();
% f2 = figure();
% f3 = figure();
% f4 = figure();
% f5 = figure();
% for ss = 1:length(tr_set)
%     % sbeg = find(all_mod_SU(tr_set) > 0,1);
%     % for ss = sbeg:length(tr_set)
%     ss
%     init_mod = it_mods{1}(tr_set(ss));
%     xtargs = [init_mod.mods(:).Xtarget];
%     kmat = [init_mod.mods(xtargs == 1).filtK];
%     figure(f1); clf
%     subplot(2,2,1)
%     imagesc(reshape(kmat(:,1),flen,use_nPix_us/add_usfac));
%     ca = max(abs(kmat(:,1))); caxis([-ca ca]);
%     for ii = 1:(size(kmat,2)-1)
%         subplot(2,2,2+ii)
%         imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us/add_usfac));
%         ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
%     end
%     colormap(gray)
%
% %     fin_mod = it_mods{2}(tr_set(ss));
% %     xtargs = [fin_mod.mods(:).Xtarget];
% %     kmat = [fin_mod.mods(xtargs == 1).filtK];
% %     figure(f2); clf
% %     subplot(2,2,1)
% %     imagesc(reshape(kmat(:,1),flen,use_nPix_us/add_usfac));
% %     ca = max(abs(kmat(:,1))); caxis([-ca ca]);
% %     for ii = 1:(size(kmat,2)-1)
% %         subplot(2,2,2+ii)
% %         imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us/add_usfac));
% %         ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
% %     end
% %     colormap(gray)
%
%     fin_mod = it_mods{end}(tr_set(ss));
%     xtargs = [fin_mod.mods(:).Xtarget];
%     kmat = [fin_mod.mods(xtargs == 1).filtK];
%     figure(f3); clf
%     subplot(2,2,1)
%     imagesc(reshape(kmat(:,1),flen,use_nPix_us/add_usfac));
%     ca = max(abs(kmat(:,1))); caxis([-ca ca]);
%     for ii = 1:(size(kmat,2)-1)
%         subplot(2,2,2+ii)
%         imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us/add_usfac));
%         ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
%     end
%     colormap(gray)
%
%     fin_mod = dit_mods{end}(tr_set(ss));
%     xtargs = [fin_mod.mods(:).Xtarget];
%     kmat = [fin_mod.mods(xtargs == 1).filtK];
%     figure(f4); clf
%     subplot(2,2,1)
%     imagesc(reshape(kmat(:,1),flen,use_nPix_us/add_usfac));
%     ca = max(abs(kmat(:,1))); caxis([-ca ca]);
%     for ii = 1:(size(kmat,2)-1)
%         subplot(2,2,2+ii)
%         imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us/add_usfac));
%         ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
%     end
%     colormap(gray)
%
% %     figure(f5); clf
% %     plot(1:n_fix_inf_it+1,it_LLimp(:,tr_set(ss)),'o-');
% %     hold on
% %     plot((n_fix_inf_it+1):(n_fix_inf_it+n_drift_inf_it+1),dit_LLimp(:,tr_set(ss)),'ro-');
%
%     fprintf('Cell %d of %d\n',ss,length(tr_set));
%     fprintf('Original: %.4f Next-fix: %.4f  Fin-fix: %.4f Fin: %.4f\n',it_LLimp(1,tr_set(ss)),it_LLimp(2,tr_set(ss)),...
%         it_LLimp(end,tr_set(ss)),dit_LLimp(end,tr_set(ss)));
%     pause
% end

