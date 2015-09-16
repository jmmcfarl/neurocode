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

init_mods_name = 'NIM_ET_initmods';
ET_results_name = 'NIM_ET';

init_mods_name = [init_mods_name sprintf('_ori%d',bar_ori)];
ET_results_name = [ET_results_name sprintf('_ori%d',bar_ori)];
if rec_number > 1
    init_mods_name = strcat(init_mods_name,sprintf('r%d',rec_number));
    ET_results_name = strcat(ET_results_name,sprintf('r%d',rec_number));
end

%%
min_trial_dur = 0.75; %minimal trial dur (s)
%exclude data at beginning and end of each trial (numbers in sec)
beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

stim_fs = 100; %in Hz
dt = 0.01; %in sec
Fr = 1;

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
full_nPix = params.full_nPix;
if use_nPix > full_nPix
    fprintf('Using npix == full_nPix\n');
    use_nPix = full_nPix;
end
full_nPix = full_nPix; 
use_nPix = use_nPix;

%for bar widths bigger than this do additional spatial upsampling (using
%interpolation, no filter up-sampling)
if mode(expt_data.expt_dw)/params.scale_fac > 1.5*.0565
    add_usfac = 2;
else
    add_usfac = 1;
end
if ~isempty(params.rpt_seeds)
    has_rpts = true;
else
    has_rpts = false;
end

base_spatial_usfac = 2; %baseline spatial up-sampling (up-sampling of filters)
spatial_usfac = add_usfac*base_spatial_usfac; %total spatial up-sampling
full_nPix_us = spatial_usfac*full_nPix;
use_nPix_us = use_nPix*spatial_usfac;
flen = 12; %number of stimulus time lags to use
klen_us = use_nPix_us*flen;

%%
%model-fiting parameters
MP.xv_frac = 0.2; %fraction of trials to use for cross-validation
MP.recompute_init_mods = 0; %use existing initial models?
MP.use_sac_kerns = 1; %use sac-modulation kernels
MP.model_pop_avg = 1; %include an explicit model of population avg (and use it as a predictor in unit models)
MP.pop_avg_sigma = 0.05; %gaussian smoothing sigma for pop avg rates
MP.init_jitter_SD = 0.075; %SD for initial random EP distribution

%saccade kernel time axis
sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_lag_set = -sac_backlag:sac_forlag;
n_sac_lags = length(sac_lag_set);

base_sp_dx = mode(expt_data.expt_dw); %size of bars in pixels
if length(unique(expt_data.expt_dw)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac/params.scale_fac; %model dx in deg

%%
HMM_params.max_fix_shift = 0.85; %maximum shift to consider when estimating fixation-corrections (in deg)
HMM_params.max_drift_shift = 0.45; %max shift for drift corrections (in deg)
HMM_params.neural_delay = 0.05; %how much to shift jump points forward in time
HMM_params.n_fix_iter = 3; %3 number of iterations to estimate fixation-based EP
HMM_params.n_drift_iter = 1; %1 number of iterations to estimate within-fixation drift
% set priors on EP (all values in deg)
HMM_params.fix_prior_sigma = 0.15; %prior sigma on fixation-based EP
HMM_params.fix_delta_noise_sigma = 0.1; %sigma on noise of eye-trackers in measuring change in EP between fixations
%this sets the posterior sigma to be ~0.004 in all cases (just a heuristic)
HMM_params.drift_noise_sigma = 0; %gaussian noise sigma on coil-measured drift
HMM_params.drift_prior_sigma = 0.004; %sigma of gaussian prior on diff of eye-pos
if all(params.use_coils) %if using both coils
    HMM_params.drift_prior_sigma = sqrt(0.004^2*3);
    HMM_params.drift_noise_sigma = sqrt(0.004^2*3);
elseif any(params.use_coils) %if just using one coil
    HMM_params.drift_prior_sigma = sqrt(0.004^2*2);
    HMM_params.drift_noise_sigma = sqrt(0.004^2*2);
end
HMM_params.drift_jump_sigma = 0.075; %gaussian prior sigma on change in EP between fixations during drift-corrections
HMM_params.drift_dsf = 3; %temporal down-sampling for estimating drift

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
if MP.use_sac_kerns %could do separate predictors for big sacs and micros sacs but they have similar effects on firing rates, so easiest to pool them
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
n_fixs = length(fix_start_inds);

%is the given fixation following a blink or at the beginning of a new trial?
fix_post_break = ismember(fix_start_inds,[saccade_stop_inds(used_is_blink); trial_end_inds+1]); 

%compute fixation id vector
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

fix_boundaries = [fix_start_inds fix_stop_inds]; %store [Nx2] mat of fixation boundaries
trial_boundaries = [trial_start_inds trial_end_inds]; %store [Nx2] mat of trial boundary indices

%% INITIALIZE STIMULUS FOR INITIAL MODEL FITTING
rand_fixpos = randn(n_fixs,1)*MP.init_jitter_SD; %draw a random position for each fixation
init_eyepos = zeros(NT,1);
all_stimmat_shift = all_stimmat_up;
if MP.init_jitter_SD > 0
    init_eyepos(used_inds(~isnan(fix_ids))) = rand_fixpos(fix_ids(~isnan(fix_ids)));
    max_sim_pos = full_nPix*sp_dx/2;%maximum initial corrections (deg)
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
n_xv_trials = round(MP.xv_frac*nuse_trials);
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

if MP.model_pop_avg
    norm_rates = Robs_mat(:,1:n_probes); %define pop rate over MUA only
    for tt = 1:length(trial_data) %loop over trials
        cur_trial_inds = find(all_trialvec(used_inds) == tt);
        if ~isempty(cur_trial_inds)
            for ss = 1:n_probes %smooth the spk cnt data with a Gaussian kernel
                norm_rates(cur_trial_inds,ss) = jmm_smooth_1d_cor(norm_rates(cur_trial_inds,ss),...
                    round(MP.pop_avg_sigma/dt));
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
mod_stim_params(1) = NIM.create_stim_params([flen use_nPix_us/add_usfac],...
    'stim_dt',dt,'stim_dx',sp_dx*add_usfac);
if MP.model_pop_avg
    X{2} = [Xblock pop_rate]; %if using a pop-rate predictor, add it to the block-index predictors
else
    X{2} = [Xblock];
end
mod_stim_params(2) = NIM.create_stim_params([n_block_filts 1]);
if MP.use_sac_kerns
    X{3} = Xsac;
end
mod_stim_params(3) = NIM.create_stim_params(n_sac_lags,'boundary_conds',[0 0 0]);

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
silent = 1;
MP.sac_d2t_lambda = 100; %temporal smoothness regularization on saccade kernels
MP.block_L2_lambda = 1; %slight penalty on block coefs to ensure convergence

%baseline regularization strengths
MP.base_lambda_d2XT = 100;
MP.base_lambda_L1 = 5;
MP.init_lambda_d2XT = 100;
%if using sparser stimuli, reduce base reg strenghts to approximately
%account for difference in stimulus marginal variance
if unique(expt_data.expt_dds) == 12
   MP.base_lambda_d2XT = MP.base_lambda_d2XT/5;
   MP.base_lambda_L1 = MP.base_lambda_L1/2;
   MP.init_lambda_d2XT = MP.init_lambda_d2XT/5;
end

%stimulus subunits
MP.n_Esquared_filts = 2; %number of excitatory squared filters
MP.n_Isquared_filts = 0; %number of inhibitory squared filters
mod_signs = [1 ones(1,MP.n_Esquared_filts) -1*ones(1,MP.n_Isquared_filts)];
NL_types = [{'lin'} repmat({'quad'},1,MP.n_Esquared_filts+MP.n_Isquared_filts)]; %NL types for full model
Xtargs = ones(1, MP.n_Esquared_filts + MP.n_Isquared_filts + 1);

%add block-predictor subunit
mod_signs = [mod_signs 1];
NL_types = [NL_types {'lin'}];
Xtargs = [Xtargs 2];

%initial regularization parameters (only smoothness reg on filters at
%first, used to estimate filter scales and adjust relative reg strengths)
init_d2XT = [MP.init_lambda_d2XT*ones(MP.n_Esquared_filts + MP.n_Isquared_filts + 1,1); 0;];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT);

%range of smoothness regularization 'rescaling' to try
% MP.poss_lambda_scales = logspace(-2,2,10);
MP.poss_lambda_scales = 1;

if ~exist(['./' init_mods_name '.mat'],'file') || MP.recompute_init_mods == 1
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
        if MP.use_sac_kerns
            null_mod = NIM(mod_stim_params,repmat({'lin'},1,2),[1 1],'Xtargets',[2 3],...
                'd2t',[0 MP.sac_d2t_lambda],'l2',[MP.block_L2_lambda 0]);
        else
            null_mod = NIM(mod_stim_params,'lin',1,'Xtargets',2); %otherwise null model just has block predictors
        end
        null_mod = null_mod.fit_filters(Robs,X,cur_tr_inds,'silent',silent);
        null_xvLL(ss) = null_mod.eval_model(Robs,X,cur_xv_inds);
        
        %% estimate initial models and adjust base reg strength
        init_gqm = NIM(mod_stim_params,NL_types,mod_signs,'Xtargets',Xtargs,...
            'd2xt',MP.init_lambda_d2XT*(Xtargs==1),'l2',MP.block_L2_lambda*(Xtargs == 2));
        init_gqm = init_gqm.fit_filters(Robs,X,all_modfit_inds,'silent',silent);
        
        %adjust regularization stregnths of each stimulus filter relative to the scale of filter outputs
        [~,~,mod_internals] = init_gqm.eval_model(Robs,X,all_modfit_inds);
        new_lambda_d2xt = MP.base_lambda_d2XT./var(mod_internals.gint(:,Xtargs == 1));
        new_lambda_l1 = MP.base_lambda_L1./std(mod_internals.gint(:,Xtargs == 1));
        init_gqm = init_gqm.set_reg_params('sub_inds',find(Xtargs == 1),'d2xt',new_lambda_d2xt,'l1',new_lambda_l1);
        
        %% scan a range of smoothness regularization strengths and pick the best using xval
        cur_lambda_d2XT = init_gqm.get_reg_lambdas('sub_inds',find(Xtargs == 1),'d2xt');
        new_gqm = init_gqm;
        new_gqm_xvLL = nan(length(MP.poss_lambda_scales),1);
        clear all_gqms
        for ll = 1:length(MP.poss_lambda_scales) %loop over possible range of smoothness lambdas
            fprintf('%d ',ll);
            new_lambda_d2XT = cur_lambda_d2XT*MP.poss_lambda_scales(ll);
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
        if MP.use_sac_kerns
            up_gqm = up_gqm.add_subunits('lin',[1],'xtargs',[3],'d2t',MP.sac_d2t_lambda);
            up_gqm = up_gqm.fit_filters(Robs,X,cur_tr_inds,'silent',silent);
        end
        xvLL = up_gqm.eval_model(Robs,X,cur_xv_inds);
        all_mod_xvLLimp(ss) = (xvLL - null_xvLL(ss))/log(2);
        
        %% refit models using all data
        up_gqm = up_gqm.fit_filters(Robs,X,all_modfit_inds,'silent',silent);
        null_mod = null_mod.fit_filters(Robs,X,all_modfit_inds,'silent',silent);
        
        %% store results
        all_nullmod(ss) = null_mod;
        all_mod_fits(ss) = up_gqm.fit_spkNL(Robs,X,all_modfit_inds,'silent',silent); %fit spkNL params
        null_LL(ss) = null_mod.eval_model(Robs,X,all_modfit_inds);
        LL = all_mod_fits(ss).eval_model(Robs,X,all_modfit_inds);
        all_mod_LLimp(ss) = (LL - null_LL(ss))/log(2); %LL improvement in bits per spike
        all_mod_fits(ss).fit_props.null_LL = null_LL(ss); %set the null using the complete null model
    end
end
    save(init_mods_name,'all_mod_*','all_nullmod','Clust_data','null_*LL','*_trials','rand_fixpos','MP');
else
    fprintf('Loading pre-computed initial models\n');
    load(init_mods_name);
end

%%
%if using jittered initialization, reconstruct unjittered stimulus
if MP.init_jitter_SD > 0
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
        cur_pinds = cur_inds + round(HMM_params.neural_delay/dt);
        if length(cur_inds) > round(HMM_params.neural_delay/dt)
            forward_measured_drift(cur_inds(round(HMM_params.neural_delay/dt)+1:end),:) = ...
                measured_drift(cur_inds(1:(end-round(HMM_params.neural_delay/dt))),:);
        end
    end
    measured_drift = [zeros(1,2); diff(forward_measured_drift)];
    clear forward_measured_drift
    
    measured_fix_deltas = nan(n_fixs,1); %measured change in eye position between adjacent fixations
    if all(params.use_coils == 1) %if using both coils
        measured_fix_deltas(2:end) = mean(diff(measured_fix_med),2); %avg change in eye position between fixation
        measured_drift = mean(measured_drift,2); %avg of measured drifts
        n_coil_samps = 2;
    else %if only using one coil
        measured_fix_deltas(2:end) = diff(measured_fix_med(:,params.use_coils == 1)); %just use the good coil
        measured_drift = measured_drift(:,params.use_coils == 1);
        n_coil_samps = 1;
    end
else
    measured_fix_deltas = [];
    measured_drift = [];
    n_coil_samps = 0;
end

ET_meas = struct('fix_deltas',measured_fix_deltas,'drift',measured_drift,'n_used_coils',n_coil_samps);

%% package data structure to pass into the EM algo
gen_data = struct('stim_dx',sp_dx,'add_usfac',add_usfac,'trial_boundaries',trial_boundaries,...
    'fix_boundaries',fix_boundaries,'used_inds',used_inds,'stim_params',stim_params,'dt',dt,...
    'use_kInds',use_kInds,'fix_post_break',fix_post_break,'modfit_inds',all_modfit_inds);

%% run EM algo using all units to infer EP
[EP_pos,mod_fits] = eyetracking_EM(...
    all_mod_fits(tr_units),Robs_mat(:,tr_units),all_stimmat_up,X(2:end),...
    tr_units,HMM_params,ET_meas,gen_data);

fix_post_mean = squeeze(EP_pos.fix_mean(end,:));
fix_post_std = squeeze(EP_pos.fix_std(end,:));
drift_post_mean = squeeze(EP_pos.drift_mean(end,:));
drift_post_std = squeeze(EP_pos.drift_std(end,:));
[fin_tot_corr,fin_tot_std] = construct_eye_position(fix_post_mean,fix_post_std,...
    drift_post_mean,drift_post_std,fix_ids,trial_start_inds,trial_end_inds,round(HMM_params.neural_delay/dt));

%% Run the algo on all sets of units (one different set for each unit we want to compute LOO estimates for)
if use_LOOXV %if doing LOO on each SU
    loo_set = find(all_mod_SU(tr_units) > 0); %subset of tr_units that are SUs
else
    loo_set = [];
end
for xv = 1:length(loo_set)
    fprintf('Running ET EM on xv set %d of %d\n',xv,length(loo_set));
    use_set = setdiff(1:length(tr_units),loo_set(xv));
    [EP_pos_xv{xv},mod_fits_xv{xv}] = eyetracking_EM(...
        all_mod_fits(tr_units),Robs_mat(:,tr_units),all_stimmat_up,X(2:end),...
        use_set,HMM_params,ET_meas,gen_data);
    fix_post_mean = squeeze(EP_pos_xv{xv}.fix_mean(end,:));
    fix_post_std = squeeze(EP_pos_xv{xv}.fix_std(end,:));
    drift_post_mean = squeeze(EP_pos_xv{xv}.drift_mean(end,:));
    drift_post_std = squeeze(EP_pos_xv{xv}.drift_std(end,:));
    [fin_tot_corr_xv(xv,:),fin_tot_std_xv(xv,:)] = construct_eye_position(fix_post_mean,fix_post_std,...
        drift_post_mean,drift_post_std,fix_ids,trial_start_inds,trial_end_inds,round(HMM_params.neural_delay/dt));
    
end



