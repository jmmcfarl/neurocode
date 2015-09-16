function [EP_post,mod_fits] = eyetracking_EM(...
    init_mods,Robs_mat,stim_mat,nonstim_X,tr_units,HMM_params,ET_meas,gen_data)
% [EP_pos,mod_fits] = eyetracking_EM(init_mods,Robs_mat,stim_mat,nonstim_X,tr_units,HMM_params,ET_meas,gen_data)
% Runs HMM-based EM algorithm iteratively on neural data to infer eye position.
% INPUTS:   
%         init_mods: [Nx1] array of initial model fits across all N units. Must be NIM objects
%         Robs_mat: [T x N] matrix of observed spiking data across all T time points used in analysis. 
%             nan values indicate times when a unit was not isolated
%         stim_mat: [T_full x K] matrix of pixel values across all T_full time points 
%             (note: T_full >= T because we might exclude some time points from analysis, but we still
%             want to have the stimulus stored at those points to reconstruct the time-embedded stimulus at the boundaries
%         nonstim_X: cell array of predictors ([T x k] matrices) other than the stimulus.
%             Assumes the models have a linear dependence on each of these additional predictors
%         tr_units: subset of units (out of the full set 1:N) used for inferring eye position
%         HMM_params: struct of parameter values for the HMM inference
%             max_fix_shift: maximum allowable inferred fixation-avg eye position (deg)
%             max_drift_shift: maximum allowable drift-correction (on top of the fixation corrections). Also in deg
%             neural_delay: fixed time delay (in sec) associated with neural response in V1 relative to stimulus on retina
%             n_fix_iter: number of EM iterations for inferring fixation-based corrections
%             n_drift_iter: number of EM iterations for inferring additional drift-corrections
%             fix_prior_sigma: SD of the gaussian prior on fixation position (deg)
%             fix_delta_noise_sigma: SD of the gaussian noise in coil-measured between-fixation changes in eye-position (deg)
%             drift_noise_sigma: SD of gaussian noise in coil-measured drift (change in position between adjacent time points). (deg)
%             drift_prior_sigma: SD of gaussian prior on drift magnitude (change in position between adjacent time points) (deg)
%             drift_jump_sigma: SD of gaussian prior on the initial drift correction at the beginning of each fixation (deg)
%             drift_dsf: temporal down-sampling factor of inferring drift corrections within each fixation
%         ET_meas: struct of measured eye position signals
%             fix_delta: [n_fixation x 1] vector of measured between-fixation changes in eye position (deg)
%             drift: [T x 1] vector of measured drifts (change in EP between adjacent time points) (deg)
%             n_used_coils: scalar indicating the number of (presumed independent) measurements used to create measured EP signals (e.g. number of coils used)
%         gen_data: struct of additional data 
%             stim_dx: pixel resolution of stimulus (deg)
%             add_usfac: spatial up-sampling of model filters (done using linear interpolation with tent-bases)
%             trial_boundaries: [n_trial x 2] matrix of start and stop indices of each trial (relative to used_inds)
%             fix_boundaries: [n_fixation x 2] matrix of start and stop indices of each fixation (relative to used_inds)
%             used_inds: Tx1 vector of the time point indices (relative to the full T_full) used in analysis
%             stim_params: stim_params struct, created by NIM.create_stim_params, used to construct the stim_mat
%             dt: time resolution (sec)
%             use_kInds: vector of indices relative to the full time-embedded stimulus matrix on which our model filters
%                 act. Usually this specifies the range of central pixels where the units have RFs\
%             fix_post_break: [n_fix x 1] bool indicating whether each fixation should ignore
%                 coil-measured information on deltaX (e.g. whether it follows a blink, or a new trial)
%             modfit_inds: vector of indices on which the model parameters should be fit (e.g. excludes repeat trials)
%               
% OUTPUTS:
%         EP_post: struct of results:
%             EP_post.fix_mean: [n_fix_iters x n_fixations] array of the posterior mean EPs at each fixation
%             EP_post.fix_std: same as above, with the posterior SDs in each fixation
%             EP_post.drift_mean: [n_drift_iters x T] array of pontwise posterior mean inferred drift
%             EP_post.drift_std: same as above, with the posterior SDs of inferred drift
%         mod_fits: struct of model fits across all units at each EM iteration
%             mod_fits.fix_iters: n_fix_iterx1 cell array, each element containing an Nx1 array of model fits at that iteration
%             mod_fits.drift_iters: same for drift-iterations

%unpack some useful variables from structs
add_usfac = gen_data.add_usfac;
%assuming all models have the same structure!
mod_stim_params = init_mods(1).stim_params;
drift_dsf = HMM_params.drift_dsf;
stim_dx = gen_data.stim_dx;
neural_delay = (HMM_params.neural_delay/gen_data.dt); %in units of time bins

n_tr_units = length(tr_units); %number of units used to infer eye-pos
[NT,n_units] = size(Robs_mat);
nPix = size(stim_mat,2); %total number of spatial pixels in the stimulus
n_add_covariates = length(nonstim_X); %number of additional "Xmats" (beyond the visual stim)
n_trials = length(gen_data.trial_boundaries);
use_coils = ET_meas.n_used_coils > 0;

stim_klen = prod(mod_stim_params(1).dims); %number of coefs in stimulus filters. visual stimulus must be the first Xtarg
flen = mod_stim_params(1).dims(1); %number of time lags
Xtargs = [init_mods(1).subunits(:).Xtarg]; %vector of Xtargets for model subunits (all models must have same structure)
stim_mod_signs = [init_mods(1).subunits(find(Xtargs == 1)).weight]; %subunit weightss
stim_NL_types = init_mods(1).get_NLtypes(find(Xtargs == 1)); %subunit NL types
n_stim_subs = sum(Xtargs == 1); %number of stimulus-targeting subunits
assert(strcmp(init_mods(1).spkNL.type,'softplus'),'only works with softplus spkNL at this point');

%% compute fixation indices (normal and forward-projected)
fix_boundaries = gen_data.fix_boundaries;
trial_start_inds = gen_data.trial_boundaries(:,1);
n_fixs = size(fix_boundaries,1);
pfix_boundaries = fix_boundaries;
for ii = 1:n_fixs
    next_trial = trial_start_inds(find(trial_start_inds >= fix_boundaries(ii,1),1,'first'));
    if next_trial > fix_boundaries(ii,1) + neural_delay %if the forward projection does not push you into another trial
        pfix_boundaries(ii,1) = fix_boundaries(ii,1) + neural_delay;
    end
    next_trial = trial_start_inds(find(trial_start_inds >= fix_boundaries(ii,2),1,'first'));
    if next_trial > fix_boundaries(ii,2) + neural_delay
        pfix_boundaries(ii,2) = fix_boundaries(ii,2) + neural_delay;
    end
end

%vectors of index values for each fixation (normal and forward-projected)
fix_ids = nan(NT,1);
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_boundaries(ii,1):fix_boundaries(ii,2);
    fix_ids(cur_inds) = ii;
    cur_inds = pfix_boundaries(ii,1):pfix_boundaries(ii,2);
    pfix_ids(cur_inds) = ii;
end

%% generate shift matrices. Must be applied to the stimulus (not the filters)
%shifts for inferring fixation corrections
max_fix_shift_pix = round(HMM_params.max_fix_shift/stim_dx); %max eye pos deviation (deg)
fix_shifts = -max_fix_shift_pix:max_fix_shift_pix; %shift-range
n_fix_shifts = length(fix_shifts); %number of shifts

%shifts for inferring drift
max_drift_shift_pix = round(HMM_params.max_drift_shift/stim_dx);
drift_shifts = -max_drift_shift_pix:max_drift_shift_pix;
n_drift_shifts = length(drift_shifts);

It = speye(flen); 

fix_shift_mats = cell(n_fix_shifts,1);
for xx = 1:n_fix_shifts
    temp = spdiags( ones(nPix,1), -fix_shifts(xx),nPix,nPix);
    temp = kron(temp,It);
    fix_shift_mats{xx} = temp(:,gen_data.use_kInds);
end

drift_shift_mats = cell(n_drift_shifts,1);
for xx = 1:n_drift_shifts
    temp = spdiags( ones(nPix,1), -drift_shifts(xx), nPix, nPix);
    temp = kron(temp,It);
    drift_shift_mats{xx} = temp(:,gen_data.use_kInds);
end

%% precompute HMM priors
%define updated posteriors (incorporating measured coil data)
if use_coils
    HMM_params.fix_delta_noise_sigma = HMM_params.fix_delta_noise_sigma/sqrt(ET_meas.n_used_coils); %if we are avging multiple measures of delta_EP, scale noise variance of the measure (assuming independent measurements)
    post_drift_var = 1/(ET_meas.n_used_coils/HMM_params.drift_noise_sigma^2 + 1/HMM_params.drift_prior_sigma^2); %variance of posterior gaussian, given N_coil_samps observations of drift, and gaussian prior
    post_drift_mean = ET_meas.n_used_coils*ET_meas.drift*post_drift_var/HMM_params.drift_noise_sigma^2; %mean of the gaussian posterior
    post_drift_sigma = sqrt(post_drift_var);
end

%log prior dist on fixations
fix_Lprior = -(fix_shifts*stim_dx).^2./(2*HMM_params.fix_prior_sigma^2); %log-prior on within-fixation eye-pos
fix_Lprior = fix_Lprior - logsumexp(fix_Lprior); %normalize

%log prior dist on drift-corrections
drift_jump_Lprior = -(drift_shifts*stim_dx).^2./(2*HMM_params.drift_jump_sigma^2); %log-prior on initial position (after a jump) during drift-inference
drift_jump_Lprior = drift_jump_Lprior - logsumexp(drift_jump_Lprior);

%log prior matrix for drift (not including coil-signals)
cdist = squareform(pdist(drift_shifts'*stim_dx));
base_LA = -cdist.^2/(2*(HMM_params.drift_prior_sigma*drift_dsf)^2); %baseline log-prior matrix on drift position changes
base_LA = bsxfun(@minus,base_LA,logsumexp(base_LA,2)); %normalize

%% construct initial stimulus Xmat
base_Xmat = NIM.create_time_embedding(stim_mat,gen_data.stim_params);
base_Xmat = base_Xmat(gen_data.used_inds,:);

%% ITERATE FIXATION-BASED CORRECTIONS

fix_it_mods{1} = init_mods;
fix_it_LLimp = nan(HMM_params.n_fix_iter+1,n_units); %initialize matrix of LL improvements
for ss = 1:length(init_mods) %store initial model LL improvemens
    fix_it_LLimp(1,ss) = (init_mods(ss).fit_props.LL - init_mods(ss).fit_props.null_LL)/log(2);
end

fix_it_post_mean = nan(HMM_params.n_fix_iter,n_fixs);
fix_it_post_std = nan(HMM_params.n_fix_iter,n_fixs);
for nn = 1:HMM_params.n_fix_iter %loop over EM iterations
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,HMM_params.n_fix_iter);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(stim_klen/add_usfac,sum(Xtargs==1),n_tr_units); %compile all stimulus filters
    lin_kerns = cell(n_add_covariates,1);
    for ii = 1:n_add_covariates %compile linear filters for all additional model covariates
        n_preds = prod(mod_stim_params(ii+1).dims);
        lin_kerns{ii} = nan(n_preds,n_tr_units);
    end
    mod_spkNL_params = nan(3,n_tr_units); %assuming softplus spkNL (with three params including offset)
    for ss = 1:n_tr_units
        filt_bank(:,:,ss) = cell2mat(fix_it_mods{nn}(tr_units(ss)).get_filtKs(find(Xtargs == 1))');
        mod_spkNL_params(:,ss) = [fix_it_mods{nn}(tr_units(ss)).spkNL.params fix_it_mods{nn}(tr_units(ss)).spkNL.theta];
        for ii = 1:n_add_covariates
           lin_kerns{ii}(:,ss) = cell2mat(fix_it_mods{nn}(tr_units(ss)).get_filtKs(find(Xtargs == ii + 1))); 
        end
    end
    lin_pred_out = zeros(NT,n_tr_units);
    for ii = 1:n_add_covariates %precompute output of additional linear terms
        lin_pred_out = lin_pred_out + nonstim_X{ii}*lin_kerns{ii};
    end
    
    %% ESTIMATE LL for each shift in each stimulus frame
    %precompute LL at all shifts summed over units
    frame_LLs = nan(NT,n_fix_shifts);
    reverseStr = '';
    for xx = 1:n_fix_shifts
        msg = sprintf('Calculating LLs for fixation shift %d of %d\n',xx,n_fix_shifts);
        fprintf([reverseStr msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        
        cur_stim_shift = base_Xmat*fix_shift_mats{xx}; %shift stimulus matrix
        
        %process stimulus with spatial tent-basis if doing additional up-sampling
        if add_usfac > 1
            cur_stim_shift = tb_proc_stim(cur_stim_shift,add_usfac,flen);
        end
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,n_tr_units);
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(3,:)); %add in constant offset terms
        for ff = 1:n_stim_subs %add in contribution from all stimulus filters
            cur_filt_outs = cur_stim_shift*squeeze(filt_bank(:,ff,:));
            if ~strcmp(stim_NL_types{ff},'lin')
                cur_filt_outs = init_mods(1).subunits(ff).apply_NL(cur_filt_outs);
            end
            gfuns = gfuns + stim_mod_signs(ff)*cur_filt_outs;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + lin_pred_out;
        
        %incorporate beta parameter in sotplus function
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(1,:));
        
        %handle numerical overflow with softplus
        too_large = gfuns > 50;
        pred_rate = log(1+exp(gfuns));
        pred_rate(too_large) = gfuns(too_large);
        
        %incorporate alpha in softplus
        pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(2,:));
        
        %enforce min predicted rate to avoid nans
        pred_rate(pred_rate < 1e-50) = 1e-50;
        
        frame_LLs(:,xx) = squeeze(nansum(Robs_mat(:,tr_units).*log(pred_rate) - pred_rate,2)); %total LL at each time point for this shift
    end
    
    %% INFER MICRO-SAC SEQUENCE    fix_LLs = nan(n_fixs,n_fix_shifts);
    fix_LLs = nan(n_fixs,n_fix_shifts);
    for ii = 1:n_fixs %compute LL of each eye position summed across data within each fixation
        cur_inds = pfix_boundaries(ii,1):pfix_boundaries(ii,2); %use forward-projected fixation indices to account for neural delay
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:)); %sum LLs for data within the fixation
    end
    
    if ~use_coils %if not using coil info we just compute the posterior at each fixation independently
        Lgamma = bsxfun(@plus,fix_LLs,fix_Lprior); %add log-likelihood and log-prior
        Lgamma = bsxfun(@minus,Lgamma,logsumexp(Lgamma,2)); %normalize
    else %if we are using coil info on the between-fixation changes in EP, need an HMM
        Lalpha=zeros(n_fixs,n_fix_shifts); %log of forward messages
        Lbeta = zeros(n_fixs,n_fix_shifts); %log of backward messages
        Lalpha(1,:) = fix_Lprior + fix_LLs(1,:); %initialize using prior
        for t=2:n_fixs %compute forward messages
            if ~gen_data.fix_post_break(t) %if this fixation is NOT following a break
                %make quadratic transition matrix centered on coil-measured deltaX, with variance
                %given by the delta_noise_sigma^2
                cdist = pdist2(fix_shifts'*stim_dx + ET_meas.fix_deltas(t),fix_shifts'*stim_dx);
                cur_LA = -cdist.^2/(2*HMM_params.fix_delta_noise_sigma^2);
                cur_LA = bsxfun(@plus,cur_LA,fix_Lprior); %multiply this by the prior over fixation positions
                cur_LA = bsxfun(@minus,cur_LA,logsumexp(cur_LA,2)); %normalize log transition-matrix
                Lalpha(t,:) = logmulexp(Lalpha(t-1,:),cur_LA) + fix_LLs(t,:); %compute forward message
            else %if it is following a break, ignore the eye-trackers
                Lalpha(t,:) = fix_Lprior + fix_LLs(t,:);
            end
        end
        
        Lbeta(n_fixs,:)=zeros(1,n_fix_shifts); %initialize log-backward message
        for t=n_fixs-1:-1:1%compute backward messages
            if ~gen_data.fix_post_break(t+1) %if the fixation isn't followed by a break
                %construct A contingent on measured deltaX between current and NEXT fixation
                cdist = pdist2(fix_shifts'*stim_dx + ET_meas.fix_deltas(t+1),fix_shifts'*stim_dx);
                cur_LA = -cdist.^2/(2*HMM_params.fix_delta_noise_sigma^2);
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
    fix_it_post_mean(nn,:) = sum(bsxfun(@times,gamma,fix_shifts),2);
    cur_diff = bsxfun(@minus,fix_it_post_mean(nn,:)',fix_shifts).^2;
    fix_it_post_std(nn,:) = sqrt(sum(cur_diff.*gamma,2));
    
    %back-project saccade-times and interpolate over saccades
    fix_post_mean_cor = nan(NT,1);
    fix_post_mean_cor(~isnan(fix_ids)) = fix_it_post_mean(nn,fix_ids(~isnan(fix_ids)));
    fix_post_mean_cor = interp1(find(~isnan(fix_ids)),fix_post_mean_cor(~isnan(fix_ids)),1:NT); %interpolate position between fixations
    fix_post_mean_cor(isnan(fix_post_mean_cor)) = 0;
    
    %% adjust stimulus matrix for fixation eye positions
    cur_fix_shifts = round(fix_post_mean_cor); %round to pixels
    shift_stim_mat = stim_mat;
    for ii=1:NT %apply inferred fixation corrections
        shift_stim_mat(gen_data.used_inds(ii),:) = shift_matrix_Nd(stim_mat(gen_data.used_inds(ii),:),-cur_fix_shifts(ii),2);
    end
    all_Xmat = NIM.create_time_embedding(shift_stim_mat,gen_data.stim_params);
    if add_usfac %project onto spatial TBs if needed
        X{1} = tb_proc_stim(all_Xmat(gen_data.used_inds,gen_data.use_kInds),add_usfac,flen);
    else
        X{1} = all_Xmat(gen_data.used_inds,gen_data.use_kInds);
    end
    X(2:1+n_add_covariates) = nonstim_X;
    
    %% REFIT ALL CELLS
    silent = 1;
    for ss = 1:n_units
        cur_fit_inds = gen_data.modfit_inds(~isnan(Robs_mat(gen_data.modfit_inds,ss))); %all non-rpt indices when current unit was isolated
        
        fix_it_mods{nn+1}(ss) = fix_it_mods{nn}(ss); %initialize model from previous iteration
        fix_it_mods{nn+1}(ss).spkNL.params = [1 1]; %set alpha and beta to their default values when fitting filters to avoid issues with regularization
        fix_it_mods{nn+1}(ss) = fix_it_mods{nn+1}(ss).fit_filters(Robs_mat(:,ss),X,cur_fit_inds,'silent',silent,'sub_inds',[]); %fit constant offset first to get back our original model (sans other spkNL params)
        fix_it_mods{nn+1}(ss) = fix_it_mods{nn+1}(ss).fit_filters(Robs_mat(:,ss),X,cur_fit_inds,'silent',silent); %fit stimulus filters
        fix_it_mods{nn+1}(ss) = fix_it_mods{nn+1}(ss).fit_spkNL(Robs_mat(:,ss),X,cur_fit_inds,'silent',silent);  %refit spk NL
        newLL = fix_it_mods{nn+1}(ss).eval_model(Robs_mat(:,ss),X,cur_fit_inds);
        fix_it_LLimp(nn+1,ss) = (newLL - init_mods(ss).fit_props.null_LL)/log(2);
        if nn > 1
            fprintf('Unit %d LLimps.  Original LL: %.4f  Prev: %.4f  New: %.4f\n',ss,fix_it_LLimp(1,ss),fix_it_LLimp(nn,ss),fix_it_LLimp(nn+1,ss));
        else
            fprintf('Unit %d LLimps. Original: %.4f  New: %.4f\n',ss,fix_it_LLimp(1,ss),fix_it_LLimp(nn+1,ss));
        end
    end
end

%% parse stimulus matrix with best-estimate of fixation correctoins applied
if HMM_params.n_fix_iter == 0 %if were not doing any fixation corrections, need to create a new Xmat
    all_Xmat = NIM.create_time_embedding(stim_mat,gen_data.stim_params); 
end
all_Xmat = all_Xmat(gen_data.used_inds,:); %grab relevant time indices
base_Xmat = all_Xmat; %this is our new 'base' Xmat

%% now apply iteractive drift-based corrections
drift_it_mods{1} = fix_it_mods{end}; %initialize models from last fixation iteration

drift_it_LLimp = nan(HMM_params.n_drift_iter + 1,n_units);
drift_it_LLimp(1,:) = fix_it_LLimp(end,:);

drift_it_post_mean = nan(HMM_params.n_drift_iter,NT);
drift_it_post_std = nan(HMM_params.n_drift_iter,NT);
for nn = 1:HMM_params.n_drift_iter
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,HMM_params.n_drift_iter);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(stim_klen/add_usfac,sum(Xtargs==1),n_tr_units); %compile all stimulus filters
    lin_kerns = cell(n_add_covariates,1);
    for ii = 1:n_add_covariates
        n_preds = prod(mod_stim_params(ii+1).dims);
        lin_kerns{ii} = nan(n_preds,n_tr_units);
    end
    mod_spkNL_params = nan(3,n_tr_units); %assuming softplus spkNL (with three params including offset)
    for ss = 1:n_tr_units
        filt_bank(:,:,ss) = cell2mat(drift_it_mods{nn}(tr_units(ss)).get_filtKs(find(Xtargs == 1))');
        mod_spkNL_params(:,ss) = [drift_it_mods{nn}(tr_units(ss)).spkNL.params drift_it_mods{nn}(tr_units(ss)).spkNL.theta];
        for ii = 1:n_add_covariates
           lin_kerns{ii}(:,ss) = cell2mat(drift_it_mods{nn}(tr_units(ss)).get_filtKs(find(Xtargs == ii + 1))); 
        end
    end
    lin_pred_out = zeros(NT,n_tr_units);
    for ii = 1:n_add_covariates
        lin_pred_out = lin_pred_out + nonstim_X{ii}*lin_kerns{ii};
    end
    
    %% ESTIMATE LL for each shift in each stimulus frame
    %precompute LL at all shifts for all units
    frame_LLs = nan(NT,n_drift_shifts);
    reverseStr = '';
    for xx = 1:n_drift_shifts
        msg = sprintf('Calculating LLs for drift shift %d of %d\n',xx,n_drift_shifts);
        fprintf([reverseStr msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        cur_stim_shift = base_Xmat*drift_shift_mats{xx};
        
        %process with spatial TB if doing additional up-sampling
        if add_usfac > 1
            cur_stim_shift = tb_proc_stim(cur_stim_shift,add_usfac,flen);
        end
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,n_tr_units);
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(3,:)); %add in constant offset terms
        for ff = 1:n_stim_subs
            cur_filt_outs = cur_stim_shift*squeeze(filt_bank(:,ff,:));
            if ~strcmp(stim_NL_types{ff},'lin')
                cur_filt_outs = init_mods(1).subunits(ff).apply_NL(cur_filt_outs);
            end
            gfuns = gfuns + stim_mod_signs(ff)*cur_filt_outs;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + lin_pred_out;
        
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
        
        frame_LLs(:,xx) = squeeze(nansum(Robs_mat(:,tr_units).*log(pred_rate) - pred_rate,2));
    end
    
    %% INFER DRIFT CORRECTIONS
    Lgamma = nan(NT,n_drift_shifts);
    reverseStr = '';
    for ff = 1:n_fixs
        if mod(ff,100)==0
        msg = sprintf('Inferring drift in fixation %d of %d\n',ff,n_fixs);
        fprintf([reverseStr msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        end
        
        tset = find(pfix_ids==ff)'; %indices in this projected fixation
        ntset = length(tset);
        if ntset > drift_dsf %if we have at least drift_dsf points in this fixation
            nt_pts = ceil(ntset/drift_dsf); %number of points to infer drift at
            [talpha,tbeta] = deal(zeros(nt_pts,n_drift_shifts)); %initialize forward and backward messages
            
            tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf); %drift inference locations
            tpt_loc(end) = ntset;
            
            if use_coils; cur_drift_mean = post_drift_mean(tset); end;
            cur_LL_set = frame_LLs(tset,:);
            if mod(ntset,drift_dsf) ~= 0 %deal with any 'dangling points'
                dangling_pts = nt_pts*drift_dsf-ntset;
                cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_drift_shifts));
                if use_coils; cur_drift_mean = cat(1,cur_drift_mean,nan(dangling_pts,1)); end;
            end
            %sum LL over points within a downsampled bin
            cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_drift_shifts]);
            cur_LL_set = squeeze(sum(cur_LL_set,1));
            
            %precompute the transition matrix at each time bin within the fixation
            if ~use_coils %if we're not using coil-measured drift,
                cur_LA = repmat(base_LA,[1 1 nt_pts]); %log-trans matrix is constant
            else
                cur_drift_mean = drift_dsf*nanmean(reshape(cur_drift_mean,[drift_dsf nt_pts])); %assume constant drift vel, so measured deltaX is scaled by drift_dsf
                cur_LA = nan(n_drift_shifts,n_drift_shifts,nt_pts);
                for iii = 1:nt_pts
                    if ~isnan(cur_drift_mean(iii))
                        %transition matrix given by gaussian posterior with mean cur_drift_mean and variance
                        %determined by post_drift_sigma (scaled by drift_dsf)
                        cdist = pdist2(drift_shifts'*stim_dx + cur_drift_mean(iii),drift_shifts'*stim_dx);
                        cur_LA(:,:,iii) = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
                    else
                        cur_LA(:,:,iii) = base_LA; %otherwise set to prior
                    end
                end
                cur_LA = bsxfun(@minus,cur_LA,logsumexp(cur_LA,2)); %normalize
            end
            
            %forward messages
            talpha(1,:) = drift_jump_Lprior + cur_LL_set(1,:); %initialize forward messages within this fixation
            for t = 2:nt_pts
                talpha(t,:) = logmulexp(talpha(t-1,:),cur_LA(:,:,t)) + cur_LL_set(t,:);
            end
            
            %backward messages
            tbeta(end,:) = zeros(1,n_drift_shifts);
            for t = (nt_pts-1):-1:1
                lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
                tbeta(t,:) = logmulexp(lf1,cur_LA(:,:,t+1)');
            end
            temp_gamma = talpha + tbeta;
            
            %interpolate the posterior onto the original time axis if needed
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
    Lgamma = bsxfun(@minus,Lgamma,logsumexp(Lgamma,2)); %normalize log posterior
    gamma = exp(Lgamma); %posterior
    
    drift_it_post_mean(nn,:) = sum(bsxfun(@times,gamma,drift_shifts),2);
    drift_it_post_std(nn,:) = sqrt(sum(bsxfun(@times,gamma,drift_shifts.^2),2) - drift_it_post_mean(nn,:)'.^2);
    
    drift_it_post_mean(nn,isnan(drift_it_post_mean(nn,:))) = 0;
    drift_it_post_std(nn,isnan(drift_it_post_std(nn,:))) = 0;
    
    %interpolate over saccades
    drift_it_post_mean(nn,:) = interp1(find(~isnan(pfix_ids)),drift_it_post_mean(nn,~isnan(pfix_ids)),1:NT);
    drift_it_post_std(nn,:) = interp1(find(~isnan(pfix_ids)),drift_it_post_std(nn,~isnan(pfix_ids)),1:NT);
    drift_it_post_mean(nn,isnan(drift_it_post_mean(nn,:))) = 0;
    
    %% construct drift-corrected X-mat   
    drift_post_cor = squeeze(drift_it_post_mean(nn,:));
    for ii = 1:n_trials %for each trial, backproject drift by sac_shift
        cur_inds = gen_data.trial_boundaries(ii,1):gen_data.trial_boundaries(ii,2);
        drift_post_cor(cur_inds(1:end - neural_delay)) = drift_post_cor(cur_inds(neural_delay + 1:end));
    end
    drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids)),1:NT); %interpolate
    drift_post_cor(isnan(drift_post_cor)) = 0;
    
    if HMM_params.n_fix_iter == 0
        fix_post_mean_cor = zeros(NT,1);
    end
    cur_fix_shifts = round(fix_post_mean_cor + drift_post_cor); %total shift to apply to raw stimulus (both fixation and drift corrections)
    
    shift_stim_mat = stim_mat;
    for ii=1:NT %apply shifts to stimulus
        shift_stim_mat(gen_data.used_inds(ii),:) = shift_matrix_Nd(stim_mat(gen_data.used_inds(ii),:),-cur_fix_shifts(ii),2);
    end
    all_Xmat = NIM.create_time_embedding(shift_stim_mat,gen_data.stim_params);
    if add_usfac %project onto spatial TBs if needed
        X{1} = tb_proc_stim(all_Xmat(gen_data.used_inds,gen_data.use_kInds),add_usfac,flen);
    else
        X{1} = all_Xmat(gen_data.used_inds,gen_data.use_kInds);
    end
    X(2:1+n_add_covariates) = nonstim_X;
    
    %% REFIT ALL CELLS
    silent = 1;
    for ss = 1:n_units
        cur_fit_inds = gen_data.modfit_inds(~isnan(Robs_mat(gen_data.modfit_inds,ss))); %all non-rpt indices when current unit was isolated
        
        drift_it_mods{nn+1}(ss) = drift_it_mods{nn}(ss);
        drift_it_mods{nn+1}(ss).spkNL.params = [1 1]; %set alpha and beta to their default values when fitting filters to avoid issues with regularization
        drift_it_mods{nn+1}(ss) = drift_it_mods{nn+1}(ss).fit_filters(Robs_mat(:,ss),X,cur_fit_inds,'silent',silent,'sub_inds',[]); %fit constant offset first to get back our original model (sans other spkNL params)
        drift_it_mods{nn+1}(ss) = drift_it_mods{nn+1}(ss).fit_filters(Robs_mat(:,ss),X,cur_fit_inds,'silent',silent); %fit stimulus filters
        drift_it_mods{nn+1}(ss) = drift_it_mods{nn+1}(ss).fit_spkNL(Robs_mat(:,ss),X,cur_fit_inds,'silent',silent);  %refit spk NL
        newLL = drift_it_mods{nn+1}(ss).eval_model(Robs_mat(:,ss),X,cur_fit_inds);
        drift_it_LLimp(nn+1,ss) = (newLL - init_mods(ss).fit_props.null_LL)/log(2);
        if nn > 1
            fprintf('Unit %d LLimps.  Original: %.4f  Prev: %.4f  New: %.4f\n',ss,drift_it_LLimp(1,ss),drift_it_LLimp(nn,ss),drift_it_LLimp(nn+1,ss));
        else
            fprintf('Unit %d LLimps. Original: %.4f  New: %.4f\n',ss,drift_it_LLimp(1,ss),drift_it_LLimp(nn+1,ss));
        end
    end
end

%%
mod_fits.fix_iters = fix_it_mods;
mod_fits.drift_iters = drift_it_mods;
EP_post.fix_mean = fix_it_post_mean;
EP_post.fix_std = fix_it_post_std;
EP_post.drift_mean = drift_it_post_mean;
EP_post.drift_std = drift_it_post_std;

