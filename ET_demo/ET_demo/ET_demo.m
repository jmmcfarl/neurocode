clear all
ET_demo_dir = '~/James_scripts/ET_demo'; %set this to the directory containing the ET_demo folder
cd(ET_demo_dir);
addpath(genpath(ET_demo_dir)); %add helper functions to path
load('simDATA.mat');

init_mods_name = 'init_mod_data.mat'; %name of file to save initial models

%% SET MODEL PARAMETERS
useFrac = 1; %(1) fraction of total data to use (1 = all; 0.5 = first 50%)
use_coil = false; %use or don't use information from the (simulated) eye coil signal

flen = 12; %(12) number of time lags to include in the time-embedded stimulus matrix (units of dt)
use_nPix = 16; %(16) number of pixels to include in the stimulus models
spatial_usfac = 2; %(2) spatial up-sample factor for initial model-fitting, and estimation of fixation-corrections. Only set up for multiples of 2.
hr_usfac = 4; % (4) spatial up-sample factor for estimating drift (helpful to use a higher resolution here).

n_fix_inf_it = 3; %(3) number of iterations for estimating fixation corrections
n_drift_inf_it = 1; %(1) number of iterations for estimating drift corrections

fix_prior_sigma = 0.15; %(0.15) SD of Gaussian prior on the fixation corrections (in deg)
drift_prior_sigma = 0.004; %(0.004) SD of Gaussian prior on the change in eye position between adjacent time points
drift_jump_sigma = 0.075;  %(0.075) SD of Gaussian prior on the eye position deviation from the fixation avg
drift_dsf = 2; %(2) temporal down-sample factor for inferring drift corrections (relative to base time resolution dt).

fix_noise_sigma = 0.1; %(0.1) assumed noise SD of coils, in deg (for measuring change in fixation-avg position)
drift_noise_sigma = 0.004; %(0.004) assumed noise SD of coils (for measuring change in position between adjacent time bins).

recompute_initial_models = false; %load in precomputed initial models if they exist in the current dir
%% PRE-COMPUTATIONS
bar_ori = 0; %orientation of the bar stimulus (so we will infer the orthoganol component of eye position)
base_dx = 0.0565; %width of bar stimuli in deg
sac_shift = round(0.05/dt); %the average time-lag of V1 responses (50ms here). This is the assumed time-lag of neural responses to a change in eye position

beg_buffer = flen*dt; %length of data segment to ignore at the beginning of each trial (to avoid transient effects)
full_nPix = size(stim_mat,2); %number of stimulus pixels per frame

%up-sampled versions of use_nPix klen and full_nPix
use_nPix_us = use_nPix*spatial_usfac;
use_nPix_hr = use_nPix*hr_usfac;
klen_us = use_nPix_us*flen;
klen_hr = use_nPix_hr*flen;
full_nPix_us = spatial_usfac*full_nPix;
full_nPix_hr = hr_usfac*full_nPix;

%up-sampled pixel widths
sp_dx = base_dx/spatial_usfac;
sp_dx_hr = base_dx/hr_usfac;

%create parameter structures for the stimulus models
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
stim_params_us = NIMcreate_stim_params([flen full_nPix_us],dt);
stim_params_hr = NIMcreate_stim_params([flen full_nPix_hr],dt);

max_shift = round(15*spatial_usfac); %maximum image translation considered for fixation corrections (in corresponding pixel units)
max_Dshift = round(8*hr_usfac); %same as above for drift corrections
dshift = 1; %lattice spacing in pixels for discretization of eye positions

%vector of possible eye positions (equivalently image shifts) for fixation
%corrections
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);

%same for drift corrections
Dshifts = -max_Dshift:dshift:max_Dshift;
n_Dshifts = length(Dshifts);

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

[Xinds_hr,~] = meshgrid(1/hr_usfac:1/hr_usfac:full_nPix,1:flen);
cur_use_pix = (1/hr_usfac:1/hr_usfac:use_nPix) + buffer_pix;
use_kInds_hr = find(ismember(Xinds_hr(:),cur_use_pix));

%generate shift matrices for fixation-corrections. Must be applied to the stimulus (not the filters)
%these are stored as a cell array of sparse matrices
It = speye(flen);
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = spdiags( ones(full_nPix_us,1), -shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    shift_mat{xx} = temp(:,use_kInds_up);
end

%generate shift matrices for drift-corrections. Must be applied to the stimulus (not the filters)
Dshift_mat = cell(n_Dshifts,1);
for xx = 1:n_Dshifts
    temp = spdiags( ones(full_nPix_hr,1), -Dshifts(xx), full_nPix_hr, full_nPix_hr);
    temp = kron(temp,It);
    Dshift_mat{xx} = temp(:,use_kInds_hr);
end

%log-prior on fixation corrections
fix_Lprior = -(shifts*sp_dx).^2./(2*fix_prior_sigma^2);
fix_Lprior = fix_Lprior - logsumexp(fix_Lprior); %normalize

%log prior on the drift corrections at boundary points
drift_jump_Lprior = -(Dshifts*sp_dx_hr).^2./(2*drift_jump_sigma^2);
drift_jump_Lprior = drift_jump_Lprior - logsumexp(drift_jump_Lprior); %normalize

%log prior on drift
cdist = squareform(pdist(Dshifts'*sp_dx_hr));
drift_lA = -cdist.^2/(2*(drift_prior_sigma*drift_dsf)^2);
drift_lA = bsxfun(@minus,drift_lA,logsumexp(drift_lA,2)); %normalize

%% Create up-sampled stimulus matrix
if spatial_usfac > 1
    stimmat_up = zeros(size(stim_mat,1),full_nPix_us);
    for ii = 1:full_nPix
        for jj = 1:spatial_usfac
            stimmat_up(:,spatial_usfac*(ii-1)+jj) = stim_mat(:,ii);
        end
    end
else
    stimmat_up = stim_mat;
end

%make time embedded X matrices
all_Xmat = create_time_embedding(stim_mat,stim_params);
all_Xmat_us = create_time_embedding(stimmat_up,stim_params_us);

%% DEFINE AND PARSE DATA USED FOR ANALYSIS
n_trials = length(trial_start_inds);
trial_ids = nan(length(sim_eyepos),1);
t_since_start = nan(length(sim_eyepos),1); %stores time since trial onset (sec)
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
    t_since_start(cur_inds) = (0:length(cur_inds)-1)*dt;
end

used_inds = find(t_since_start >= beg_buffer); %exclude initial beg_buffer sec from each trial
NT = length(used_inds);

%if restricting to subset of data, find next trial end point to stop on
useLen = ceil(NT*useFrac); %number of samples to use
if NT > useLen
    endpt = trial_end_inds(find(trial_end_inds >= useLen,1,'first'));
    used_inds((endpt+1):end) = [];
    NT = length(used_inds);
end

%compute saccade start/stop indices and trial start/stop indices relative
%to new vector (after exclusion)
[C1,IA1,IB1] = intersect(used_inds,saccade_start_inds);
[C2,IA2,IB2] = intersect(used_inds,saccade_stop_inds);

not_present = find(~ismember(IB2,IB1)); %these are instances where a saccade started during the excluded portion of a trial but ended during included portion
sac_beg_trials = trial_ids(used_inds(IA2(not_present))); %these are the trials where this happens
new_tstart = IA2(not_present);%will make the new trial onset at the time of the saccade ending
IA2(not_present) = []; %eliminate these saccade end times

%saccade start/stop inds relative to the used data vector
saccade_stop_inds = IA2;
saccade_start_inds = IA1;

%trial start/stop inds relative to the used data vector
trial_ids = trial_ids(used_inds);
trial_start_inds = [1; find(diff(trial_ids) ~= 0) + 1];
trial_end_inds = [find(diff(trial_ids) ~= 0); NT];

%adjust trial onsets when there was a saccade in process.
trial_start_inds(sac_beg_trials) = new_tstart;

%define start/stop indices for each fixation
fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds;
%get rid of double-counted fixations (due to trial onsets coincident with
%saccade endings
fix_start_inds(fix_durs <= 0) = [];
fix_stop_inds(fix_durs <= 0) = [];
n_fixs = length(fix_start_inds);

%push the fixation markers forward in time to account for neural time lags
%(pfix* will be shifted start/stop inds)
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

%create vectors of fixation indices, shifted and un-shifted
fix_ids = nan(NT,1);
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
    cur_inds = pfix_start_inds(ii):(pfix_stop_inds(ii));
    pfix_ids(cur_inds) = ii;
end

%% Create simulated eye coil signal if needed
if use_coil
   coil_fix_noise_SD = 0.075; %noise SD for the simulated coil (one noise term per fixation). 
   coil_noise = randn(n_fixs,1)*coil_fix_noise_SD;
   coil_EP = sim_eyepos(used_inds);
   coil_EP(~isnan(fix_ids)) = coil_EP(~isnan(fix_ids)) + coil_noise(fix_ids(~isnan(fix_ids)));
end

%% PROCESS DATA FOR MODEL FITTING
Robs_mat = double(Robs_mat(used_inds,:)); %convert spike count data to doubles
n_units = size(Robs_mat,2); %number of units

%keep only used time bins
all_Xmat = all_Xmat(used_inds,:);
all_Xmat_us = all_Xmat_us(used_inds,:);

%create Xmats using only central pixels
Xmat = all_Xmat(:,use_kInds);
Xmat_us = all_Xmat_us(:,use_kInds_up);

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
silent = 1; %suppress output of model-estimation routine
base_lambda_d2XT = 10; %[10] spatiotemporal smoothness regularization strength
base_lambda_L1 = 2; %[2] sparseness regularization strength

%parameter structures defining the stimulus
init_stim_params = NIMcreate_stim_params([flen use_nPix],dt);
fin_stim_params = NIMcreate_stim_params([flen use_nPix_us],dt);

n_squared_filts = 2; %number of squared stimulus filters to use
mod_signs = ones(1,n_squared_filts+1); %assume all squared filters are 'excitatory' here (sign doesn't matter for linear filter)
NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts)]; %one linear filter and n_squared_filts squared filters
init_d2XT = 1; %[1] weak initial smoothness regularization
init_reg_params = NIMcreate_reg_params('lambda_d2XT',init_d2XT);

%works best to set the initial optimization tolerances lower so initial
%estimates of filters aren't too far off before adjusting regularization
%parameters
init_optim.optTol = 1e-6; init_optim.progTol = 1e-8;

if ~exist(init_mods_name,'file') || recompute_initial_models %if computing initial models

%initialization
init_filts = zeros(klen_us,length(mod_signs));
null_LL = nan(n_units,1);
all_mod_R2 = nan(n_units,1);
all_mod_LLimp = nan(n_units,1);
for ss = 1:n_units
    fprintf('Computing base model for Unit %d of %d\n',ss,n_units);
    cur_Robs = Robs_mat(:,ss); %binned spikes for current unit
    
    %initial model estimate at base stimulus resolution
    gqm1 = NIMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params);
    gqm1 = NIMfit_filters(gqm1,cur_Robs,Xmat,[],[],silent,init_optim);
    
    %spatial up-sampling of initial filter estimates
    base_filts = reshape([gqm1.mods(:).filtK],[flen use_nPix n_squared_filts+1]);
    if spatial_usfac > 1
        base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
        for ii = 1:use_nPix
            for jj = 1:spatial_usfac
                base_filts_up(:,spatial_usfac*(ii-1)+jj,:) = 0.5*base_filts(:,ii,:);
            end
        end
    else
        base_filts_up = base_filts;
    end
    base_filts_up = reshape(base_filts_up,use_nPix_us*flen,n_squared_filts+1);
    
    %initialize up-sampled stimulus model
    for ii = 1:n_squared_filts+1
        init_filts(:,ii) = base_filts_up(:,ii);
    end
    gqm2 = NIMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_filts);
    gqm2.spk_NL_params(1) = gqm1.spk_NL_params(1);
    
    %adjust regularization strength to match estimated scale of filters
    [LL, penLL, pred_rate, G, gint] = NIMmodel_eval(gqm2,cur_Robs,Xmat_us);
    std_gint = std(gint); %scale of filter outputs
    for ii = 1:length(mod_signs)
        gqm2 = NIMadjust_regularization(gqm2,ii,'lambda_d2XT',base_lambda_d2XT/std_gint(ii)^2,'lambda_L1',base_lambda_L1/std_gint(ii));
    end
    
    %estimate model params
    gqm2 = NIMfit_filters(gqm2,cur_Robs,Xmat_us,[],[],silent);
    all_mod_fits(ss) = gqm2; %store model
    all_mod_fits_spkNL(ss) = NIMfit_logexp_spkNL(gqm2,cur_Robs,Xmat_us); %estimate spiking NL parameters
    
    %compute model performance measures
    [LL,~,pred_rate] = NIMmodel_eval(all_mod_fits_spkNL(ss),cur_Robs,Xmat_us);
    null_prate = ones(NT,1)*mean(cur_Robs); %prediction of the null model is a constant rate in this case
    null_LL(ss) = sum(cur_Robs.*log(null_prate) - null_prate)/sum(cur_Robs);
    
    all_mod_R2(ss) = pseudo_r2(cur_Robs,pred_rate,null_prate); %'pseudo-R2', a likelihood based R2 measure
    all_mod_LLimp(ss) = (LL-null_LL(ss))/log(2); %log-likelihood improvement over null model (in bits/spk)
end

save(init_mods_name,'all_mod_*','null_LL');
else
    fprintf('Loading pre-computed initial models\n');
    load(init_mods_name);
end
clear Xmat all_Xmat Xmat_us
%% COMPUTE INFO ABOUT CHANGE IN COIL-MEASURED EYE POSITION BETWEEN FIXATIONS (IF USING COIL)
if use_coil
    coil_fix_avg = nan(n_fixs,1); %average coil-measured eye-position within each fixation
    for ii = 1:n_fixs
        cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
        coil_fix_avg(ii) = mean(coil_EP(cur_inds));
    end
    
    %change in fixation-avg eye position between adjacent fixations
    coil_fix_deltas = nan(n_fixs,1);
    coil_fix_deltas(2:end) = diff(coil_fix_avg);
end
%% ITERATIVE ESTIMATION OF FIXATION-BASED CORRECTIONS
%initialization
it_mods{1} = all_mod_fits;
it_mods_spkNL{1} = all_mod_fits_spkNL;
it_LLimp(1,:) = all_mod_LLimp;
it_R2(1,:) = all_mod_R2;
it_fix_post_mean = nan(n_fix_inf_it,n_fixs);
it_fix_post_std = nan(n_fix_inf_it,n_fixs);
for nn = 1:n_fix_inf_it
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,n_fix_inf_it);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_units,klen_us,n_squared_filts+1);
    mod_spkNL_params = nan(n_units,3);
    for ss = 1:n_units
        cur_k = [it_mods{nn}(ss).mods(:).filtK];
        filt_bank(ss,:,:) = cur_k; %store all stimulus filters
        mod_spkNL_params(ss,:) = it_mods_spkNL{nn}(ss).spk_NL_params; %store all spiking NL parameters
    end
    filt_bank = permute(filt_bank,[2 1 3]);
    
    %% ESTIMATE LL for each shift in each stimulus frame
    frame_LLs = nan(NT,n_shifts); %total LL (summed across units) for each shift
    for xx = 1:n_shifts
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_stim_shift = all_Xmat_us*shift_mat{xx}; %shifted stimulus matrix
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = bsxfun(@times,ones(NT,n_units),mod_spkNL_params(:,1)'); %constant offset
        gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1)); %output of stimulus filters
        for ff = 2:(n_squared_filts+1)
            gfuns = gfuns + mod_signs(ff)*(cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2; %square the output of the quadratic stim filters
        end
        
        %incorporate spkNL parameter beta
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
        
        %handle numerical overflow with log(1+exp)
        too_large = gfuns > 50;
        pred_rate = log(1+exp(gfuns));
        pred_rate(too_large) = gfuns(too_large);
        
        %incorporate spkNL parameter alpha
        pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
        
        %enforce minimum predicted rate
        pred_rate(pred_rate < 1e-50) = 1e-50;
        
        %compute LLs and sum across units at each time
        frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
    end
    
    %% INFER SEQUENCE OF FIXATION CORRECTIONS
    fix_LLs = nan(n_fixs,n_shifts); %LL for each shift summed across all time points in each fixation
    for ii = 1:n_fixs
        cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii); %time points of the current fixation (projected forward in time to account for neural delay)
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
    end
    
    if ~use_coil %if not using coil info, dont need to use forward-backward algo
        lgamma = bsxfun(@plus,fix_LLs,fix_Lprior); %log-likelihood plus log-prior
        
    else %if using coil info need to implement forward-backward algo
        
        %initialize forward and backward variables
        lalpha=zeros(n_fixs,n_shifts);
        lbeta = zeros(n_fixs,n_shifts);
        
        %compute forward messages
        lalpha(1,:) = fix_Lprior + fix_LLs(1,:);
        for t=2:n_fixs
            %compute state-transition matrix given coil-measured fix delta
            cdist = pdist2(shifts'*sp_dx + coil_fix_deltas(t),shifts'*sp_dx);
            cur_lA = -cdist.^2/(2*fix_noise_sigma^2);
            cur_lA = bsxfun(@plus,cur_lA,fix_Lprior);
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
            
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + fix_LLs(t,:);
        end
        
        %compute backward messages
        lbeta(n_fixs,:)=log(ones(1,n_shifts));
        for t=n_fixs-1:-1:1
            cdist = pdist2(shifts'*sp_dx + coil_fix_deltas(t+1),shifts'*sp_dx);
            cur_lA = -cdist.^2/(2*fix_noise_sigma^2);
            cur_lA = bsxfun(@plus,cur_lA,fix_Lprior);
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
            
            lf1 = lbeta(t+1,:) + fix_LLs(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA');
        end
        lgamma= lalpha + lbeta;
    end
        
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2)); %normalize
    gamma = exp(lgamma); %probability of each eye position at each time
    
    %estimate posterior mean and SD
    it_fix_post_mean(nn,:) = sum(bsxfun(@times,gamma,shifts),2);
    cur_diff = bsxfun(@minus,it_fix_post_mean(nn,:)',shifts).^2;
    it_fix_post_std(nn,:) = sqrt(sum(cur_diff.*gamma,2));
    
    %% COMPUTE ESTIMATED EYE POSITION AND CORRECTED RETINAL STIMULUS
    %estimated eye position given fixation corrections
    fix_EP = nan(NT,1);
    fix_EP(~isnan(fix_ids)) = it_fix_post_mean(nn,fix_ids(~isnan(fix_ids)));
    fix_EP = round(interp1(find(~isnan(fix_ids)),fix_EP(~isnan(fix_ids)),1:NT)); %interpolate through saccade times
    
    %now correct for inferred eye positions
    shift_stimmat_up = stimmat_up;
    for i=1:NT
        shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(shift_stimmat_up(used_inds(i),:),-fix_EP(i),2);
    end
    
    %create corrected time-embedded stim-mat
    Xmat_up_fixcor = create_time_embedding(shift_stimmat_up,stim_params_us);
    Xmat_up_fixcor = Xmat_up_fixcor(used_inds,use_kInds_up);
    
    %% REESTIMATE MODELS FOR ALL UNITS
    silent = 1;
    for ss = 1:n_units
        fprintf('Refitting model for unit %d of %d\n',ss,n_units);
        cur_Robs = Robs_mat(:,ss);
        
        %initialize from previous iteration model and reestimate given
        %new Xmat
        it_mods{nn+1}(ss) = it_mods{nn}(ss);
        it_mods{nn+1}(ss) = NIMfit_filters(it_mods{nn+1}(ss),cur_Robs,Xmat_up_fixcor,[],[],silent); %fit stimulus filters
        
        %refit spk NL
        it_mods_spkNL{nn+1}(ss) = NIMfit_logexp_spkNL(it_mods{nn+1}(ss),cur_Robs,Xmat_up_fixcor);
        
        %compute performance
        [newLL,~,new_prate] = NIMmodel_eval(it_mods_spkNL{nn+1}(ss),cur_Robs,Xmat_up_fixcor);
        null_prate = ones(NT,1)*mean(cur_Robs);
        it_R2(nn+1,ss) = pseudo_r2(cur_Robs,new_prate,null_prate);
        it_LLimp(nn+1,ss) = (newLL - null_LL(ss))/log(2);
        
        %display LL improvements
        if nn > 1
            fprintf('Original: %.4f  Prev: %.4f  New: %.4f\n',all_mod_LLimp(ss),it_LLimp(nn,ss),it_LLimp(nn+1,ss));
        else
            fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(ss),it_LLimp(nn+1,ss));
        end
    end
end

clear Xmat_up_fixcor all_Xmat_us

%% COMPUTE UP-SAMPLED STIMULUS MATRIX
if hr_usfac > spatial_usfac %if using further spatial up-sampling
    stimmat_hr = zeros(size(stim_mat,1),full_nPix_hr);
    for ii = 1:full_nPix
        for jj = 1:hr_usfac
            stimmat_hr(:,hr_usfac*(ii-1)+jj) = stim_mat(:,ii);
        end
    end
else
    stimmat_hr = stimmat_up;
end

%% INCORPORATE BEST ESTIMATES OF FIXATION CORRECTIONS
%estimated eye position given fixation corrections
fix_EP = nan(NT,1);
fix_EP(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fix_EP = interp1(find(~isnan(fix_ids)),fix_EP(~isnan(fix_ids)),1:NT);

%estimated eye position at new spatial resolution
fix_EP = fix_EP*hr_usfac/spatial_usfac;

%correct new stimulus for estimated eye positions
shift_stimmat_up = stimmat_hr;
for i=1:NT
    shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(shift_stimmat_up(used_inds(i),:),-round(fix_EP(i)),2);
end

%create new time embedded Xmat
all_Xmat_us = create_time_embedding(shift_stimmat_up,stim_params_hr);
all_Xmat_us = all_Xmat_us(used_inds,:);
Xmat_us = all_Xmat_us(:,use_kInds_hr);

%% ESTIMATE MODEL PARAMETERS FOR NEW STIMULUS
%initialize up-sampled data
all_mod_fits_hr = it_mods{end};
all_mod_fits_hr_spkNL = it_mods_spkNL{end};
all_mod_R2_hr = it_R2(end,:);
all_mod_LLimp_hr = it_LLimp(end,:);

init_filts = zeros(klen_hr,length(mod_signs));
fin_stim_params = NIMcreate_stim_params([flen use_nPix_hr],dt);
if hr_usfac > spatial_usfac %if using further spatial up-sampling
    for ss = 1:n_units
        fprintf('Computing model for Unit %d of %d\n',ss,n_units);
        cur_Robs = Robs_mat(:,ss);
        
        %spatial up-sampling of filter estimates
        base_filts = reshape([it_mods{end}(ss).mods(:).filtK],[flen use_nPix_us n_squared_filts+1]);
        base_filts_up = zeros(flen,use_nPix_hr,n_squared_filts+1);
        for ii = 1:use_nPix_us
            for jj = 1:hr_usfac/spatial_usfac
                base_filts_up(:,hr_usfac/spatial_usfac*(ii-1)+jj,:) = 0.5*base_filts(:,ii,:);
            end
        end
        base_filts_up = reshape(base_filts_up,use_nPix_hr*flen,n_squared_filts+1);
        
        %initialize up-sampled stimulus filters
        for ii = 1:n_squared_filts+1
            init_filts(:,ii) = base_filts_up(:,ii);
        end
        gqm2 = NIMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_filts);
        
        %adjust reglurization strength to new resolution
        [LL, penLL, pred_rate, G, gint] = NIMmodel_eval(gqm2,cur_Robs,Xmat_us);
        sgint = std(gint); sgint(sgint == 0) = 1;
        for ii = 1:length(mod_signs)
            gqm2 = NIMadjust_regularization(gqm2,ii,'lambda_d2XT',base_lambda_d2XT/sgint(ii)^2,'lambda_L1',base_lambda_L1/sgint(ii));
        end
        all_mod_fits_hr(ss) = gqm2;
        all_mod_fits_hr(ss).spk_NL_params(1) = it_mods{end}(ss).spk_NL_params(1);
        
        %estimate model parameters
        all_mod_fits_hr(ss) = NIMfit_filters(all_mod_fits_hr(ss),cur_Robs,Xmat_us);
        
        %estimate spk NL parameters
        all_mod_fits_hr_spkNL(ss) = NIMfit_logexp_spkNL(all_mod_fits_hr(ss),cur_Robs,Xmat_us);
        
        %model improvements
        [LL,~,pred_rate] = NIMmodel_eval(all_mod_fits_hr_spkNL(ss),cur_Robs,Xmat_us);
        all_mod_R2_hr(ss) = pseudo_r2(cur_Robs,pred_rate,null_prate);
        all_mod_LLimp_hr(ss) = (LL-null_LL(ss))/log(2);
    end
end

clear Xmat_us
%% COMPUTE INFO ABOUT COIL MEASURED DRIFT (IF USING COIL)
if use_coil
    coil_drift = nan(NT,1); %average coil-measured eye-position within each fixation
    for ii = 1:n_fixs
        cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
        coil_drift(cur_inds) = coil_EP(cur_inds) - mean(coil_EP(cur_inds));
    end
    
    %coil measured drift velocity
    coil_drift_vel = [0; diff(coil_drift)];
    
    %compute posterior mean and SD of the drift velocity, given the
    %measured coil drift, assumed noise variance, and prior
    post_drift_sigma = sqrt(1/(1/drift_noise_sigma^2 + 1/drift_prior_sigma^2));
    post_mean_drift = post_drift_sigma^2*coil_drift_vel/drift_noise_sigma^2;
    
    %forward-project the drift velocities to align with the time-lagged
    %neural responses
    for ii = 1:n_fixs
        cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
        if length(cur_inds) > sac_shift
            post_mean_drift(cur_inds(sac_shift+1:end)) = post_mean_drift(cur_inds(1:(end-sac_shift)));
        end
    end
end

%% NOW INFER DRIFT CORRECTIONS

dit_mods{1} = all_mod_fits_hr;
dit_mods_spkNL{1} = all_mod_fits_hr_spkNL;
dit_LLimp(1,:) = all_mod_LLimp_hr;
dit_R2(1,:) = all_mod_R2_hr;
drift_post_mean = nan(n_drift_inf_it,NT);
drift_post_std = nan(n_drift_inf_it,NT);
for nn = 1:n_drift_inf_it
    fprintf('Inferring drift corrections, iter %d of %d\n',nn,n_drift_inf_it);
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_units,klen_hr,n_squared_filts+1);
    mod_spkNL_params = nan(n_units,3);
    for ss = 1:n_units
        cur_k = [dit_mods{nn}(ss).mods(:).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = dit_mods_spkNL{nn}(ss).spk_NL_params;
    end
    filt_bank = permute(filt_bank,[2 1 3]);
    
    %% ESTIMATE LL for each shift in each stimulus frame
    
    frame_LLs = nan(NT,n_Dshifts);
    for xx = 1:length(Dshifts)
        fprintf('Dshift %d of %d\n',xx,n_Dshifts);
        cur_stim_shift = all_Xmat_us*Dshift_mat{xx}; %shifted stimulus matrix
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = bsxfun(@times,ones(NT,n_units),mod_spkNL_params(:,1)');
        gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1));
        for ff = 2:(n_squared_filts+1)
            gfuns = gfuns + mod_signs(ff)*(cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
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
    for ff = 1:n_fixs %estimate drift within each fixation
        if mod(ff,100)==0
            fprintf('Fixation %d of %d\n',ff,n_fixs);
        end
        
        %indices corresponding to the current fixation (projected forward
        %in time to account for neural delay)
        tset = find(pfix_ids==ff)';
        ntset = length(tset);
        
        if ntset > drift_dsf
            %down-sample in time by drift_dsf to estimate drift corrections
            nt_pts = ceil(ntset/drift_dsf); %number of corrections
            tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf); %index values of down-sampled time bins
            tpt_loc(end) = ntset;
            
            %initialize forward and backward variables
            talpha = zeros(nt_pts,n_Dshifts);
            tbeta = zeros(nt_pts,n_Dshifts);
            
            cur_LL_set = frame_LLs(tset,:); %set of LLs for given time bins
            if mod(ntset,drift_dsf) ~= 0 %if the last time bin doesn't have drift_dsf points add some 'zero-padding'
                dangling_pts = nt_pts*drift_dsf-ntset;
                cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_Dshifts));
            end
            cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_Dshifts]);
            
            %sum LLs within each NEW time bin
            cur_LL_set = squeeze(sum(cur_LL_set,1));
            
            %if using coil info, compute the average coil-measured drift
            %vel, at the lower time resolution
            if use_coil
                cur_drift_vel = post_mean_drift(tset);
                if mod(ntset,drift_dsf) ~= 0 %if the last time bin doesn't have drift_dsf points add some 'zero-padding'
                    cur_drift_vel = cat(1,cur_drift_vel,nan(dangling_pts,1));
                end
                cur_drift_vel = drift_dsf*nanmean(reshape(cur_drift_vel,[drift_dsf nt_pts]));
            end
            
            if ~use_coil  %if not using coil then use a fixed state-transition matrix
                cur_drift_lA = repmat(drift_lA,[1 1 nt_pts]);
            else %if using coil, compute the state-transition matrix for each time bin
                cur_drift_lA = nan(n_Dshifts,n_Dshifts,nt_pts);
                for iii = 1:nt_pts
                    if ~isnan(cur_drift_vel(iii))
                        %log of gaussian state-transition matrix, with
                        %estimated posterior mean and SD
                        cdist = pdist2(Dshifts'*sp_dx_hr + cur_drift_vel(iii),Dshifts'*sp_dx_hr);
                        cur_drift_lA(:,:,iii) = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
                    else
                        cur_drift_lA(:,:,iii) = drift_lA;
                    end
                end
                
                %normalize state transition matrices
                cur_drift_lA = bsxfun(@minus,cur_drift_lA,logsumexp(cur_drift_lA,2));
            end
            
            %compute forward variables
            talpha(1,:) = drift_jump_Lprior + cur_LL_set(1,:); %initial time bin is given by a sum of the log-prior and the LL (independent of previous times)
            for t = 2:nt_pts
                talpha(t,:) = logmulexp(talpha(t-1,:),cur_drift_lA(:,:,t)) + cur_LL_set(t,:);
            end
            
            %compute backward variables
            tbeta(end,:)=log(ones(1,n_Dshifts));
            for t = (nt_pts-1):-1:1
                lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
                tbeta(t,:) = logmulexp(lf1,cur_drift_lA(:,:,t+1)');
            end
            temp_gamma = talpha + tbeta; %combine forward and backward variables
            
            if drift_dsf > 1
                if nt_pts > 1
                    %interpolate the log-posterior onto the original time
                    %axis
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
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2)); %normalize log posterior
    gamma = exp(lgamma); %compute posterior
    
    %compute posterior mean and SD
    drift_post_mean(nn,:) = sum(bsxfun(@times,gamma,Dshifts),2);
    drift_post_std(nn,:) = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2) - drift_post_mean(nn,:)'.^2);
        
    %% construct drift-corrected X-mat
    drift_post_cor = drift_post_mean(nn,:);
    
    %back-project drift corrections within each trial by sac_shift
    for ii = 1:length(trial_start_inds)
        cur_inds = trial_start_inds(ii):trial_end_inds(ii);
        drift_post_cor(cur_inds(1:end-sac_shift)) = drift_post_cor(cur_inds(sac_shift+1:end));
    end
    %interpolate through saccades
    drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids)),1:NT);
    drift_post_cor(isnan(drift_post_cor)) = 0;
    
    %total eye position is sum of fixation corrections and drift
    %corrections
    all_post_EP = fix_EP+drift_post_cor;
    
    %Compute corrected Xmat
    for i=1:NT
        shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(stimmat_hr(used_inds(i),:),-round(all_post_EP(i)),2);
    end
    
    %create corrected time-embedded stim mat
    cur_X = create_time_embedding(shift_stimmat_up,stim_params_hr);
    cur_X = cur_X(used_inds,use_kInds_hr);
    
    %% REFIT MODELS
    silent = 1;
    for ss = 1:n_units
        fprintf('Refitting model for tr cell %d of %d\n',ss,n_units);
        cur_Robs = Robs_mat(:,ss);
        
        %initialize new model and estimate 
        dit_mods{nn+1}(ss) = dit_mods{nn}(ss);
        dit_mods{nn+1}(ss) = NIMfit_filters(dit_mods{nn+1}(ss),cur_Robs,cur_X,[],[],silent); %fit stimulus filters
        dit_mods_spkNL{nn+1}(ss) = NIMfit_logexp_spkNL(dit_mods{nn+1}(ss),cur_Robs,cur_X); %estimate spk NL params
        
        [newLL,~,new_prate] = NIMmodel_eval(dit_mods_spkNL{nn+1}(ss),cur_Robs,cur_X);
        null_prate = ones(NT,1)*mean(cur_Robs);
        dit_R2(nn+1,ss) = pseudo_r2(cur_Robs,new_prate,null_prate);
        dit_LLimp(nn+1,ss) = (newLL - null_LL(ss))/log(2);
        fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(ss),dit_LLimp(nn,ss),dit_LLimp(nn+1,ss));
    end
end

%%
save demo_output_fulldata2 dit_* it_* drift_* spatial_usfac hr_usfac sp_dx sp_dx_hr fix_ids trial_*inds NT sac_shift use_coil

%% COMPUTE FINAL EYE POSITION ESTIMATE
%initialize
fin_fix_mean = nan(NT,1);
fin_fix_std = nan(NT,1);

%best estimate of fixation corrections
fin_fix_mean(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));

%interpolate through saccades
fin_fix_mean = interp1(find(~isnan(fix_ids)),fin_fix_mean(~isnan(fix_ids)),1:NT);
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

%convert to degrees (accounting for potentially different spatial
%resolutions used to estimate drift and fixation corrections)
fin_fix_mean = fin_fix_mean*sp_dx_hr*hr_usfac/spatial_usfac;
fin_fix_std = fin_fix_std*sp_dx_hr*hr_usfac/spatial_usfac;

%initialize drift estimates
fin_drift_mean = drift_post_mean(1,:)*sp_dx_hr;
fin_drift_std = drift_post_std(1,:)*sp_dx_hr;

%back-project eye positions by neural delay within each trial
fix_inds = [];
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fix_inds = [fix_inds cur_inds];
    fin_drift_mean(cur_inds(1:end-sac_shift)) = fin_drift_mean(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end

%interpolate through saccades
fin_drift_mean = interp1(find(~isnan(fix_ids)),fin_drift_mean(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

%combine fixation and drift correction estimates
fin_tot_corr = fin_fix_mean + fin_drift_mean;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);

%% COMPARE ESTIMATED AND TRUE EYE POSITION DATA
close all
yl = [-0.4 0.4]; %range of y-axis (in deg)

t_axis = (1:NT)*dt;
H = figure(); 
for tt = 1:n_trials
    fprintf('Trial %d of %d\n',tt,n_trials);
    uu = find(trial_ids == tt);
    if ~isempty(uu)
        
        %beginning and end times of current trial
        bt = t_axis(uu(1));
        et = t_axis(uu(end));
        
        %saccade times 
        cur_sac_inds = find(ismember(saccade_start_inds,uu));
        rel_sac_start_times = t_axis(saccade_start_inds(cur_sac_inds)) - bt;
        
        h1=shadedErrorBar(t_axis(uu)-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','k'},0);
        hold on;
        h2=plot(t_axis(uu)-bt,sim_eyepos(used_inds(uu)),'r');
        line([0 et-bt],[0 0],'color','k','linestyle','--');
        xlim([0 et-bt]);
        ylim(yl);
        xlabel('Time (s)');
        ylabel('Orthoganol position (deg)');
        title(sprintf('Trial %d',tt));
        
        %make a vertical line at the start time of each saccade
        for ii = 1:length(rel_sac_start_times)
            line(rel_sac_start_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
        end
        
        delete(h1.edge);
        pause
        clf(H);
    end
end
