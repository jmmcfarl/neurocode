function [fix_it_mods,fix_post_data,drift_it_mods,drift_post_data] = eyetracking_EM(...
    init_mods,Robs_mat,stim_mat,nonstim_X,tr_units,HMM_data,ET_meas,gen_data)


% init_mods: array of initial NIM model objects
% Robs_mat: matrix of observed spike counts across all units (Nans when a unit was not isolated)
% stim_mat: raw stimulus matrix
% Xcell: cell array of xmatrices used to train models (stimulus must be first element of cell array)
% fit_data: struct containing general data
%     tr_units: the set of units used to infer eye position
%     used_inds: the set of time-indices from the raw stim_mat used for model-fitting. 
% HMM_data: struct containing parameters for HMM eye-position inference, as well as precomputed shift-matrices
% ET_meas: information from measured eye position signals
% fixation_data: information needed to parse the data into separate fixations
%
% NOTES: stimulus must be first xtarget of models. any additional covariates must be associated with linear
% filters

add_usfac = gen_data.add_usfac;
[NT,n_units] = size(Robs_mat);
n_tr_units = length(tr_units); %number of units used to infer eye-pos

%assuming all models have the same structure!
mod_stim_params = init_mods(1).stim_params;
stim_klen = prod(mod_stim_params(1).dims); %visual stimulus must be the first Xtarg
flen = mod_stim_params(1).dims(1);
Xtargs = [init_mods(1).subunits(:).Xtarg];
stim_mod_signs = [init_mods(1).subunits(find(Xtargs == 1)).weight];
stim_NL_types = init_mods(1).get_NLtypes(find(Xtargs == 1));
n_stim_subs = sum(Xtargs == 1);
assert(strcmp(init_mods(1).spkNL.type,'softplus'),'only works with softplus spkNL at this point');

n_add_covariates = length(nonstim_X); %number of additional "Xmats" (beyond the visual stim)
fix_shifts = HMM_data.fix_shifts;
n_fix_shifts = length(fix_shifts);
drift_shifts = HMM_data.drift_shifts;
n_drift_shifts = length(drift_shifts);
drift_dsf = HMM_data.drift_dsf;

stim_dx = gen_data.stim_dx;

n_trials = length(gen_data.trial_boundaries);

use_coils = any(HMM_data.use_coils);
neural_delay = gen_data.neural_delay;
%%
base_Xmat = NIM.create_time_embedding(stim_mat,gen_data.stim_params);
base_Xmat = base_Xmat(gen_data.used_inds,:);

%%
%compute fixation IDS (normal and forward-projected)
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

fix_ids = nan(NT,1);
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_boundaries(ii,1):fix_boundaries(ii,2);
    fix_ids(cur_inds) = ii;
    cur_inds = pfix_boundaries(ii,1):pfix_boundaries(ii,2);
    pfix_ids(cur_inds) = ii;
end

%% ITERATE FIXATION-BASED CORRECTIONS

fix_it_mods{1} = init_mods;
fix_it_LLimp = nan(HMM_data.n_fix_iters+1,n_units);
for ss = 1:length(init_mods)
    fix_it_LLimp(1,ss) = (init_mods(ss).fit_props.LL - init_mods(ss).fit_props.null_LL)/log(2);
end

fix_it_post_mean = nan(HMM_data.n_fix_iters+1,n_fixs);
fix_it_post_std = nan(HMM_data.n_fix_iters+1,n_fixs);
for nn = 1:HMM_data.n_fix_iters
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,HMM_data.n_fix_iters);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(stim_klen/add_usfac,sum(Xtargs==1),n_tr_units); %compile all stimulus filters
    lin_kerns = cell(n_add_covariates,1);
    for ii = 1:n_add_covariates
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
    for ii = 1:n_add_covariates
        lin_pred_out = lin_pred_out + nonstim_X{ii}*lin_kerns{ii};
    end
    
    %% ESTIMATE LL for each shift in each stimulus frame
    %precompute LL at all shifts for all units
    frame_LLs = nan(NT,n_fix_shifts);
    reverseStr = '';
    for xx = 1:n_fix_shifts
        msg = sprintf('Calculating LLs for fixation shift %d of %d\n',xx,n_fix_shifts);
        fprintf([reverseStr msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        cur_stim_shift = base_Xmat*HMM_data.fix_shift_mats{xx};
        
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
        
        frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
    end
    
    %% INFER MICRO-SAC SEQUENCE    fix_LLs = nan(n_fixs,n_fix_shifts);
    fix_LLs = nan(n_fixs,n_fix_shifts);
    for ii = 1:n_fixs %compute total LL of each eye position for all data within each fixation
        cur_inds = pfix_boundaries(ii,1):pfix_boundaries(ii,2); %use forward-projected fixation indices to account for neural delay
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:)); %sum LLs for data within the fixation
    end
    
    if ~use_coils %if not using coil info we just compute the posterior at each fixation independently
        Lgamma = bsxfun(@plus,fix_LLs,HMM_data.fix_Lprior); %add log-likelihood and log-prior
        Lgamma = bsxfun(@minus,Lgamma,logsumexp(Lgamma,2)); %normalize
    else
        Lalpha=zeros(n_fixs,n_fix_shifts); %log of forward messages
        Lbeta = zeros(n_fixs,n_fix_shifts); %log of backward messages
        Lalpha(1,:) = HMM_data.fix_Lprior + fix_LLs(1,:); %initialize using prior
        for t=2:n_fixs %compute forward messages
            if ~gen_data.fix_post_blink(t) %if this fixation is NOT following a blink
                %make quadratic transition matrix centered on coil-measured deltaX, with variance
                %given by the delta_noise_sigma^2
                cdist = pdist2(fix_shifts'*stim_dx + ET_meas.fix_deltas(t),fix_shifts'*stim_dx);
                cur_LA = -cdist.^2/(2*HMM_data.fix_delta_noise_sigma^2);
                cur_LA = bsxfun(@plus,cur_LA,HMM_data.fix_Lprior); %multiply this by the prior over fixation positions
                cur_LA = bsxfun(@minus,cur_LA,logsumexp(cur_LA,2)); %normalize log transition-matrix
                Lalpha(t,:) = logmulexp(Lalpha(t-1,:),cur_LA) + fix_LLs(t,:); %compute forward message
            else %if it is following a blink, ignore the eye-trackers
                Lalpha(t,:) = HMM_data.fix_Lprior + fix_LLs(t,:);
            end
        end
        
        Lbeta(n_fixs,:)=zeros(1,n_fix_shifts); %initialize log-backward message
        for t=n_fixs-1:-1:1%compute backward messages
            if ~gen_data.fix_post_blink(t+1) %if the fixation isn't followed by a blink
                %construct A contingent on measured deltX between current and NEXT fixation
                cdist = pdist2(fix_shifts'*stim_dx + ET_meas.fix_deltas(t+1),fix_shifts'*stim_dx);
                cur_LA = -cdist.^2/(2*HMM_data.fix_delta_noise_sigma^2);
                cur_LA = bsxfun(@plus,cur_LA,HMM_data.fix_Lprior); %multiply by log-prior
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
    cur_fix_shifts = round(fix_post_mean_cor);
    shift_stim_mat = stim_mat;
    for ii=1:NT
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
        
        fix_it_mods{nn+1}(ss) = fix_it_mods{nn}(ss);
        fix_it_mods{nn+1}(ss).spkNL.params = [1 1]; %set alpha and beta to their default values when fitting filters to avoid issues with regularization
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

%%
if HMM_data.n_fix_iters == 0
    all_Xmat = NIM.create_time_embedding(stim_mat,gen_data.stim_params); %if we havent already made this
end
all_Xmat = all_Xmat(gen_data.used_inds,:); %grab relevant time indices
base_Xmat = all_Xmat;
%%
drift_it_mods{1} = fix_it_mods{end};
drift_it_LLimp = nan(HMM_data.n_drift_iters+1,n_units);
drift_it_LLimp(1,:) = fix_it_LLimp(end,:);

drift_it_post_mean = nan(HMM_data.n_drift_iters+1,NT);
drift_it_post_std = nan(HMM_data.n_drift_iters+1,NT);
for nn = 1:HMM_data.n_drift_iters
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,HMM_data.n_drift_iters);
    
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
        cur_stim_shift = base_Xmat*HMM_data.drift_shift_mats{xx};
        
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
        
        frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
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
            talpha=zeros(nt_pts,n_drift_shifts);
            tbeta = zeros(nt_pts,n_drift_shifts);
            
            tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf); %drift inference locations
            tpt_loc(end) = ntset;
            
            if use_coils; cur_drift_mean = ET_meas.post_drift_mean(tset); end;
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
                cur_LA = repmat(HMM_data.base_LA,[1 1 nt_pts]); %log-trans matrix is constant
            else
                cur_drift_mean = drift_dsf*nanmean(reshape(cur_drift_mean,[drift_dsf nt_pts])); %assume constant drift vel, so measured deltaX is scaled by drift_dsf
                cur_LA = nan(n_drift_shifts,n_drift_shifts,nt_pts);
                for iii = 1:nt_pts
                    if ~isnan(cur_drift_mean(iii))
                        %transition matrix given by gaussian posterior with mean cur_drift_mean and variance
                        %determined by post_drift_sigma (scaled by drift_dsf)
                        cdist = pdist2(drift_shifts'*stim_dx + cur_drift_mean(iii),drift_shifts'*stim_dx);
                        cur_LA(:,:,iii) = -cdist.^2/(2*(HMM_data.post_drift_sigma*drift_dsf)^2);
                    else
                        cur_LA(:,:,iii) = HMM_data.base_LA; %otherwise set to prior
                    end
                end
                cur_LA = bsxfun(@minus,cur_LA,logsumexp(cur_LA,2)); %normalize
            end
            
            %forward messages
            talpha(1,:) = HMM_data.drift_jump_Lprior + cur_LL_set(1,:); %initialize forward messages within this fixation
            for t = 2:nt_pts
                talpha(t,:) = logmulexp(talpha(t-1,:),cur_LA(:,:,t)) + cur_LL_set(t,:);
            end
            
            %backward messages
            tbeta(end,:)=log(ones(1,n_drift_shifts));
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
    
    %%
    if HMM_data.n_fix_iters == 0
        fix_post_mean_cor = zeros(NT,1);
    end
    cur_fix_shifts = round(fix_post_mean_cor + drift_post_cor);
    
    shift_stim_mat = stim_mat;
    for ii=1:NT
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
fit_post_data.mean = fix_it_post_mean;
fix_post_data.SD = fix_it_post_std;
drift_post_data.mean = drift_it_post_mean;
drift_post_data.SD = drift_it_post_std;
