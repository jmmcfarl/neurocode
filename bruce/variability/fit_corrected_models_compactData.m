clear all
close all

global Expt_name bar_ori monk_name rec_type fit_unCor use_MUA

Expt_name = 'M005';
monk_name = 'jbe';
bar_ori = 50; %bar orientation to use (only for UA recs)
fit_unCor = false; %use eye correction
use_MUA = false; %use MUA in model-fitting

save_name = 'corrected_models_comp';

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

%load in packaged data
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
fprintf('Loading %s\n',data_name);
load(data_name);


%%
Expt_num = str2num(Expt_name(2:end));

%load in array RF position data
load ~/Data/bruce/general_array_data/array_pos_data.mat
interp_ecc = sqrt(interp_x.^2 + interp_y.^2);

cd(data_dir);
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
save_dir = ['~/Analysis/bruce/' Expt_name '/models'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';

%if using coil info
if any(params.use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
save_name = [save_name sprintf('_ori%d',bar_ori)];

%% LOAD EYE-TRACKING DATA
cd(et_dir)
load(et_mod_data_name,'all_mod*');
load(et_anal_name,'drift*','it_*','et_tr_set','et_saccades','et_params');
tr_set = et_tr_set;

%%
if ~isnan(params.rpt_seeds)
    xv_type = 'rpt';
else
    xv_type = 'uni';
end
xv_frac = 0.2; %fraction of trials to use for XVAL

flen = 15; %time lags for ST filters
spatial_usfac = et_params.spatial_usfac; %spatial up-sampling factor

use_nPix = et_params.use_nPix; %number of pixels (bars) used in models
full_nPix = params.full_nPix; %total number of bars to keep track of in stimulus
dt = params.dt;

use_nPix_us = use_nPix*spatial_usfac; %number of pixels in stimulus filters
klen_us = use_nPix_us*flen; %number of parameters in stim-filters
sp_dx = et_params.sp_dx; %pixel size (deg)

use_LOOXV = 1; %[0 is no LOO; 1 is SUs only; 2 is SU + MU]
all_stim_mat = decompressTernNoise(stimComp);

%% spatial up-sampling of stimulus
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

%select subset of pixels used for model fitting
buffer_pix = floor((full_nPix - use_nPix)/2);
[Xinds_up,~] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
cur_use_pix = (1/spatial_usfac:1/spatial_usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

stim_params_full = NMMcreate_stim_params([flen full_nPix_us]);

%% BIN SPIKES FOR MU AND SU
all_binned_mua = spikes_int82double(spike_data.binned_mua);
all_binned_sua = spikes_int82double(spike_data.binned_sua);
Clust_data = spike_data.Clust_data;
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

%% extract time series to keep track of trial number and block number
NT = length(used_inds);
fullNT = size(all_binned_mua,1);
n_trials = length(time_data.trial_flip_ids);
n_blocks = length(expt_data.used_blocks);

all_t_axis = time_data.t_axis;
trial_start_inds = [1+time_data.trial_flip_inds];
trial_end_inds = [time_data.trial_flip_inds(2:end); fullNT];
all_trialvec = nan(fullNT,1); %trial index vector
for ii = 1:n_trials
    all_trialvec(trial_start_inds(ii):trial_end_inds(ii)) = time_data.trial_flip_ids(ii);
end

block_start_inds = [1+time_data.block_flip_inds];
block_end_inds = [time_data.block_flip_inds(2:end); fullNT];
all_blockvec = nan(fullNT,1); %block index vector
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

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

%% DEFINE FIXATIONS
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink));
n_fixs = length(fix_start_inds);

%index values numbering the different fixations
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

%% Recon retinal stim for non LOO data
sac_shift = et_params.sac_shift; %forward projection of saccade start times

cur_fix_post_mean = squeeze(it_fix_post_mean(end,:));
cur_fix_post_std = squeeze(it_fix_post_std(end,:));
cur_drift_post_mean = squeeze(drift_post_mean(end,:));
cur_drift_post_std = squeeze(drift_post_std(end,:));
[fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
    cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);

fin_shift_cor = round(fin_tot_corr);

%RECOMPUTE XMAT
best_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
end
all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_full);
all_Xmat_shift = all_Xmat_shift(used_inds,use_kInds_up);

if fit_unCor
    all_Xmat_unCor = create_time_embedding(all_stimmat_up,stim_params_full);
    all_Xmat_unCor = all_Xmat_unCor(used_inds,use_kInds_up);
end

%% COMBINE SUA AND MUA and make a matrix out of the binned spike data
n_units = size(all_binned_mua,2) + size(all_binned_sua,2);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_units);
for ss = 1:n_units
    if ss > params.n_probes %these are the SUs
        su_probe_ind = ss-params.n_probes;
        Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end

%% IDENTIFY REPEAT TRIALS
%don't use repeat trials for training or cross-validation here.
if strcmp(xv_type,'rpt')
    rpt_trials = find(ismember([trial_data(:).se],params.rpt_seeds));
    rpt_inds = find(ismember(all_trialvec(used_inds),rpt_trials));
end

%%
%set of units where we have LOOXV on eye-tracking
if use_LOOXV == 2
    loo_set = tr_set; %all usable units
elseif use_LOOXV == 1
    loo_set = tr_set(all_mod_SU(tr_set) > 0); % just SUs
else
    loo_set = [];
end

%for model fitting
if use_MUA
    targs = 1:n_units; %SU and MU
else
    targs = setdiff(1:n_units,1:params.n_probes); %SU only
end
%%
for cc = targs
    fprintf('Starting model fits for unit %d\n',cc);
    
    loo_cc = find(loo_set == cc); %index within the LOOXV set
    cur_Robs = Robs_mat(:,cc);
    cc_uinds = find(~isnan(cur_Robs)); %usable time points for this unit
        
    if ~isempty(cc_uinds)
        use_trials = unique(all_trialvec(used_inds(cc_uinds))); %set of potentially usable trials
        use_trials(ismember(use_trials,rpt_trials)) = []; %dont use repeat trials
        
        %create sets of training and XV trials
        nuse_trials = length(use_trials);
        n_xv_trials = round(xv_frac*nuse_trials);
        xv_trials = randperm(nuse_trials);
        xv_trials(n_xv_trials+1:end) = [];
        xv_trials = use_trials(xv_trials);
        tr_trials = setdiff(use_trials,xv_trials);
        n_tr_trials = length(tr_trials);
        fprintf('Initializing models with %d training trials and %d xval trials\n',n_tr_trials,n_xv_trials);
        
        cur_tr_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),tr_trials));
        cur_xv_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),xv_trials));
        cur_full_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),use_trials));
        
        %% COMPUTE UNIT DATA
        unit_data.isLOO = ismember(cc,loo_set);
        unit_data.avg_rate = mean(cur_Robs)/dt;
        unit_data.tot_spikes = sum(cur_Robs);
        unit_data.N_used_samps = length(cur_Robs);
        used_blocks = unique(all_blockvec(used_inds(cc_uinds)));
        unit_data.n_used_blocks = length(used_blocks);
        block_rates = nan(unit_data.n_used_blocks,1);
        for ii = 1:unit_data.n_used_blocks
            block_rates(ii) = mean(cur_Robs(all_blockvec(used_inds(cc_uinds)) == used_blocks(ii)));
        end
        unit_data.rate_stability_cv = std(block_rates)/mean(block_rates);
        unit_data.block_rates = block_rates/dt;
        if cc > params.n_probes
            su_Ind = cc - params.n_probes;
            unit_data.SU_number = spike_data.Clust_data.SU_numbers(su_Ind);
            unit_data.probe_number = spike_data.Clust_data.SU_probes(su_Ind);
            unit_data.SU_Lratio = nanmean(spike_data.Clust_data.SU_Lratios(su_Ind,used_blocks),2);
            unit_data.SU_isodist = nanmean(spike_data.Clust_data.SU_isodists(su_Ind,used_blocks(spike_data.Clust_data.SU_isoRel(su_Ind,used_blocks)==1)),2);
            unit_data.SU_refract = nanmean(spike_data.Clust_data.SU_refract(su_Ind,used_blocks),2);
            unit_data.SU_dprime = nanmean(spike_data.Clust_data.SU_dprimes(su_Ind,used_blocks),2);
            unit_data.dprime_stability_cv = nanstd(spike_data.Clust_data.SU_dprimes(su_Ind,used_blocks),[],2)/nanmean(spike_data.Clust_data.SU_dprimes(su_Ind,used_blocks),2);
        else
            unit_data.SU_number = nan;
            unit_data.probe_number = cc;
            unit_data.SU_Lratio = nan;
            unit_data.SU_isodist = nan;
            unit_data.SU_refract = nan;
            unit_data.SU_dprime = nan;
        end
        
        %% RECON RETINAL STIM
        if ~isempty(loo_cc)
            fprintf('Reconstructing retinal stim for unit %d\n',cc);
            cur_fix_post_mean = squeeze(it_fix_post_mean_LOO(loo_cc,end,:));
            cur_fix_post_std = squeeze(it_fix_post_std_LOO(loo_cc,end,:));
            cur_drift_post_mean = squeeze(drift_post_mean_LOO(loo_cc,end,:));
            cur_drift_post_std = squeeze(drift_post_std_LOO(loo_cc,end,:));
            [fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
                cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
            
            fin_shift_cor = round(fin_tot_corr);
            
            %RECOMPUTE XMAT
            all_shift_stimmat_up = all_stimmat_up;
            for i=1:NT
                all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
            end
            all_Xmat_shift = create_time_embedding(all_shift_stimmat_up,stim_params_full);
            all_Xmat_shift = all_Xmat_shift(used_inds,use_kInds_up);
        end
        
        %% FIT (eye-corrected) STIM-PROCESSING MODEL
        fprintf('Fitting stim models for unit %d\n',cc);
        silent = 1;
        if mode(expt_data.expt_dds) == 12
            base_lambda_d2XT = 25; %target ST smoothness strength
            base_lambda_L1 = 2.5; %sparseness strength
            init_lambda_d2XT = 2.5; %initial ST smoothness strength
        else
            base_lambda_d2XT = 100;
            base_lambda_L1 = 10;
            init_lambda_d2XT = 10;
        end
        init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_lambda_d2XT,'boundary_conds',[0 0 0]);
        mod_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
        
        if cc <= n_probes %for MUA just fit a 3-filter quad model
            
            %start out with just a linear filter
            nEfilts = 2;
            nIfilts = 0;
            
            stim_mod_signs = [1 ones(1,nEfilts) -1*ones(1,nIfilts)];
            stim_NL_types = [{'lin'} repmat({'quad'},1,nEfilts) repmat({'quad'},1,nIfilts)];
            
            init_mod = NMMinitialize_model( mod_stim_params, stim_mod_signs, stim_NL_types, init_reg_params);
            init_mod = NMMfit_filters(init_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds,silent);
            [~, ~, ~, ~, gint] = NMMeval_model(init_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds);
            vgint = var(gint);sgint = std(gint);
            init_mod = NMMadjust_regularization(init_mod,[],'lambda_d2XT',base_lambda_d2XT./vgint);
            init_mod = NMMadjust_regularization(init_mod,[],'lambda_L1',base_lambda_L1./sgint);
            
            bestGQM = NMMfit_filters(init_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds,silent);
            bestGQM_spkNL = NMMfit_logexp_spkNL(bestGQM,cur_Robs,all_Xmat_shift,[],cur_tr_inds);
            if isempty(cur_xv_inds)
                error('Need xval trials');
            end
            [bestGQM_xvLL, nullxvLL] = NMMeval_model(bestGQM_spkNL, cur_Robs, all_Xmat_shift,[],cur_xv_inds);
            bestGQM.xvLLimp = (bestGQM_xvLL-nullxvLL)/log(2);
            
        else %for SUA try a number of different model structures and pick the best by xvalLL
            
            max_Emods = 3; %max excitatory quad filters
            max_Imods = 3; %max inhibitory quad filters
            
            %start out with just a linear filter
            nEfilts = 0;
            nIfilts = 0;
            
            stim_mod_signs = [1 ones(1,nEfilts) -1*ones(1,nIfilts)];
            stim_NL_types = [{'lin'} repmat({'quad'},1,nEfilts) repmat({'quad'},1,nIfilts)];
            
            init_mod = NMMinitialize_model( mod_stim_params, stim_mod_signs, stim_NL_types, init_reg_params);
            init_mod = NMMfit_filters(init_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds,silent);
            [~, ~, ~, ~, gint] = NMMeval_model(init_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds);
            vgint = var(gint);sgint = std(gint);
            init_mod = NMMadjust_regularization(init_mod,[],'lambda_d2XT',base_lambda_d2XT./vgint);
            init_mod = NMMadjust_regularization(init_mod,[],'lambda_L1',base_lambda_L1./sgint);
            
            bestGQM = NMMfit_filters(init_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds,silent);
            bestGQM_spkNL = NMMfit_logexp_spkNL(bestGQM,cur_Robs,all_Xmat_shift,[],cur_tr_inds);
            if isempty(cur_xv_inds)
                error('Need xval trials');
            end
            [bestGQM_xvLL, nullxvLL] = NMMeval_model(bestGQM_spkNL, cur_Robs, all_Xmat_shift,[],cur_xv_inds);
            
            %now try adding up to max_Emods exc squared filters. Stop when xval
            %LL doesn't improve
            cur_imp = Inf;
            while nEfilts < max_Emods && cur_imp > 0
                cur_mod = bestGQM;
                cur_mod = NMMadd_NLinput(cur_mod,'quad',1,1,0.1*randn(flen*use_nPix_us,1));
                cur_mod.mods(end).reg_params = init_reg_params;
                nEfilts = nEfilts + 1;
                
                fprintf('Fitting model with LIN + %dE and %dI\n',nEfilts,nIfilts);
                cur_mod = NMMfit_filters(cur_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds,silent);
                [~, ~, ~, ~, gint] = NMMeval_model(cur_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds);
                vgint = var(gint); sgint = std(gint);
                cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_d2XT',base_lambda_d2XT./vgint(end));
                cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_L1',base_lambda_L1./sgint(end));
                cur_mod = NMMfit_filters(cur_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds,silent);
                cur_mod_spkNL = NMMfit_logexp_spkNL(cur_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds);
                cur_mod_xvLL = NMMeval_model(cur_mod_spkNL, cur_Robs, all_Xmat_shift,[],cur_xv_inds);
                
                cur_imp = cur_mod_xvLL - bestGQM_xvLL;
                if cur_imp > 0
                    fprintf('Keeping new model\n');
                    bestGQM = cur_mod;
                    bestGQM_spkNL = cur_mod_spkNL;
                    bestGQM_xvLL = cur_mod_xvLL;
                end
            end
            
            %now try adding up to max_Imods inh squared filters. Stop when xval
            %LL doesn't improve
            cur_imp = Inf;
            while nIfilts < max_Imods && cur_imp > 0
                cur_mod = bestGQM;
                cur_mod = NMMadd_NLinput(cur_mod,'quad',-1,1,0.1*randn(flen*use_nPix_us,1));
                cur_mod.mods(end).reg_params = init_reg_params;
                nIfilts = nIfilts + 1;
                
                fprintf('Fitting model with LIN + %dE and %dI\n',nEfilts,nIfilts);
                cur_mod = NMMfit_filters(cur_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:),[],[],silent);
                [~, ~, ~, ~, gint] = NMMeval_model(cur_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds);
                vgint = var(gint); sgint = std(gint);
                cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_d2XT',base_lambda_d2XT./vgint(end));
                cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_L1',base_lambda_L1./sgint(end));
                cur_mod = NMMfit_filters(cur_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds,silent);
                cur_mod_spkNL = NMMfit_logexp_spkNL(cur_mod,cur_Robs,all_Xmat_shift,[],cur_tr_inds);
                cur_mod_xvLL = NMMeval_model(cur_mod_spkNL, cur_Robs, all_Xmat_shift,[],cur_xv_inds);
                
                cur_imp = cur_mod_xvLL - bestGQM_xvLL;
                if cur_imp > 0
                    fprintf('Keeping new model\n');
                    bestGQM = cur_mod;
                    bestGQM_spkNL = cur_mod_spkNL;
                    bestGQM_xvLL = cur_mod_xvLL;
                end
            end
            
            bestGQM.xvLLimp = (bestGQM_xvLL-nullxvLL)/log(2);
        end
        %% Fit a GQM with the linear filter split into two rectified subunits
        %convert linear filter into separate thresh-lin filters and refit
        rectGQM = bestGQM;
        rectGQM.mods = [rectGQM.mods(:); rectGQM.mods(1)];
        rectGQM.mods(end).sign = -1;
        rectGQM.mods(end).filtK = -rectGQM.mods(end).filtK; %initialize filter to be negative of the linear filter
        rectGQM.mods(1).NLtype = 'threshlin'; rectGQM.mods(end).NLtype = 'threshlin';
        rectGQM = NMMfit_filters(rectGQM,cur_Robs,all_Xmat_shift,[],cur_tr_inds,silent);
        rectGQM_spkNL = NMMfit_logexp_spkNL(rectGQM,cur_Robs,all_Xmat_shift,[],cur_tr_inds);
        rectGQM_xvLL = NMMeval_model(rectGQM_spkNL, cur_Robs, all_Xmat_shift,[],cur_xv_inds);
        rectGQM.xvLLimp = (rectGQM_xvLL-nullxvLL)/log(2);
        
        %% Fit uncorrected model
        if fit_unCor
            bestGQM_unCor = bestGQM;
            bestGQM_unCor = NMMfit_filters(bestGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_tr_inds);
            cur_mod_spkNL = NMMfit_logexp_spkNL(bestGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_tr_inds);
            cur_xvLL = NMMeval_model(cur_mod_spkNL,cur_Robs,all_Xmat_unCor,[],cur_xv_inds);
            bestGQM_unCor.xvLLimp = (cur_xvLL - nullxvLL)/log(2);
            
            rectGQM_unCor = rectGQM;
            rectGQM_unCor = NMMfit_filters(rectGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_tr_inds);
            cur_mod_spkNL = NMMfit_logexp_spkNL(rectGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_tr_inds);
            cur_xvLL = NMMeval_model(cur_mod_spkNL,cur_Robs,all_Xmat_unCor,[],cur_xv_inds);
            rectGQM_unCor.xvLLimp = (cur_xvLL - nullxvLL)/log(2);
        end
        
        %% Refit model parameters on the full data set (excluding repeats)
        if ~isempty(cur_xv_inds)
            bestGQM = NMMfit_filters(bestGQM,cur_Robs,all_Xmat_shift,[],cur_full_inds,silent);
            bestGQM_spkNL = NMMfit_logexp_spkNL(bestGQM,cur_Robs,all_Xmat_shift,[],cur_full_inds);
            [bestGQM_LL, nullLL, ~, ~, ~, fgint] = NMMeval_model(bestGQM_spkNL, cur_Robs, all_Xmat_shift,[],cur_full_inds);
            rel_filt_weights = std(fgint);
            bestGQM.rel_filt_weights = rel_filt_weights/sum(rel_filt_weights);
            bestGQM.LLimp = (bestGQM_LL - nullLL)/log(2);
            
            if fit_unCor
                bestGQM_unCor = NMMfit_filters(bestGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_full_inds);
                bestGQM_unCor_spkNL = NMMfit_logexp_spkNL(bestGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_full_inds);
                [~,~,~,~,~,fgint] = NMMeval_model(bestGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_full_inds);
                rel_filt_weights = std(fgint);
                bestGQM_unCor.rel_filt_weights = rel_filt_weights/sum(rel_filt_weights);
                bestGQM_unCor.LLimp = (bestGQM_unCor_spkNL.LL_seq(end)-nullLL)/log(2);
            end
            rectGQM = NMMfit_filters(rectGQM,cur_Robs,all_Xmat_shift,[],cur_full_inds,silent);
            [~,~,~,~,~,fgint] = NMMeval_model(rectGQM,cur_Robs,all_Xmat_shift,[],cur_full_inds);
            rel_filt_weights = std(fgint);
            rectGQM.rel_filt_weights = rel_filt_weights/sum(rel_filt_weights);
            rectGQM_spkNL = NMMfit_logexp_spkNL(rectGQM,cur_Robs,all_Xmat_shift,[],cur_full_inds);
            rectGQM.LLimp = (rectGQM_spkNL.LL_seq(end)-nullLL)/log(2);
            
            if fit_unCor
                rectGQM_unCor = NMMfit_filters(rectGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_full_inds);
                [~,~,~,~,~,fgint] = NMMeval_model(rectGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_full_inds);
                rel_filt_weights = std(fgint);
                rectGQM_unCor.rel_filt_weights = rel_filt_weights/sum(rel_filt_weights);
                rectGQM_unCor_spkNL = NMMfit_logexp_spkNL(rectGQM_unCor,cur_Robs,all_Xmat_unCor,[],cur_full_inds);
                rectGQM_unCor.LLimp = (rectGQM_unCor_spkNL.LL_seq(end)-nullLL)/log(2);
            end
        end
        
        %% Store data
        fullLL = nansum(cur_Robs(cur_full_inds).*log(cur_Robs(cur_full_inds)) - cur_Robs(cur_full_inds));
        fullxvLL = nansum(cur_Robs(cur_xv_inds).*log(cur_Robs(cur_xv_inds)) - cur_Robs(cur_xv_inds));
        unit_data.fullLL = fullLL;
        unit_data.nullLL = nullLL*nansum(cur_Robs(cur_full_inds));
        unit_data.fullxvLL = fullxvLL;
        unit_data.nullxvLL = nullxvLL*nansum(cur_Robs(cur_xv_inds));
        unit_data.tr_trials = tr_trials;
        unit_data.xv_trials = xv_trials;
        
        ModData(cc).unit_data = unit_data;
        ModData(cc).bestGQM = bestGQM;
        if fit_unCor
            ModData(cc).bestGQM_unCor = bestGQM_unCor;
        end
        ModData(cc).rectGQM = rectGQM;
        if fit_unCor
            ModData(cc).rectGQM_unCor = rectGQM_unCor;
        end
        
        %% CALCULATE TUNING PROPERTIES
        clear tune_props
        rectGQM = ModData(cc).rectGQM;
        
        stim_filters = [rectGQM.mods(:).filtK];
        stim_mod_signs = [rectGQM.mods(:).sign];
        mod_stim_params = rectGQM.stim_params;
        filt_data = get_filter_properties_v2(stim_filters,mod_stim_params,sp_dx);
        tune_props.filt_data = filt_data;
        
        [~,~,best_pred_rate,~,~,filt_outs] = NMMeval_model(rectGQM,cur_Robs,all_Xmat_shift,[],cur_full_inds);
        [~,~,rev_pred_rate] = NMMeval_model(rectGQM,cur_Robs,-all_Xmat_shift,[],cur_full_inds);
        tune_props.PRM = mean(abs(best_pred_rate - rev_pred_rate))/mean(best_pred_rate);
        
        %only use E/lin filters for these calculations
        non_supp_filts = find(stim_mod_signs ~= -1);
        filt_out_weights = std(filt_outs);
        filt_out_weights = filt_out_weights(non_supp_filts)/nansum(filt_out_weights(non_supp_filts));
        
        tune_props.RF_mean = nansum(filt_out_weights(non_supp_filts).*filt_data.gest(non_supp_filts,1)') - use_nPix_us*sp_dx/2 -sp_dx;
        tune_props.RF_sigma = nansum(filt_out_weights(non_supp_filts).*filt_data.gest(non_supp_filts,2)');
        tune_props.RF_gSF = nansum(filt_out_weights(non_supp_filts).*filt_data.gest(non_supp_filts,3)');
        tune_props.RF_dirsel = nansum(filt_out_weights(non_supp_filts).*filt_data.dir_selectivity(non_supp_filts)');
        tune_props.RF_FTF = nansum(filt_out_weights(non_supp_filts).*filt_data.FFt(non_supp_filts)');
        tune_props.RF_FSF = nansum(filt_out_weights(non_supp_filts).*filt_data.FFx(non_supp_filts)');
        
        rect_stim_filters = [rectGQM.mods([1 end]).filtK];
        avg_rect_filts = mean(rect_stim_filters);
        absavg_rect_filts = mean(abs(rect_stim_filters(:)));
        tune_props.net_phase_polarity = (avg_rect_filts(1) - avg_rect_filts(end))/absavg_rect_filts;
        
        screen_X =  -tune_props.RF_mean'.*sind(bar_ori) + ov_RF_pos(1);
        screen_Y =  tune_props.RF_mean'.*cosd(bar_ori) + ov_RF_pos(2);
        if strcmp(rec_type,'UA')
            if bar_ori == 0
                screen_X = interp_x(unit_data.probe_number);
            elseif bar_ori == 90
                screen_Y = interp_y(unit_data.probe_number);
            else
                error('Code doesnt work with this condition!');
            end
        end
        
        tune_props.screen_X = screen_X;
        tune_props.screen_Y = screen_Y;
        tune_props.RF_ecc = sqrt(screen_X.^2 + screen_Y.^2);
        
        ModData(cc).tune_props = tune_props;
        
        %% CALCULATE UNCORRECTED TUNING PROPERTIES
        clear tune_props_unCor
        if fit_unCor
            rectGQM_unCor = ModData(cc).rectGQM_unCor;
            
            stim_filters = [rectGQM_unCor.mods(:).filtK];
            stim_mod_signs = [rectGQM_unCor.mods(:).sign];
            mod_stim_params = rectGQM_unCor.stim_params;
            filt_data = get_filter_properties_v2(stim_filters,mod_stim_params,sp_dx);
            tune_props_unCor.filt_data = filt_data;
            
            [~,~,best_pred_rate,~,~,filt_outs] = NMMeval_model(rectGQM_unCor,cur_Robs,all_Xmat_shift,[],cur_full_inds);
            [~,~,rev_pred_rate] = NMMeval_model(rectGQM_unCor,cur_Robs,-all_Xmat_shift,[],cur_full_inds);
            tune_props_unCor.PRM = mean(abs(best_pred_rate - rev_pred_rate))/mean(best_pred_rate);
            
            %only use E/lin filters for these calculations
            non_supp_filts = find(stim_mod_signs ~= -1);
            filt_out_weights = std(filt_outs);
            filt_out_weights = filt_out_weights(non_supp_filts)/nansum(filt_out_weights(non_supp_filts));
            
            tune_props_unCor.RF_mean = filt_out_weights(non_supp_filts)*filt_data.gest(non_supp_filts,1) - use_nPix_us*sp_dx/2 -sp_dx;
            tune_props_unCor.RF_sigma = filt_out_weights(non_supp_filts)*filt_data.gest(non_supp_filts,2);
            tune_props_unCor.RF_gSF = filt_out_weights(non_supp_filts)*filt_data.gest(non_supp_filts,3);
            tune_props_unCor.RF_dirsel = filt_out_weights(non_supp_filts)*filt_data.dir_selectivity(non_supp_filts);
            tune_props_unCor.RF_FTF = filt_out_weights(non_supp_filts)*filt_data.FFt(non_supp_filts);
            tune_props_unCor.RF_FSF = filt_out_weights(non_supp_filts)*filt_data.FFx(non_supp_filts);
            
            rect_stim_filters = [rectGQM_unCor.mods([1 end]).filtK];
            avg_rect_filts = mean(rect_stim_filters);
            absavg_rect_filts = mean(abs(rect_stim_filters(:)));
            tune_props_unCor.net_phase_polarity = (avg_rect_filts(1) - avg_rect_filts(end))/absavg_rect_filts;
            
            screen_X =  -tune_props_unCor.RF_mean'.*sind(bar_ori) + ov_RF_pos(1);
            screen_Y =  tune_props_unCor.RF_mean'.*cosd(bar_ori) + ov_RF_pos(2);
            if strcmp(rec_type,'UA')
                if bar_ori == 0
                    screen_X = interp_x(unit_data.probe_number);
                elseif bar_ori == 90
                    screen_Y = interp_y(unit_data.probe_number);
                else
                    error('Code doesnt work with this condition!');
                end
            end
            
            tune_props_unCor.screen_X = screen_X;
            tune_props_unCor.screen_Y = screen_Y;
            tune_props_unCor.RF_ecc = sqrt(screen_X.^2 + screen_Y.^2);
            
            ModData(cc).tune_props_unCor = tune_props_unCor;
        end
    end
end

%%
modFitParams = struct('bar_ori',bar_ori,'use_MUA',use_MUA,'fit_uncor',fit_unCor,...
    'xv_frac',xv_frac,'xv_type',xv_type,'use_nPix_us',use_nPix_us,'flen',flen,'spatial_usfac',spatial_usfac,'sp_dx',sp_dx,'dt',dt,...
    'base_lambda_d2XT',base_lambda_d2XT,'base_lambda_L1',base_lambda_L1,...
    'init_lambda_d2XT',init_lambda_d2XT);

cd(save_dir)
save(save_name,'ModData','modFitParams');