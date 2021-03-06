clear all
close all
% 
global Expt_name bar_ori monk_name rec_type rec_number

Expt_name = 'M320';
monk_name = 'lem';
bar_ori = 100; %bar orientation to use (only for UA recs)
rec_number = 1;

fit_unCor = false; %also fit models without eye corrections?
use_MUA = false; %use MUA in model-fitting
fit_rect = false; %split quad linear filter into two rectified

save_name = 'corrected_models_comp_PARDRIFT';
if fit_unCor
    save_name = strcat(save_name,'_unCor');
end

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
ff = find(cellfun(@(x) ~isempty(x),Expts),1);
if strcmp(Expts{ff}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{ff}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

%load in packaged data
data_name = sprintf('%s/PARDRIFT_packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
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
et_mod_data_name = 'PARDRIFT_eyetrack_initmods_Rinit';
et_anal_name = 'PARDRIFT_eyetrack_Rinit';

%if using coil info
if any(params.use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
save_name = [save_name sprintf('_ori%d',bar_ori)];

if rec_number > 1
    et_mod_data_name = strcat(et_mod_data_name,sprintf('r%d',rec_number));
    et_anal_name = strcat(et_anal_name,sprintf('r%d',rec_number));
    save_name = strcat(save_name,sprintf('_r%d',rec_number));
end

%% LOAD EYE-TRACKING DATA
cd(et_dir)
load(et_mod_data_name,'all_mod*');
load(et_anal_name,'drift*','it_*','dit_*','et_tr_set','et_saccades','et_params');
tr_set = et_tr_set;

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
if ~isnan(params.rpt_seeds)
    xv_type = 'rpt';
else
    xv_type = 'uni';
end
xv_frac = 0.2; %fraction of trials to use for XVAL

flen = 15; %time lags for ST filters
spatial_usfac = et_params.spatial_usfac;%spatial up-sampling factor
add_usfac = et_params.add_usfac; %additional (tent-basis) spatial up-sampling
base_usfac = spatial_usfac/add_usfac; %number of pixels used in model fits

use_nPix = et_params.use_nPix; %number of (bars) used in models
full_nPix = params.full_nPix; %total number of bars to keep track of in stimulus
dt = params.dt;

use_nPix_us = use_nPix*spatial_usfac; %number of pixels in stimulus filters
klen_us = use_nPix_us*flen; %number of parameters in stim-filters
sp_dx = et_params.sp_dx; %pixel size (deg)
mod_dx = sp_dx*add_usfac; %model pixel size

use_LOOXV = 0; %[0 is no LOO; 1 is SUs only; 2 is SU + MU]
all_stim_mat = decompressTernNoise(stimComp);

%RF center position on screen
ov_RF_pos = Expts{expt_data.used_blocks(1)}.Stimvals.rf(1:2)/params.scale_fac;
fix_point = [Expts{expt_data.used_blocks(1)}.Stimvals.fx Expts{expt_data.used_blocks(1)}.Stimvals.fy];
ov_RF_pos = ov_RF_pos - fix_point; %account for non-zero fixation point location

%% model-fitting parameters
% poss_smoothreg_scalefacs = logspace(-2,2,10); %possible scale factors to apply to smoothness reg strength
poss_smoothreg_scalefacs = 1; %possible scale factors to apply to smoothness reg strength

base_lambda_d2XT = 100;
base_lambda_L1 = 5;
init_lambda_d2XT = 100;

%number of quad filters to use for MUA
MUA_nEfilts = 2;
MUA_nIfilts = 0;

%max number of quad filter to use for fitting SUA
max_Emods = 4; %max excitatory quad filters
max_Imods = 4; %max inhibitory quad filters
           
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

%param struct for creating time-embedded stimulus
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
n_blocks = length(time_data.block_flip_ids);

all_t_axis = time_data.t_axis;
trial_start_inds = [1; 1+time_data.trial_flip_inds(2:end)];
trial_end_inds = [time_data.trial_flip_inds(2:end); fullNT];
all_trialvec = nan(fullNT,1); %trial index vector
for ii = 1:n_trials
    all_trialvec(trial_start_inds(ii):trial_end_inds(ii)) = time_data.trial_flip_ids(ii);
end

block_start_inds = [1; 1+time_data.block_flip_inds(2:end)];
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

% cur_fix_post_mean = squeeze(it_fix_post_mean(end,:));
% cur_fix_post_std = squeeze(it_fix_post_std(end,:));
% cur_drift_post_mean = squeeze(drift_post_mean(end,:));
% cur_drift_post_std = squeeze(drift_post_std(end,:));
% [fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
%     cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
% 
% fin_shift_cor = round(fin_tot_corr); %round to pixel shifts
% 
% %RECOMPUTE XMAT
% if ~fit_unCor
%     best_shift_stimmat_up = all_stimmat_up;
%     for i=1:NT
%         best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
%     end
%     all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_full);
% else
%     all_Xmat_shift = create_time_embedding(all_stimmat_up,stim_params_full);
% end

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
else
    rpt_trials = [];
    rpt_inds = [];
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
mod_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us/add_usfac],dt);
mod_stim_params(2) = NMMcreate_stim_params([n_blocks],dt);
% Xmat{1} = all_Xmat_shift(used_inds,use_kInds_up);
Xmat{2} = Xblock(used_inds,:);

%%
% cd(save_dir)
% load(save_name);
targs = [11 16]+params.n_probes;
%%
silent = 1;

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
        
        cur_tr_inds = cc_uinds(ismember(all_trialvec(used_inds(cc_uinds)),tr_trials));
        cur_xv_inds = cc_uinds(ismember(all_trialvec(used_inds(cc_uinds)),xv_trials));
        cur_full_inds = union(cur_tr_inds,cur_xv_inds);
        cur_rpt_inds = rpt_inds(ismember(rpt_inds,cc_uinds));
        
        %% COMPUTE UNIT DATA
        unit_data.isLOO = ismember(cc,loo_set);
        unit_data.avg_rate = nanmean(cur_Robs)/dt;
        unit_data.tot_spikes = nansum(cur_Robs);
        unit_data.N_used_samps = length(cur_Robs);
        used_blocks = unique(all_blockvec(used_inds(cc_uinds)));
        unit_data.n_used_blocks = length(used_blocks);
        block_rates = nan(unit_data.n_used_blocks,1);
        for ii = 1:unit_data.n_used_blocks
            block_rates(ii) = nanmean(cur_Robs(all_blockvec(used_inds(cc_uinds)) == used_blocks(ii)));
        end
        unit_data.rate_stability_cv = std(block_rates)/mean(block_rates);
        unit_data.block_rates = block_rates/dt;
        unit_data.used_blocks = used_blocks;
        unit_data.poss_used_blocks = unique(all_blockvec(used_inds));
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
        fprintf('Reconstructing retinal stim for unit %d\n',cc);
        if use_LOOXV
        cur_fix_post_mean = squeeze(it_fix_post_mean_LOO(loo_cc,end,:));
        cur_fix_post_std = squeeze(it_fix_post_std_LOO(loo_cc,end,:));
        cur_drift_post_mean = squeeze(drift_post_mean_LOO(loo_cc,end,:));
        cur_drift_post_std = squeeze(drift_post_std_LOO(loo_cc,end,:));
        else
        cur_fix_post_mean = squeeze(it_fix_post_mean(end,:));
        cur_fix_post_std = squeeze(it_fix_post_std(end,:));
        cur_drift_post_mean = squeeze(drift_post_mean(end,:));
        cur_drift_post_std = squeeze(drift_post_std(end,:));
        end
        [fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
            cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
        
        fin_shift_cor = round(fin_tot_corr);
        
        %RECOMPUTE XMAT
        all_shift_stimmat_up = all_stimmat_up;
        for i=1:NT
            all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
        end
        all_Xmat_shift = create_time_embedding(all_shift_stimmat_up,stim_params_full);
        
        %create X matrix
        if add_usfac > 1
            %if doing additional spatial up-sampling use tent basis functions
            Xmat{1} = tb_proc_stim(all_Xmat_shift(used_inds,use_kInds_up),add_usfac,flen);
        else
            Xmat{1} = all_Xmat_shift(used_inds,use_kInds_up);
        end
                
        %% compute null model
        nullMod = NMMinitialize_model(mod_stim_params,1,{'lin'},[],2);
        nullMod = NMMfit_filters(nullMod,cur_Robs,Xmat,[],cur_tr_inds,silent);

        nullMod_full = NMMfit_filters(nullMod,cur_Robs,Xmat,[],cur_full_inds,silent);
        nullMod = NMMfit_logexp_spkNL(nullMod,cur_Robs,Xmat,[],cur_tr_inds);
        nullMod_full = NMMfit_logexp_spkNL(nullMod_full,cur_Robs,Xmat,[],cur_full_inds);
        [nullMod_xvLL,null_xvLL] = NMMeval_model(nullMod,cur_Robs,Xmat,[],cur_xv_inds);
        nullMod.xvLLimp = (nullMod_xvLL - null_xvLL)/log(2);
        [nullMod_LL,null_LL] = NMMeval_model(nullMod_full,cur_Robs,Xmat,[],cur_full_inds);
        [nullMod_rptLL,null_rptLL] = NMMeval_model(nullMod_full,cur_Robs,Xmat,[],cur_rpt_inds);
        
        %% FIT (eye-corrected) STIM-PROCESSING MODEL
        fprintf('Fitting stim models for unit %d\n',cc);
        if unique(expt_data.expt_dds) == 12
            cur_base_lambda_d2XT = base_lambda_d2XT/5;
            cur_base_lambda_L1 = base_lambda_L1/2;
            cur_init_lambda_d2XT = init_lambda_d2XT/5;
        else
            cur_base_lambda_d2XT = base_lambda_d2XT;
            cur_base_lambda_L1 = base_lambda_L1;
            cur_init_lambda_d2XT = init_lambda_d2XT;
        end
        init_reg_params = NMMcreate_reg_params('lambda_d2XT',cur_init_lambda_d2XT,'boundary_conds',[0 0 0]);
        
        if cc <= params.n_probes %for MUA just fit a 3-filter quad model
            
            %start out with 3-filt quad (all exc)
            stim_mod_signs = [1 1 ones(1,MUA_nEfilts) -1*ones(1,MUA_nIfilts)]; %create quad model with one block-filter
            stim_NL_types = [{'lin','lin'} repmat({'quad'},1,MUA_nEfilts) repmat({'quad'},1,MUA_nIfilts)];
            stim_Xtargs = [2 ones(1,MUA_nEfilts+MUA_nIfilts+1)];
            init_mod = NMMinitialize_model( mod_stim_params, stim_mod_signs, stim_NL_types, init_reg_params,stim_Xtargs);
            init_mod.mods(1).reg_params = NMMcreate_reg_params(); %no regularization on block coefs
            init_mod = NMMfit_filters(init_mod,cur_Robs,Xmat,[],cur_tr_inds,silent);
            [~, ~, ~, ~, gint] = NMMeval_model(init_mod,cur_Robs,Xmat,[],cur_tr_inds);
            vgint = var(gint(:,stim_Xtargs==1));sgint = std(gint(:,stim_Xtargs==1));
            init_mod = NMMadjust_regularization(init_mod,find(stim_Xtargs == 1),'lambda_d2XT',cur_base_lambda_d2XT./vgint);
            init_mod = NMMadjust_regularization(init_mod,find(stim_Xtargs == 1),'lambda_L1',cur_base_lambda_L1./sgint);
            
            bestGQM = NMMfit_filters(init_mod,cur_Robs,Xmat,[],cur_tr_inds,silent);
            bestGQM_spkNL = NMMfit_logexp_spkNL(bestGQM,cur_Robs,Xmat,[],cur_tr_inds);
            if isempty(cur_xv_inds)
                error('Need xval trials');
            end
            bestGQM_xvLL = NMMeval_model(bestGQM_spkNL, cur_Robs, Xmat,[],cur_xv_inds);
            bestGQM.xvLLimp = (bestGQM_xvLL-nullMod_xvLL)/log(2);
            
        else %for SUA try a number of different model structures and pick the best by xvalLL
            
            %start out with just a linear filter
            nEfilts = 0;
            nIfilts = 0;
            
            stim_mod_signs = [1 1 ones(1,nEfilts) -1*ones(1,nIfilts)];
            stim_NL_types = [{'lin','lin'} repmat({'quad'},1,nEfilts) repmat({'quad'},1,nIfilts)];
            stim_Xtargs = [2 ones(1,nEfilts+nIfilts+1)];
            
            init_mod = NMMinitialize_model( mod_stim_params, stim_mod_signs, stim_NL_types, init_reg_params,stim_Xtargs);
            init_mod.mods(1).reg_params = NMMcreate_reg_params(); %no regularization on block coefs
            init_mod = NMMfit_filters(init_mod,cur_Robs,Xmat,[],cur_tr_inds,silent);
            [~, ~, ~, ~, gint] = NMMeval_model(init_mod,cur_Robs,Xmat,[],cur_tr_inds);
            vgint = var(gint(:,stim_Xtargs ==1));sgint = std(gint(:,stim_Xtargs==1));
            init_mod = NMMadjust_regularization(init_mod,find(stim_Xtargs==1),'lambda_d2XT',cur_base_lambda_d2XT./vgint);
            init_mod = NMMadjust_regularization(init_mod,find(stim_Xtargs==1),'lambda_L1',cur_base_lambda_L1./sgint);
            
            bestGQM = NMMfit_filters(init_mod,cur_Robs,Xmat,[],cur_tr_inds,silent);
            bestGQM_spkNL = NMMfit_logexp_spkNL(bestGQM,cur_Robs,Xmat,[],cur_tr_inds);
            if isempty(cur_xv_inds)
                error('Need xval trials');
            end
            bestGQM_xvLL = NMMeval_model(bestGQM_spkNL, cur_Robs, Xmat,[],cur_xv_inds);
            
            %now try adding up to max_Emods exc squared filters. Stop when xval
            %LL doesn't improve
            cur_imp = Inf;
            while nEfilts < max_Emods && cur_imp > 0
                cur_mod = bestGQM;
                cur_mod = NMMadd_NLinput(cur_mod,'quad',1,1,0.1*randn(flen*use_nPix_us/add_usfac,1)); %add a quadratic filter
                cur_mod.mods(end).reg_params = init_reg_params;
                nEfilts = nEfilts + 1;
                
                fprintf('Fitting model with LIN + %dE and %dI\n',nEfilts,nIfilts);
                cur_mod = NMMfit_filters(cur_mod,cur_Robs,Xmat,[],cur_tr_inds,silent);
                [~, ~, ~, ~, gint] = NMMeval_model(cur_mod,cur_Robs,Xmat,[],cur_tr_inds);
                vgint = var(gint); sgint = std(gint);
                cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_d2XT',cur_base_lambda_d2XT./vgint(end));
                cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_L1',cur_base_lambda_L1./sgint(end));
                cur_mod = NMMfit_filters(cur_mod,cur_Robs,Xmat,[],cur_tr_inds,silent);
                cur_mod_spkNL = NMMfit_logexp_spkNL(cur_mod,cur_Robs,Xmat,[],cur_tr_inds);
                cur_mod_xvLL = NMMeval_model(cur_mod_spkNL, cur_Robs, Xmat,[],cur_xv_inds);
                
                cur_imp = cur_mod_xvLL - bestGQM_xvLL;
                if cur_imp > 0
                    fprintf('Keeping new model\n');
                    bestGQM = cur_mod;
                    bestGQM_spkNL = cur_mod_spkNL;
                    bestGQM_xvLL = cur_mod_xvLL;
                else
                    nEfilts = nEfilts - 1;
                end
            end
            
            %now try adding up to max_Imods inh squared filters. Stop when xval
            %LL doesn't improve
            cur_imp = Inf;
            while nIfilts < max_Imods && cur_imp > 0
                cur_mod = bestGQM;
                cur_mod = NMMadd_NLinput(cur_mod,'quad',-1,1,0.1*randn(flen*use_nPix_us/add_usfac,1)); %add inhibitory quad filter
                cur_mod.mods(end).reg_params = init_reg_params;
                nIfilts = nIfilts + 1;
                
                fprintf('Fitting model with LIN + %dE and %dI\n',nEfilts,nIfilts);
                cur_mod = NMMfit_filters(cur_mod,cur_Robs,Xmat,[],cur_tr_inds,silent);
                [~, ~, ~, ~, gint] = NMMeval_model(cur_mod,cur_Robs,Xmat,[],cur_tr_inds);
                vgint = var(gint); sgint = std(gint);
                cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_d2XT',cur_base_lambda_d2XT./vgint(end));
                cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_L1',cur_base_lambda_L1./sgint(end));
                cur_mod = NMMfit_filters(cur_mod,cur_Robs,Xmat,[],cur_tr_inds,silent);
                cur_mod_spkNL = NMMfit_logexp_spkNL(cur_mod,cur_Robs,Xmat,[],cur_tr_inds);
                cur_mod_xvLL = NMMeval_model(cur_mod_spkNL, cur_Robs, Xmat,[],cur_xv_inds);
                
                cur_imp = cur_mod_xvLL - bestGQM_xvLL;
                if cur_imp > 0
                    fprintf('Keeping new model\n');
                    bestGQM = cur_mod;
                    bestGQM_spkNL = cur_mod_spkNL;
                    bestGQM_xvLL = cur_mod_xvLL;
                end
            end
            
            bestGQM.xvLLimp = (bestGQM_xvLL-nullMod_xvLL)/log(2);
        end
        %% Fit a GQM with the linear filter split into two rectified subunits
        if fit_rect
            %convert linear filter into separate thresh-lin filters and refit
            rectGQM = bestGQM;
            stim_filt_set = find([rectGQM.mods(:).Xtarget] == 1);
            rectGQM.mods = [rectGQM.mods(:); rectGQM.mods(stim_filt_set(1))]; %add a copy of the lin filter
            rectGQM.mods(end).sign = -1; %set it to be inhibitory
            rectGQM.mods(end).filtK = -rectGQM.mods(end).filtK; %initialize filter to be negative of the linear filter
            stim_filt_set = find([rectGQM.mods(:).Xtarget] == 1);
            rectGQM.mods(stim_filt_set(1)).NLtype = 'threshlin'; rectGQM.mods(stim_filt_set(end)).NLtype = 'threshlin'; 
            rectGQM = NMMfit_filters(rectGQM,cur_Robs,Xmat,[],cur_tr_inds,silent);
            rectGQM_spkNL = NMMfit_logexp_spkNL(rectGQM,cur_Robs,Xmat,[],cur_tr_inds);
            rectGQM_xvLL = NMMeval_model(rectGQM_spkNL, cur_Robs, Xmat,[],cur_xv_inds);
            rectGQM.xvLLimp = (rectGQM_xvLL-nullMod_xvLL)/log(2);
        end
        %% try different smoothness reg strengths
        if ~isempty(poss_smoothreg_scalefacs)
            base_mod = bestGQM;
            reg_optim_params.maxIter = 250;
            stim_filt_set = find([base_mod.mods(:).Xtarget] == 1);
            base_d2XT_vals = arrayfun(@(x) x.reg_params.lambda_d2XT,base_mod.mods(stim_filt_set));
            reg_val_xvLL = nan(length(poss_smoothreg_scalefacs),1); %initialize vector of xvLL at each reg strength
            clear all_reg_mods
            cur_mod = base_mod;
            for pp = 1:length(poss_smoothreg_scalefacs)
                fprintf('Fitting reg scale %d of %d\n',pp,length(poss_smoothreg_scalefacs));
                cur_reg_scale = poss_smoothreg_scalefacs(pp);
                cur_mod = NMMadjust_regularization(cur_mod,stim_filt_set,'lambda_d2XT',base_d2XT_vals*cur_reg_scale);
                cur_mod = NMMfit_filters(cur_mod,cur_Robs,Xmat,[],cur_tr_inds,silent,reg_optim_params);
                cur_mod_spkNL = NMMfit_logexp_spkNL(cur_mod,cur_Robs,Xmat,[],cur_tr_inds);
                reg_val_xvLL(pp) = NMMeval_model(cur_mod_spkNL,cur_Robs,Xmat,[],cur_xv_inds);
                all_reg_mods(pp) = cur_mod;
            end
            [best_xvLL,best_reg_ind] = max(reg_val_xvLL); %find reg strength with best xvLL performance
            best_reg_scale = poss_smoothreg_scalefacs(best_reg_ind);
            if best_xvLL > bestGQM_xvLL %if we did better than our original model (with reg scale 1)
                bestGQM = all_reg_mods(best_reg_ind);
                bestGQM.xvLLimp = (best_xvLL - nullMod_xvLL)/log(2);
            else
                best_reg_scale = 1; %if the original model (at scale == 1) was a better fit, keep it
            end
        else
            all_reg_mods = [];
            reg_val_xvLL = [];
            best_reg_scale = nan;
        end
        
        %% Refit model parameters on the full data set (excluding repeats)
        if ~isempty(cur_xv_inds)
            bestGQM = NMMfit_filters(bestGQM,cur_Robs,Xmat,[],cur_full_inds,silent);
            bestGQM_spkNL = NMMfit_logexp_spkNL(bestGQM,cur_Robs,Xmat,[],cur_full_inds);
            [bestGQM_LL, nullLL, ~, ~, ~, fgint] = NMMeval_model(bestGQM_spkNL, cur_Robs, Xmat,[],cur_full_inds);
            stim_filt_set = find([bestGQM.mods(:).Xtarget] == 1);
            rel_filt_weights = std(fgint(:,stim_filt_set));
            bestGQM_spkNL.rel_filt_weights = rel_filt_weights/sum(rel_filt_weights);
            bestGQM_spkNL.LLimp = (bestGQM_LL - nullMod_LL)/log(2);
            bestGQM_rptLL = NMMeval_model(bestGQM_spkNL, cur_Robs, Xmat,[],cur_rpt_inds);
            bestGQM_spkNL.rptLLimp = (bestGQM_rptLL - nullMod_rptLL)/log(2); %store LL imp on rpt trials
            
            if fit_rect
                rectGQM = NMMfit_filters(rectGQM,cur_Robs,Xmat,[],cur_full_inds,silent);
                [~,~,~,~,~,fgint] = NMMeval_model(rectGQM,cur_Robs,Xmat,[],cur_full_inds);
                stim_filt_set = find([rectGQM.mods(:).Xtarget] == 1);
                rel_filt_weights = std(fgint(:,stim_filt_set));
                rectGQM.rel_filt_weights = rel_filt_weights/sum(rel_filt_weights);
                rectGQM_spkNL = NMMfit_logexp_spkNL(rectGQM,cur_Robs,Xmat,[],cur_full_inds);
                rectGQM.LLimp = (rectGQM_spkNL.LL_seq(end)-nullMod_LL)/log(2);
                rectGQM_rptLL = NMMeval_model(rectGQM_spkNL, cur_Robs, Xmat,[],cur_rpt_inds);
                rectGQM.rptLLimp = (rectGQM_rptLL - nullMod_rptLL)/log(2); %store LL imp on rpt trials
            end
        end
        
        %% Store data
        fullLL = nansum(cur_Robs(cur_full_inds).*log(cur_Robs(cur_full_inds)) - cur_Robs(cur_full_inds));
        fullxvLL = nansum(cur_Robs(cur_xv_inds).*log(cur_Robs(cur_xv_inds)) - cur_Robs(cur_xv_inds));
        unit_data.fullLL = fullLL;
        unit_data.null_LL = null_LL*nansum(cur_Robs(cur_full_inds));
        unit_data.nullMod_LL = nullMod_LL*nansum(cur_Robs(cur_full_inds));
        unit_data.fullxvLL = fullxvLL;
        unit_data.nullxvLL = null_xvLL*nansum(cur_Robs(cur_xv_inds));
        unit_data.nullModxvLL = nullMod_xvLL*nansum(cur_Robs(cur_xv_inds));
        unit_data.tr_trials = tr_trials;
        unit_data.xv_trials = xv_trials;
        
        ModData(cc).unit_data = unit_data;
        ModData(cc).bestGQM = bestGQM_spkNL;
        ModData(cc).bestGQM.xvLLimp = bestGQM.xvLLimp;
        ModData(cc).nullMod = nullMod_full;
        if ~isempty(poss_smoothreg_scalefacs)
            ModData(cc).all_reg_mods = all_reg_mods;
            ModData(cc).reg_val_xvLL = reg_val_xvLL;
            ModData(cc).best_reg_scale = best_reg_scale;
            ModData(cc).poss_smoothreg_scalefacs = poss_smoothreg_scalefacs;
        end
        if fit_rect
            ModData(cc).rectGQM = rectGQM_spkNL;
        end
        %% CALCULATE TUNING PROPERTIES
        clear tune_props
        if fit_rect
            use_mod = ModData(cc).rectGQM;
        else
            use_mod = ModData(cc).bestGQM;
        end
        
        %remove the block-by-block variability in the model by
        %incorporating the avg output of the block-filter
        cur_block_filt = use_mod.mods(1).filtK; %filter applied to the block index
        cur_used_blocks = ModData(cc).unit_data.used_blocks; %which blocks was this neuron isolated during
        poss_used_blocks = ModData(cc).unit_data.poss_used_blocks; %total set of used blocks
        cur_used_blocks = find(ismember(poss_used_blocks,cur_used_blocks)); %indices of blocks where this neuron was isolated
        use_mod_no_block = use_mod;
        use_mod_no_block.spk_NL_params(1) = use_mod_no_block.spk_NL_params(1) + mean(cur_block_filt(cur_used_blocks)); %add the average output of the block-filter to the spkNL offset
        use_mod_no_block.mods(1) = []; %eliminate block filter
        
        [~,~,best_pred_rate,~,~,filt_outs] = NMMeval_model(use_mod_no_block,cur_Robs,Xmat,[],cur_full_inds);
        tempX = Xmat; tempX{1} = -tempX{1};
        [~,~,rev_pred_rate] = NMMeval_model(use_mod_no_block,cur_Robs,tempX,[],cur_full_inds);
        clear tempX
        tune_props.PRM = mean(abs(best_pred_rate - rev_pred_rate))/mean(best_pred_rate); %mean absolute deviation caused by phase-reversal divided by mean rate
        tune_props.PRI = std(best_pred_rate - rev_pred_rate)/std(best_pred_rate); %SD of phase-reversal modulation divided by total rate SD
        
        %get filters and filter signs
         stim_filters = [use_mod_no_block.mods(:).filtK];
        stim_mod_signs = [use_mod_no_block.mods(:).sign];

        filt_out_weights = std(filt_outs);
        filt_data = get_filter_properties_v2(stim_filters,mod_stim_params(1),mod_dx,filt_out_weights,stim_mod_signs);
        tune_props.filt_data = filt_data;
        
        %only use E/lin filters for these calculations
        non_supp_filts = find(stim_mod_signs ~= -1);
        filt_out_weights = filt_out_weights(non_supp_filts)/nansum(filt_out_weights(non_supp_filts));

        tune_props.avgRF_mean = filt_data.avg_gauss_mean - use_nPix_us*sp_dx/2 -sp_dx;
        tune_props.avgRF_sigma = filt_data.avg_gauss_std;
        tune_props.RF_mean = nansum(filt_out_weights(non_supp_filts).*filt_data.gest(non_supp_filts,1)') - use_nPix_us*sp_dx/2 -sp_dx;
        tune_props.RF_sigma = nansum(filt_out_weights(non_supp_filts).*filt_data.gest(non_supp_filts,2)');
        tune_props.RF_gSF = nansum(filt_out_weights(non_supp_filts).*filt_data.gest(non_supp_filts,3)');
        tune_props.RF_dirsel = nansum(filt_out_weights(non_supp_filts).*filt_data.dir_selectivity(non_supp_filts)');
        tune_props.RF_FTF = nansum(filt_out_weights(non_supp_filts).*filt_data.FFt(non_supp_filts)');
        tune_props.RF_FSF = nansum(filt_out_weights(non_supp_filts).*filt_data.FFx(non_supp_filts)');
        
        if fit_rect
            rect_stim_filters = [use_mod_no_block.mods([1 end]).filtK];
            avg_rect_filts = mean(rect_stim_filters);
            absavg_rect_filts = mean(abs(rect_stim_filters(:)));
            tune_props.net_phase_polarity = (avg_rect_filts(1) - avg_rect_filts(end))/absavg_rect_filts;
        end
        
        screen_X =  -tune_props.RF_mean'.*sind(bar_ori) + ov_RF_pos(1); %based on weighted avg of individual filter gabor fits
        screen_Y =  tune_props.RF_mean'.*cosd(bar_ori) + ov_RF_pos(2);
        screen_X_avg =  -tune_props.avgRF_mean'.*sind(bar_ori) + ov_RF_pos(1); %based on gaussian fit to weighted avg filter
        screen_Y_avg =  tune_props.avgRF_mean'.*cosd(bar_ori) + ov_RF_pos(2);
%         if strcmp(rec_type,'UA')
%             if bar_ori == 0
%                 screen_X = interp_x(unit_data.probe_number);
%             elseif bar_ori == 90
%                 screen_Y = interp_y(unit_data.probe_number);
%             else
%                 error('Code doesnt work with this condition!');
%             end
%         end
        
        tune_props.screen_X = screen_X;
        tune_props.screen_Y = screen_Y;
        tune_props.screen_X_avg = screen_X_avg;
        tune_props.screen_Y_avg = screen_Y_avg;
        tune_props.RF_ecc = sqrt(screen_X.^2 + screen_Y.^2);
        tune_props.RF_ecc_avg = sqrt(screen_X_avg.^2 + screen_Y_avg.^2);
        tune_props.filt_out_weights = filt_out_weights;
        ModData(cc).tune_props = tune_props;
    end
end

%%
interp_EP = ET_data.interp_eye_pos(used_inds,:);
par_pos = mean(interp_EP(:,[1 3]),2);


%%
cc = 35;
cur_mod = ModData(cc).bestGQM;
loo_cc = find(loo_set == cc); %index within the LOOXV set
cur_Robs = Robs_mat(:,cc);
cc_uinds = find(~isnan(cur_Robs)); %usable time points for this unit

[LL,nullLL,pred_rate,G,gint,fgint] = NMMeval_model(cur_mod, cur_Robs, Xmat,[],cc_uinds);


Xtick = linspace(-6,6,50);
n_Gbins = 35;
n_Tbins = length(Xtick); %number of time bins (with buffer)
G_lambda = 10; %d2T reg on G dimension
P_lambda = 10;
silent = 0;

%initialize 2D TB data using t-since-sac and normalized G
TB_stim = [par_pos G];

%set G-bins based on prctiles
%         %equispaced binning
Ytick = my_prctile(G(cc_uinds),linspace(0.5,99.5,n_Gbins)); %equipopulated binning

%in some cases there are an excess of G==0 values, which makes
%some adjacent bins identical. Need to cut these bin edges.
nd_bins = find(diff(Ytick) <= 0);
Ytick(nd_bins) = [];
cur_nGbins = length(Ytick);

%initialize TBs
TB = TentBasis2D(Xtick, Ytick);

%this is data within range of the TB centers
TB_used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
    TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
% TB_used_data = TB_used_data(ismember(TB_used_data,fit_inds));

null_rate = mean(cur_Robs(TB_used_data));
nullLL = nansum(cur_Robs(TB_used_data).*log2(null_rate) - null_rate)/sum(cur_Robs(TB_used_data));

%process data with TBs
[TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(TB_used_data,:));
%this gives fixed reg strength in G-direction, but
%cur_sac_lambda in T-direction
L2_params = create_L2_params([],[1 n_Tbins*cur_nGbins],[n_Tbins cur_nGbins],2,3,[Inf Inf],[1 1/P_lambda*G_lambda]);
TB_sacmod = regGLM_fit(TB_Xmat,cur_Robs(TB_used_data),L2_params,P_lambda,[],[],silent);
[LL, penLL, TB_pred_rate] = regGLM_eval(TB_sacmod,cur_Robs(TB_used_data),TB_Xmat);

%compute output of TB model
TB_K = reshape(TB_sacmod.K,n_Tbins,cur_nGbins)';
bin_areas = TB.GetBinAreas();
TB_dist = TB_counts./bin_areas;
normFac = trapz(trapz(TB_dist));
TB_dist = TB_dist'/normFac;
TB_rate = log(1 + exp(TB_K + TB_sacmod.theta));

%%
% modFitParams = struct('bar_ori',bar_ori,'use_MUA',use_MUA,'fit_uncor',fit_unCor,...
%     'xv_frac',xv_frac,'xv_type',xv_type,'use_nPix_us',use_nPix_us,'flen',flen,...
%     'spatial_usfac',spatial_usfac,'sp_dx',sp_dx,'mod_dx',mod_dx,'add_usfac',add_usfac,'dt',dt,...
%     'base_lambda_d2XT',base_lambda_d2XT,'base_lambda_L1',base_lambda_L1,'init_lambda_d2XT',init_lambda_d2XT);
% 
% cd(save_dir)
% save(save_name,'ModData','modFitParams');