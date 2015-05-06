% clear all
% close all

global Expt_name bar_ori monk_name rec_type 

% Expt_name = 'M012';
% monk_name = 'jbe';
% bar_ori = 0; %bar orientation to use (only for UA recs)
% rec_number = 1;

use_MUA = false; %use MUA in model-fitting
fit_rect = false; %split quad linear filter into two rectified

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
ff = find(cellfun(@(x) ~isempty(x),Expts),1);
if strcmp(Expts{ff}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{ff}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

%load in packaged data
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
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
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';

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

%% LOAD IN PREVIOUS MODEL FITS
cd(save_dir)
load(save_name);

%%
xv_type = modFitParams.xv_type;
xv_frac = modFitParams.xv_frac;

flen = modFitParams.flen; %time lags for ST filters
spatial_usfac = modFitParams.spatial_usfac; %spatial up-sampling factor

full_nPix = params.full_nPix; %total number of bars to keep track of in stimulus
dt = params.dt;

use_nPix_us = modFitParams.use_nPix_us; %number of pixels in stimulus filters
klen_us = use_nPix_us*flen; %number of parameters in stim-filters
sp_dx = modFitParams.sp_dx; %pixel size (deg)
use_nPix = use_nPix_us/spatial_usfac;

all_stim_mat = decompressTernNoise(stimComp);

ov_RF_pos = Expts{expt_data.used_blocks(1)}.Stimvals.rf(1:2)/params.scale_fac;

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

all_Xmat_unCor = create_time_embedding(all_stimmat_up,stim_params_full);

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
% n_blocks = length(expt_data.used_blocks);
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
end

%%

%for model fitting
if use_MUA
    targs = 1:n_units; %SU and MU
else
    targs = setdiff(1:n_units,1:params.n_probes); %SU only
end
%%
mod_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
mod_stim_params(2) = NMMcreate_stim_params([n_blocks],dt);
Xmat{1} = all_Xmat_unCor(used_inds,use_kInds_up);
Xmat{2} = Xblock(used_inds,:);

%%
for cc = targs
    fprintf('Starting model fits for unit %d\n',cc);
    
    cur_Robs = Robs_mat(:,cc);
    cc_uinds = find(~isnan(cur_Robs)); %usable time points for this unit
    
    if ~isempty(cc_uinds)
        
        use_trials = unique(all_trialvec(used_inds(cc_uinds))); %set of potentially usable trials
        use_trials(ismember(use_trials,rpt_trials)) = []; %dont use repeat trials
        
        tr_trials = ModData(cc).unit_data.tr_trials;
        xv_trials = ModData(cc).unit_data.xv_trials;
                
        cur_tr_inds = cc_uinds(ismember(all_trialvec(used_inds(cc_uinds)),tr_trials));
        cur_xv_inds = cc_uinds(ismember(all_trialvec(used_inds(cc_uinds)),xv_trials));
        cur_full_inds = cc_uinds(ismember(all_trialvec(used_inds(cc_uinds)),use_trials));
        
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
        
        %%
        nullMod = NMMinitialize_model(mod_stim_params,1,{'lin'},[],2);
        nullMod = NMMfit_filters(nullMod,cur_Robs,Xmat,[],cur_tr_inds,silent);
        nullMod = NMMfit_logexp_spkNL(nullMod,cur_Robs,Xmat,[],cur_tr_inds);
        [nullMod_xvLL,null_xvLL] = NMMeval_model(nullMod,cur_Robs,Xmat,[],cur_xv_inds);
        nullMod.xvLLimp = (nullMod_xvLL - null_xvLL)/log(2);
                
        %% Fit uncorrected model
        bestGQM_unCor = ModData(cc).bestGQM;
        bestGQM_unCor.spk_NL_params = [0 1 1 0]; %reset spk NL params to initial
        bestGQM_unCor = NMMfit_filters(bestGQM_unCor,cur_Robs,Xmat,[],cur_full_inds);
        bestGQM_unCor = NMMfit_logexp_spkNL(bestGQM_unCor,cur_Robs,Xmat,[],cur_full_inds);
        [cur_LL,nullLL] = NMMeval_model(bestGQM_unCor,cur_Robs,Xmat,[],cur_full_inds);
        bestGQM_unCor.LLimp = (cur_LL - nullLL)/log(2);
                
        %% Store data
        ModData(cc).bestGQM_unCor = bestGQM_unCor;
        ModData(cc).nullMod = nullMod;
        %% CALCULATE UNCORRECTED TUNING PROPERTIES
        clear tune_props_unCor
        use_mod = ModData(cc).bestGQM_unCor;

        stim_filt_set = find([use_mod.mods(:).Xtarget] == 1);
        stim_filters = [use_mod.mods(stim_filt_set).filtK];
        stim_mod_signs = [use_mod.mods(stim_filt_set).sign];
        filt_data = get_filter_properties_v2(stim_filters,mod_stim_params(1),sp_dx);
        tune_props_unCor.filt_data = filt_data;
        
        [~,~,best_pred_rate,~,~,filt_outs] = NMMeval_model(use_mod,cur_Robs,Xmat,[],cur_full_inds);
        tempX = Xmat; tempX{1} = -tempX{1};
        [~,~,rev_pred_rate] = NMMeval_model(use_mod,cur_Robs,tempX,[],cur_full_inds);
        tune_props_unCor.PRM = mean(abs(best_pred_rate - rev_pred_rate))/mean(best_pred_rate);
        
        %only use E/lin filters for these calculations
        non_supp_filts = find(stim_mod_signs ~= -1);
        filt_out_weights = std(filt_outs(:,stim_filt_set));
        filt_out_weights = filt_out_weights(non_supp_filts)/nansum(filt_out_weights(non_supp_filts));
        
        tune_props_unCor.RF_mean = filt_out_weights(non_supp_filts)*filt_data.gest(non_supp_filts,1) - use_nPix_us*sp_dx/2 -sp_dx;
        tune_props_unCor.RF_sigma = filt_out_weights(non_supp_filts)*filt_data.gest(non_supp_filts,2);
        tune_props_unCor.RF_gSF = filt_out_weights(non_supp_filts)*filt_data.gest(non_supp_filts,3);
        tune_props_unCor.RF_dirsel = filt_out_weights(non_supp_filts)*filt_data.dir_selectivity(non_supp_filts);
        tune_props_unCor.RF_FTF = filt_out_weights(non_supp_filts)*filt_data.FFt(non_supp_filts);
        tune_props_unCor.RF_FSF = filt_out_weights(non_supp_filts)*filt_data.FFx(non_supp_filts);
        
        rect_stim_filters = [use_mod.mods(stim_filt_set([1 end])).filtK];
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
        tune_props_unCor.filt_out_weights = filt_out_weights;
        
        ModData(cc).tune_props_unCor = tune_props_unCor;
    end
end

%%
fit_unCor = true;
modFitParams = struct('bar_ori',bar_ori,'use_MUA',use_MUA,'fit_uncor',fit_unCor,...
    'xv_frac',xv_frac,'xv_type',xv_type,'use_nPix_us',use_nPix_us,'flen',flen,'spatial_usfac',spatial_usfac,'sp_dx',sp_dx,'dt',dt,...
    'base_lambda_d2XT',base_lambda_d2XT,'base_lambda_L1',base_lambda_L1,...
    'init_lambda_d2XT',init_lambda_d2XT);

cd(save_dir)
save(save_name,'ModData','modFitParams');