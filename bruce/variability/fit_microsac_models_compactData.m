% clear all
% close all

global Expt_name bar_ori monk_name rec_type

% Expt_name = 'M011';
% monk_name = 'jbe';
% bar_ori = 160; %bar orientation to use (only for UA recs)
% rec_number = 1;

fit_unCor = false; %use eye correction
use_MUA = false; %use MUA in model-fitting
include_bursts = true;
fit_preFilter = false;

save_name = 'msac_models_compData';


mod_name = 'corrected_models_comp';
if fit_unCor
    mod_name = strcat(mod_name,'_unCor');
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
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
fprintf('Loading %s\n',data_name);
load(data_name);

%%
sac_burst_isi = 0.15; %minimum inter-saccade interval for micros (to eliminate 'bursts')
micro_thresh = 1; %max amp of microsac (deg)
gsac_thresh = 1;
max_gsac_dur = 0.1; %maximum saccade duration before we call it a likely blink

backlag = round(0.1/params.dt);
forlag = round(0.3/params.dt);
slags = -backlag:forlag;
n_sac_bins = length(slags);

off_d2T = 100;
gain_d2T = 50;

sacParams.fit_unCor = fit_unCor;
sacParams.fit_preFilter = fit_preFilter;
sacParams.include_bursts = include_bursts;
sacParams.off_d2T = off_d2T;
sacParams.gain_d2T = gain_d2T;
sacParams.slags = slags;
sacParams.sac_burst_isi = sac_burst_isi;
sacParams.micro_thresh = micro_thresh;
sacParams.gsac_thresh = gsac_thresh;
sacParams.max_gsac_dur = max_gsac_dur;
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
mod_name = [mod_name sprintf('_ori%d',bar_ori)];
save_name = [save_name sprintf('_ori%d',bar_ori)];

if rec_number > 1
    et_mod_data_name = strcat(et_mod_data_name,sprintf('r%d',rec_number));
    et_anal_name = strcat(et_anal_name,sprintf('r%d',rec_number));
    mod_name = strcat(mod_name,sprintf('_r%d',rec_number));
    save_name = strcat(save_name,sprintf('_r%d',rec_number));
end

%% load model fits
cd(save_dir)
load(mod_name);

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
spatial_usfac = et_params.spatial_usfac; %spatial up-sampling factor

use_nPix = et_params.use_nPix; %number of pixels (bars) used in models
full_nPix = params.full_nPix; %total number of bars to keep track of in stimulus
dt = params.dt;

use_nPix_us = use_nPix*spatial_usfac; %number of pixels in stimulus filters
klen_us = use_nPix_us*flen; %number of parameters in stim-filters
sp_dx = et_params.sp_dx; %pixel size (deg)

use_LOOXV = 1; %[0 is no LOO; 1 is SUs only; 2 is SU + MU]
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
saccades = ET_data.saccades(used_saccade_set);
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);
used_is_blink = ET_data.is_blink(used_saccade_set);

sac_amps = [saccades(:).amplitude];
sac_direction = [saccades(:).direction];
sac_durs = [saccades(:).duration];
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

sacburst_set = find([saccades(:).isi] < sac_burst_isi | [saccades(:).next_isi] < sac_burst_isi);
micro_sacs = find([saccades(:).amplitude] < micro_thresh & ~used_is_blink');

msac_bursts = micro_sacs(ismember(micro_sacs,sacburst_set));
if ~include_bursts
    micro_sacs(ismember(micro_sacs,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'
end

%guided saccades are those whose parallel component is large enough and
%that aren't blinks (and whose duration is not too long to be suspicious
big_sacs = find(abs(sac_deltaX) > gsac_thresh & ~used_is_blink' & sac_durs <= max_gsac_dur);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

%% MAKE TIME-EMBEDDED SAC PREDICTOR MATS
Xsac = zeros(NT,length(slags));
Xmsac = zeros(NT,length(slags));
for ii = 1:n_sac_bins
    cur_sac_target = saccade_start_inds(big_sacs) + slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(big_sacs(uu))) = [];
    Xsac(cur_sac_target,ii) = 1;
    
    cur_sac_target = saccade_start_inds(micro_sacs) + slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(micro_sacs(uu))) = [];
    Xmsac(cur_sac_target,ii) = 1;
end

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
if ~fit_unCor
    best_shift_stimmat_up = all_stimmat_up;
    for i=1:NT
        best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
    end
    all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_full);
else
    all_Xmat_shift = create_time_embedding(all_stimmat_up,stim_params_full);
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
mod_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
mod_stim_params(2) = NMMcreate_stim_params([n_blocks],dt);
mod_stim_params(3) = NMMcreate_stim_params([n_sac_bins],dt);
Xmat{1} = all_Xmat_shift(used_inds,use_kInds_up);
Xmat{2} = Xblock(used_inds,:);
Xmat{3} = Xmsac;

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
        cur_full_inds = cc_uinds(ismember(all_trialvec(used_inds(cc_uinds)),use_trials));
        
        cur_Xsac = Xmsac; %saccade indicator Xmat
        
        %only use indices during guided saccade expts here
        any_sac_inds = cc_uinds(any(cur_Xsac(cc_uinds,:) > 0,2));
        
        
        %% RECON RETINAL STIM
        if ~fit_unCor & ~isempty(loo_cc)
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
            Xmat{1} = all_Xmat_shift(used_inds,use_kInds_up);
        end
        
        %% Fit spk NL params and refit scale of each filter using target data (within trange of sacs)
        cur_GQM = ModData(cc).bestGQM;
        cur_GQM.mods(1) = []; %eliminate block filter
        cur_GQM = NMMfit_logexp_spkNL(cur_GQM,cur_Robs,Xmat,[],any_sac_inds); %refit spk NL
        
        stim_filters = find([cur_GQM.mods(:).Xtarget] == 1);
        stim_mod_signs = [cur_GQM.mods(:).sign];
        [~,~,basemod_pred_rate,~,filt_outs,fgint] = NMMeval_model(cur_GQM,cur_Robs,Xmat);
        fgint = bsxfun(@times,fgint,stim_mod_signs);
        stimG = sum(fgint(:,stim_filters),2);
        
        sacMod(cc).base_mod = cur_GQM;
        sacMod(cc).msac_ovavg_rate = mean(cur_Robs(any_sac_inds));
        [baseLL,~,~,~,~,~,base_nullLL] = NMMeval_model(cur_GQM,cur_Robs,Xmat,[],any_sac_inds);
        sacMod(cc).msac_base_LLimp = (baseLL-base_nullLL)/log(2);
        
        %% fit msac gain/offset models (post-filtering)
        %fit offset only model
        [sacMod(cc).msac_off_mod] = fit_sacMod_basic(cur_GQM,cur_Robs,cur_Xsac,stimG,any_sac_inds,off_d2T,gain_d2T,true);
        
        %single post-gain filter with offset
        [sacMod(cc).msac_post_mod] = fit_sacMod_basic(cur_GQM,cur_Robs,cur_Xsac,stimG,any_sac_inds,off_d2T,gain_d2T);
        
        %% fit pre-filtering model
        if fit_preFilter
            Xsac_mat = cur_Xsac(any_sac_inds,:);
            [sacMod(cc).msacPreGainMod] = fit_pre_gainmodel_basic...
                (cur_GQM,cur_Robs(any_sac_inds),Xmat{1}(any_sac_inds,:),Xsac_mat,off_d2T,gain_d2T);
        end
        
        %% compute basic saccade-triggered stats
        sacMod(cc).msac_avgrate = nan(length(slags),1);
        sacMod(cc).msac_Nspks = nan(length(slags),1);
        for ii = 1:length(slags)
            temp = find(cur_Xsac(:,ii) == 1);
            sacMod(cc).msac_avgrate(ii) = mean(cur_Robs(temp));
            sacMod(cc).msac_Nspks(ii) = sum(cur_Robs(temp));
        end
        
    end
end

%%
cd(save_dir)
save(save_name,'sacMod','sacParams');


