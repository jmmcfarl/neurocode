clear all
close all

global Expt_name bar_ori monk_name rec_type

Expt_name = 'M012';
monk_name = 'jbe';
bar_ori = 0; %bar orientation to use (only for UA recs)
rec_number = 1;

use_MUA = false;
fit_unCor = false;
use_sacMods = false;

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

sname = 'model_variability_compact';
if fit_unCor
    sname = strcat(sname,'_unCor');
end

%%
Expt_num = str2num(Expt_name(2:end));

%load in array RF position data
load ~/Data/bruce/general_array_data/array_pos_data.mat
interp_ecc = sqrt(interp_x.^2 + interp_y.^2);

cd(data_dir);
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
model_dir = ['~/Analysis/bruce/' Expt_name '/models'];
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';

if rec_number > 1
    cluster_dir = [cluster_dir sprintf('/rec%d',rec_number)];
end

mod_name = 'corrected_models_comp';
sacMod_name = 'msac_models_compData';
if fit_unCor
    mod_name = strcat(mod_name,'_unCor');
    sacMod_name = strcat(sacMod_name,'_unCor');
end


%if using coil info
if any(params.use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
mod_name = [mod_name sprintf('_ori%d',bar_ori)];
sacMod_name = [sacMod_name sprintf('_ori%d',bar_ori)];

if rec_number > 1
    mod_name = strcat(mod_name,sprintf('_r%d',rec_number));
    sacMod_name = strcat(sacMod_name,sprintf('_r%d',rec_number));
    et_mod_data_name = strcat(et_mod_data_name,sprintf('r%d',rec_number));
    et_anal_name = strcat(et_anal_name,sprintf('r%d',rec_number));
end

load([model_dir '/' mod_name]);
if use_sacMods
    load([model_dir '/' sacMod_name]);
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
saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));
used_is_blink = ET_data.is_blink(used_saccade_set);

if use_sacMods
sac_amps = [saccades(:).amplitude];
sac_direction = [saccades(:).direction];
sac_durs = [saccades(:).duration];
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

sacburst_set = find([saccades(:).isi] < sacParams.sac_burst_isi | [saccades(:).next_isi] < sacParams.sac_burst_isi);
micro_sacs = find([saccades(:).amplitude] < sacParams.micro_thresh & ~used_is_blink');

msac_bursts = micro_sacs(ismember(micro_sacs,sacburst_set));
if ~sacParams.include_bursts
    micro_sacs(ismember(micro_sacs,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'
end

%guided saccades are those whose parallel component is large enough and
%that aren't blinks (and whose duration is not too long to be suspicious
big_sacs = find(abs(sac_deltaX) > sacParams.gsac_thresh & ~used_is_blink' & sac_durs <= sacParams.max_gsac_dur);


%% MAKE TIME-EMBEDDED SAC PREDICTOR MATS
slags = sacParams.slags;
n_sac_bins = length(slags);
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

%% Recon retinal stim
sac_shift = et_params.sac_shift; %forward projection of saccade start times

cur_fix_post_mean = squeeze(it_fix_post_mean(end,:));
cur_fix_post_std = squeeze(it_fix_post_std(end,:));
cur_drift_post_mean = squeeze(drift_post_mean(end,:));
cur_drift_post_std = squeeze(drift_post_std(end,:));
[fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
    cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);

fin_shift_cor = round(fin_tot_corr);

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
targs(targs > length(ModData)) = [];

%% MARK INDICES DURING BLINKS AND SACCADES
sac_buff = round(0.06/dt);
sac_delay = round(0.04/dt);
blink_buff = round(0.1/dt);

in_sac_inds = zeros(NT,1);
nblink_start_inds = saccade_start_inds(~used_is_blink);
for ii = 1:(sac_buff+1)
    cur_inds = nblink_start_inds + sac_delay + (ii-1);
    uu = find(cur_inds <= NT);
    uu(all_trialvec(used_inds(cur_inds(uu))) ~= all_trialvec(used_inds(nblink_start_inds(uu)))) = [];
    in_sac_inds(cur_inds(uu)) = 1;
end
in_sac_inds = logical(in_sac_inds);

blink_start_inds = saccade_start_inds(used_is_blink);
in_blink_inds = zeros(NT,1);
for ii = 1:(blink_buff+1)
    cur_inds = blink_start_inds + sac_delay + (ii-1);
    uu = find(cur_inds <= NT);
    uu(all_trialvec(used_inds(cur_inds(uu))) ~= all_trialvec(used_inds(blink_start_inds(uu)))) = [];
    in_blink_inds(cur_inds(uu)) = 1;
end
in_blink_inds = logical(in_blink_inds);

%% get set of trials used for calcs, and construct TBT matrices
use_trials = unique(all_trialvec(used_inds)); %set of potentially usable trials
use_trials(ismember(use_trials,rpt_trials)) = []; %dont use repeat trials

%use only completed trials
target_uf = 400 - (params.beg_buffer + params.end_buffer)/params.dt;
T = tabulate(all_trialvec(used_inds));
complete_trials = find(T(:,2) >= target_uf);
use_trials(~ismember(use_trials,complete_trials)) = [];

full_uinds = find(ismember(all_trialvec(used_inds),use_trials));
n_utrials = length(full_uinds)/target_uf;

%matrix of trial-by-trial EP estimates
base_EP_est = fin_tot_corr(full_uinds) - nanmedian(fin_tot_corr(full_uinds));
EP_tbt = reshape(base_EP_est,target_uf,n_utrials);

%TBT mats for sac and blink indicators
inblink_tbt = reshape(in_blink_inds(full_uinds),target_uf,n_utrials);
insac_tbt = reshape(in_sac_inds(full_uinds),target_uf,n_utrials);

if use_sacMods
Xmsac = Xmsac(full_uinds,:);
Xmsac_tbt = reshape(Xmsac,target_uf,n_utrials,length(slags));
end
%% absorb block filter into spkNL offset parameter
    
%COMPUTE XMAT
cur_shift_stimmat_up = all_stimmat_up;
for i=1:length(full_uinds)
    cur_shift_stimmat_up(used_inds(full_uinds(i)),:) = shift_matrix_Nd(all_stimmat_up(used_inds(full_uinds(i)),:),-fin_shift_cor(i),2);
end
all_Xmat_shift = create_time_embedding(cur_shift_stimmat_up,stim_params_full);
all_Xmat_shift = all_Xmat_shift(used_inds(full_uinds),use_kInds_up);


for cc = targs
    fprintf('Refitting spk NLs for model %d/%d\n',cc,max(targs));
    cur_Robs = Robs_mat(full_uinds,cc);
    uinds = find(~isnan(cur_Robs));
    if ~isempty(ModData(cc).bestGQM)
        cur_mod = ModData(cc).bestGQM;
        cur_mod.mods(1) = []; %eliminate block filter
        cur_mod = NMMfit_logexp_spkNL(cur_mod,cur_Robs,all_Xmat_shift,[],uinds); %refit spk NL
        
%         %absorb block-by-block offsets into overall spkNL offset param
%         cur_block_filt = cur_mod.mods(1).filtK;
%         cur_used_blocks = ModData(cc).unit_data.used_blocks;
%         poss_used_blocks = ModData(cc).unit_data.poss_used_blocks;
%         cur_used_blocks = find(ismember(cur_used_blocks,poss_used_blocks));
%         cur_mod.spk_NL_params(1) = cur_mod.spk_NL_params(1) + mean(cur_block_filt(cur_used_blocks));
%         cur_mod.mods(1) = [];
        GQM_mod{cc} = cur_mod;
        
        if use_sacMods
        %refit spk NL for saccade model
        [~,~,~,~,~,fgint] = NMMmodel_eval(cur_mod,[],all_Xmat_shift);
        stimG = sum(fgint,2);
        tr_stim{1} = stimG;
        tr_stim{2} = Xmsac; %saccade timing indicator matrix
        tr_stim{3} = reshape(bsxfun(@times,Xmsac,reshape(stimG,[],1)), size(stimG,1),[]);
        cur_sac_mod = sacMod(cc).msac_post_mod;
        cur_sac_mod = NMMfit_logexp_spkNL(cur_sac_mod,cur_Robs,tr_stim,[],uinds);
        sac_mod{cc} = cur_sac_mod;
        end
    else
        GQM_mod{cc} = [];
        
    end
end

%%
% poss_SDs = [0:0.025:0.2 robust_std_dev(base_EP_est*sp_dx)];
poss_SDs = [0:0.05:0.2 robust_std_dev(base_EP_est*sp_dx)];
max_shift = round(full_nPix_us*0.8); %maximum shift size (to avoid going trying to shift more than the n

%create a set of temporally shifted versions of the spiking data
max_tlag = 0; %max time lag for computing autocorrs
tlags = [-max_tlag:max_tlag];

poss_ubins = [1 2 5 10 20 50 100];

ep_rates = nan(length(full_uinds),length(targs));
ep_rates2 = nan(length(full_uinds),length(targs));
ep_rate_cov = nan(length(targs),length(targs),length(tlags),length(poss_SDs));
true_rate_cov = nan(length(targs),length(targs),length(tlags),length(poss_SDs));
ep_spk_cov = nan(length(targs),length(targs),length(poss_ubins),length(poss_SDs));
true_spk_cov = nan(length(targs),length(targs),length(poss_ubins),length(poss_SDs));
% ep_FF_est = nan(length(targs),length(poss_ubins),length(poss_SDs));
% true_FF_est = nan(length(targs),length(poss_ubins),length(poss_SDs));
binned_PSTH_vars = nan(length(targs),length(poss_ubins),length(poss_SDs));
binned_tot_vars = nan(length(targs),length(poss_ubins),length(poss_SDs));
binned_tot_avgs = nan(length(targs),length(poss_ubins),length(poss_SDs));
simRpt_FF_ests = nan(length(targs),length(poss_ubins),length(poss_SDs));

if use_sacMods
Sep_rate_cov = nan(length(targs),length(targs),length(tlags),length(poss_SDs));
Sep_rates = nan(length(full_uinds),length(targs));
Sep_rates2 = nan(length(full_uinds),length(targs));
Strue_rate_cov = nan(length(targs),length(targs),length(tlags),length(poss_SDs));
end

for sd = 1:length(poss_SDs)
    fprintf('SD %d of %d\n',sd,length(poss_SDs));
    
    %rescale EP data to have specified (robust) SD
    cur_EP = base_EP_est./robust_std_dev(base_EP_est)*poss_SDs(sd)/sp_dx;
    fin_shift_cor = round(cur_EP);
    fin_shift_cor(fin_shift_cor > max_shift) = max_shift;
    fin_shift_cor(fin_shift_cor < -max_shift) = -max_shift;
    uset1 = ~inblink_tbt(:);
    
    %RECOMPUTE XMAT
    cur_shift_stimmat_up = all_stimmat_up;
    for i=1:length(full_uinds)
        cur_shift_stimmat_up(used_inds(full_uinds(i)),:) = shift_matrix_Nd(all_stimmat_up(used_inds(full_uinds(i)),:),-fin_shift_cor(i),2);
    end
    all_Xmat_shift = create_time_embedding(cur_shift_stimmat_up,stim_params_full);
    all_Xmat_shift = all_Xmat_shift(used_inds(full_uinds),use_kInds_up);
    
    %compute model-predicted rates given this retinal stim
    for ss = 1:length(targs)
        cc = targs(ss);
        if isstruct(GQM_mod{cc})
            [~,~,ep_rates(:,ss)] = NMMmodel_eval(GQM_mod{cc},[],all_Xmat_shift);
            
            if use_sacMods
            [~,~,~,~,~,fgint] = NMMmodel_eval(sacMod(cc).base_mod,[],all_Xmat_shift);
            stimG = sum(fgint,2);
            tr_stim{1} = stimG;
            tr_stim{2} = Xmsac; %saccade timing indicator matrix
            tr_stim{3} = reshape(bsxfun(@times,Xmsac,reshape(stimG,[],1)), size(stimG,1),[]);
            [~,~,Sep_rates(:,ss)] = NMMmodel_eval(sac_mod{cc},[],tr_stim);
            end
        end
    end
    ep_spks = poissrnd(ep_rates);
    cur_avg_rates = nanmean(ep_rates(uset1,:));
    ep_rates = bsxfun(@minus,ep_rates,cur_avg_rates);
    if use_sacMods
    cur_Savg_rates = nanmean(Sep_rates(uset1,:));
    Sep_rates = bsxfun(@minus,Sep_rates,cur_Savg_rates);
    end
    
    %sample again from the same ET distribution but randomly permute trial
    %assignments
    tperm = randperm(n_utrials);
    cur_EP2 = EP_tbt(:,tperm);
    cur_EP2 = cur_EP2(:)./robust_std_dev(base_EP_est)*poss_SDs(sd)/sp_dx; %rescale to specified SD
    
    cur_inblink = inblink_tbt(:,tperm);
    cur_inblink = cur_inblink(:);
    cur_insac = insac_tbt(:,tperm);
    cur_insac = cur_insac(:);
    uset2 = ~cur_inblink;
    
    if use_sacMods
        cur_Xmsac = Xmsac_tbt(:,tperm,:);
    cur_Xmsac = reshape(cur_Xmsac,[],length(slags));
    end
    
    fin_shift_cor2 = round(cur_EP2);
    fin_shift_cor2(fin_shift_cor2 > max_shift) = max_shift;
    fin_shift_cor2(fin_shift_cor2 < -max_shift) = -max_shift;
    cur_shift_stimmat_up = all_stimmat_up;
    for i=1:length(full_uinds)
        cur_shift_stimmat_up(used_inds(full_uinds(i)),:) = shift_matrix_Nd(all_stimmat_up(used_inds(full_uinds(i)),:),-fin_shift_cor2(i),2);
    end
    all_Xmat_shift = create_time_embedding(cur_shift_stimmat_up,stim_params_full);
    all_Xmat_shift = all_Xmat_shift(used_inds(full_uinds),use_kInds_up);
    
    %compute model-rates for the second version of the retinal stim
    for ss = 1:length(targs)
        cc = targs(ss);
        if ~isempty(GQM_mod{cc})
            [~,~,ep_rates2(:,ss)] = NMMmodel_eval(GQM_mod{cc},[],all_Xmat_shift);
            
            if use_sacMods
            [~,~,~,~,~,fgint] = NMMmodel_eval(sacMod(cc).base_mod,[],all_Xmat_shift);
            stimG = sum(fgint,2);
            tr_stim{1} = stimG;
            tr_stim{2} = cur_Xmsac; %saccade timing indicator matrix
            tr_stim{3} = reshape(bsxfun(@times,cur_Xmsac,reshape(stimG,[],1)), size(stimG,1),[]);
            [~,~,Sep_rates2(:,ss)] = NMMmodel_eval(sac_mod{cc},[],tr_stim);
            end
        end
    end
    ep_spks2 = poissrnd(ep_rates2);
    ep_rates2 = bsxfun(@minus,ep_rates2,nanmean(ep_rates2(uset2,:)));
    if use_sacMods
    Sep_rates2 = bsxfun(@minus,Sep_rates2,nanmean(Sep_rates2(uset2,:)));
    end
    %shiftify the second rate matrix
    mod_rates_shifted = nan(length(full_uinds),length(targs),length(tlags));
    for tt = 1:length(tlags)
        mod_rates_shifted(:,:,tt) = shift_matrix_Nd(ep_rates2,-tlags(tt),1);
    end
    if use_sacMods
    Smod_rates_shifted = nan(length(full_uinds),length(targs),length(tlags));
    for tt = 1:length(tlags)
        Smod_rates_shifted(:,:,tt) = shift_matrix_Nd(Sep_rates2,-tlags(tt),1);
    end
    end
    
    cur_uset = uset1 & uset2;
    %compute EP-contaminated covariances
    for ll = 1:length(tlags)
        ep_rate_cov(:,:,ll,sd) = squeeze(mod_rates_shifted(cur_uset,:,ll))'*ep_rates(cur_uset,:)/length(cur_uset);
        if use_sacMods
        Sep_rate_cov(:,:,ll,sd) = squeeze(Smod_rates_shifted(cur_uset,:,ll))'*Sep_rates(cur_uset,:)/length(cur_uset);
        end
   end
    %compute true covariances
    for ll = 1:length(tlags)
        true_rate_cov(:,:,ll,sd) = squeeze(mod_rates_shifted(cur_uset,:,ll))'*ep_rates2(cur_uset,:)/length(cur_uset);
        if use_sacMods
        Strue_rate_cov(:,:,ll,sd) = squeeze(Smod_rates_shifted(cur_uset,:,ll))'*Sep_rates2(cur_uset,:)/length(cur_uset);
        end
   end
%     
%     cur_ep_rates = bsxfun(@plus,ep_rates,cur_avg_rates); %add back in avg rates
%     cur_ep_rates2 = bsxfun(@plus,ep_rates2,cur_avg_rates);
    cur_ep_rates = ep_rates; %add back in avg rates
    cur_ep_rates2 = ep_rates2;
    ep_spks = reshape(ep_spks,target_uf,n_utrials,length(targs));
    ep_spks2 = reshape(ep_spks2,target_uf,n_utrials,length(targs));
    cur_ep_rates = reshape(cur_ep_rates,target_uf,n_utrials,length(targs));
    cur_ep_rates2 = reshape(cur_ep_rates2,target_uf,n_utrials,length(targs));
    ep_uset = reshape(uset1,target_uf,n_utrials) & reshape(uset2,target_uf,n_utrials);
    
    for pp = 1:length(poss_ubins)
        bin_usfac = poss_ubins(pp);
        n_newbins = floor(target_uf/bin_usfac);
        new_ep_spks = zeros(n_newbins,n_utrials,length(targs));
        new_ep_spks2 = zeros(n_newbins,n_utrials,length(targs));
%         new_ep_avgrates = zeros(n_newbins,n_utrials,length(targs));
        new_ep_rates = zeros(n_newbins,n_utrials,length(targs));
        new_ep_rates2 = zeros(n_newbins,n_utrials,length(targs));
        new_include = zeros(n_newbins,n_utrials);
        for ii = 1:bin_usfac
            new_ep_spks = new_ep_spks + ep_spks(ii:bin_usfac:(ii+bin_usfac*(n_newbins-1)),:,:);
            new_ep_spks2 = new_ep_spks2 + ep_spks2(ii:bin_usfac:(ii+bin_usfac*(n_newbins-1)),:,:);
            new_ep_rates = new_ep_rates + cur_ep_rates(ii:bin_usfac:(ii+bin_usfac*(n_newbins-1)),:,:);
            new_ep_rates2 = new_ep_rates2 + cur_ep_rates2(ii:bin_usfac:(ii+bin_usfac*(n_newbins-1)),:,:);
            new_include = new_include + ep_uset(ii:bin_usfac:(ii+bin_usfac*(n_newbins-1)),:);
        end
        new_ep_rates = reshape(new_ep_rates,[],length(targs));
        new_ep_rates2 = reshape(new_ep_rates2,[],length(targs));
        new_ep_avgrates = cur_avg_rates*bin_usfac;
        new_ep_spks = reshape(new_ep_spks,[],length(targs));
        new_ep_spks2 = reshape(new_ep_spks2,[],length(targs));
        new_include = logical(new_include(:));
        
        %spk count variance for a rate-modulated PP is (across-trial rate
        %var + mean rate)
        binned_PSTH_vars(:,pp,sd) = nanmean(new_ep_rates(new_include,:).*new_ep_rates2(new_include,:));
        binned_tot_vars(:,pp,sd) = nanvar(new_ep_rates(new_include,:));
%         ep_FF_est(:,pp,sd) = (cur_rate_vars - cur_PSTH_vars + new_ep_avgrates)./new_ep_avgrates;
        binned_tot_avgs(:,pp,sd) = new_ep_avgrates;
        
        new_ep_spks = bsxfun(@minus,new_ep_spks,nanmean(new_ep_spks(new_include,:)));
        new_ep_spks2 = bsxfun(@minus,new_ep_spks2,nanmean(new_ep_spks2(new_include,:)));
        ep_spk_cov(:,:,pp,sd) = squeeze(new_ep_spks(new_include,:)'*new_ep_spks2(new_include,:))/sum(new_include);
        true_spk_cov(:,:,pp,sd) = squeeze(new_ep_spks(new_include,:)'*new_ep_spks(new_include,:))/sum(new_include);
    end
    
    
    N_FF_simtrials = 10;
    N_sim_rpts = 50;
    cur_full_uinds = full_uinds(1:N_FF_simtrials*target_uf);
    base_stim = all_stimmat_up(1:max(used_inds(cur_full_uinds)),:);
    rpt_rate_sim = nan(length(cur_full_uinds),N_sim_rpts,length(targs));
    for rr = 1:N_sim_rpts
        
        tperm = randperm(n_utrials);
        cur_EP2 = EP_tbt(:,tperm);
        cur_EP2 = cur_EP2(cur_full_uinds)./robust_std_dev(base_EP_est)*poss_SDs(sd)/sp_dx; %rescale to specified SD
        
        cur_inblink = inblink_tbt(:,tperm);
        cur_inblink = cur_inblink(cur_full_uinds);
        cur_insac = insac_tbt(:,tperm);
        cur_insac = cur_insac(cur_full_uinds);
        cur_uset = ~cur_inblink;
        
        fin_shift_cor2 = round(cur_EP2);
        fin_shift_cor2(fin_shift_cor2 > max_shift) = max_shift;
        fin_shift_cor2(fin_shift_cor2 < -max_shift) = -max_shift;
        cur_shift_stimmat_up = base_stim;
        for i=1:length(cur_full_uinds)
            cur_shift_stimmat_up(used_inds(cur_full_uinds(i)),:) = shift_matrix_Nd(base_stim(used_inds(cur_full_uinds(i)),:),-fin_shift_cor2(i),2);
        end
        all_Xmat_shift = create_time_embedding(cur_shift_stimmat_up,stim_params_full);
        all_Xmat_shift = all_Xmat_shift(used_inds(cur_full_uinds),use_kInds_up);
        
        %compute model-rates for the second version of the retinal stim
        for ss = 1:length(targs)
            cc = targs(ss);
            if ~isempty(GQM_mod{cc})
                [~,~,rpt_rate_sim(:,rr,ss)] = NMMmodel_eval(GQM_mod{cc},[],all_Xmat_shift);
            end
        end
        rpt_rate_sim(~cur_uset,:,:) = nan;
    end
    rpt_rate_sim = reshape(rpt_rate_sim,target_uf,N_FF_simtrials,N_sim_rpts,length(targs));
    
    for pp = 1:length(poss_ubins)
        bin_usfac = poss_ubins(pp);
        n_newbins = floor(target_uf/bin_usfac);
        
        new_ep_rates = nan(n_newbins,N_FF_simtrials,N_sim_rpts,length(targs),bin_usfac);
        for ii = 1:bin_usfac
           new_ep_rates(:,:,:,:,ii) = rpt_rate_sim(ii:bin_usfac:(ii+bin_usfac*(n_newbins-1)),:,:,:); 
        end
        new_ep_rates = squeeze(nanmean(new_ep_rates,5))*bin_usfac;
        
        acrossTrial_var = squeeze(nanvar(new_ep_rates,[],3));
        acrossTrial_avg = squeeze(nanmean(new_ep_rates,3));
        all_FF_ests = (acrossTrial_avg + acrossTrial_var)./acrossTrial_avg;
        simRpt_FF_ests(:,pp,sd) = squeeze(nanmean(reshape(all_FF_ests,[],length(targs))));
    end
end

%%
ep_noise_cov = true_rate_cov - ep_rate_cov;

true_sig_vars = nan(length(targs),length(poss_SDs));
ep_alpha_funs = nan(length(targs),length(poss_SDs));
for ss = 1:length(targs)
    true_sig_vars(ss,:) = squeeze(true_rate_cov(ss,ss,max_tlag+1,:));
    ep_alpha_funs(ss,:) = squeeze(ep_rate_cov(ss,ss,max_tlag+1,:))./true_sig_vars(ss,:)';
end
if use_sacMods
Strue_sig_vars = nan(length(targs),length(poss_SDs));
Sep_alpha_funs = nan(length(targs),length(poss_SDs));
for ss = 1:length(targs)
    Strue_sig_vars(ss,:) = squeeze(Strue_rate_cov(ss,ss,max_tlag+1,:));
    Sep_alpha_funs(ss,:) = squeeze(Sep_rate_cov(ss,ss,max_tlag+1,:))./Strue_sig_vars(ss,:)';
end
end

uset = full_uinds(~in_blink_inds(full_uinds));
ov_avg_rates = nanmean(Robs_mat(uset,targs));

%normalize to get measured correlations
ep_sig_corr = nan(size(true_rate_cov));
ep_noise_corr = nan(size(true_rate_cov));
ep_psth_corr = nan(size(true_rate_cov));
for ss = 1:length(poss_SDs)
    corr_norm = sqrt(true_sig_vars(:,ss)*true_sig_vars(:,ss)');
    ep_sig_corr(:,:,:,ss) = bsxfun(@rdivide,true_rate_cov(:,:,:,ss),corr_norm);
    ep_psth_corr(:,:,:,ss) = bsxfun(@rdivide,ep_rate_cov(:,:,:,ss),corr_norm);
    ep_noise_corr(:,:,:,ss) = bsxfun(@rdivide,ep_noise_cov(:,:,:,ss),corr_norm);
end

%%
ep_noise_cov_spk = true_spk_cov - ep_spk_cov;

true_spk_vars = nan(length(targs),length(poss_ubins));
for ss = 1:length(targs)
    true_spk_vars(ss,:) = squeeze(true_spk_cov(ss,ss,:,end));
end
ep_noise_corr_spk = nan(size(true_spk_cov));
for sd = 1:length(poss_SDs)
    for pp = 1:length(poss_ubins)
        corr_norm = sqrt(true_spk_vars(:,pp)*true_spk_vars(:,pp)');
        ep_noise_corr_spk(:,:,pp,sd) = squeeze(ep_noise_cov_spk(:,:,pp,sd))./corr_norm;
    end
end

%% COMPUTE EP AUTOCORR FUNCTION
tbt_EP_data = nan(n_trials,400-(params.beg_buffer+params.end_buffer)/params.dt);
for ii = 1:n_trials
   cur_inds = find(all_trialvec(used_inds) == ii);
   tbt_EP_data(ii,1:length(cur_inds)) = fin_tot_corr(cur_inds);
end
tbt_EP_data = (tbt_EP_data - nanmean(tbt_EP_data(:)))/nanstd(tbt_EP_data(:));

max_lag = 50;
EP_lags = 0:max_lag;

z_pad_tbt = cat(2,tbt_EP_data,nan(n_trials,max_lag));

EP_xcov = nan(length(EP_lags),1);
for tt = 1:length(EP_lags)
   shifted = shift_matrix_Nd(z_pad_tbt,EP_lags(tt),2); 
   EP_xcov(tt) = nanmean(z_pad_tbt(:).*shifted(:));
end

ov_EP_data.EP_xcov = EP_xcov;
ov_EP_data.EP_lags = EP_lags;

%%
for ss = 1:length(targs)
    cc = targs(ss);
    EP_data(cc).ModData = ModData(cc);
    EP_data(cc).poss_SDs = poss_SDs;
    EP_data(cc).poss_ubins = poss_ubins;
    EP_data(cc).base_vars = true_sig_vars(ss,:);
    EP_data(cc).alpha_funs = ep_alpha_funs(ss,:);
    EP_data(cc).sig_corr_mat = squeeze(ep_sig_corr(ss,:,:,:));
    EP_data(cc).psth_corr_mat = squeeze(ep_psth_corr(ss,:,:,:));
    EP_data(cc).noise_corr_mat = squeeze(ep_noise_corr(ss,:,:,:));
    EP_data(cc).ep_noise_cov = squeeze(ep_noise_cov(ss,:,:,:));
    EP_data(cc).true_rate_cov = squeeze(true_rate_cov(ss,:,:,:));
    EP_data(cc).ep_rate_cov = squeeze(ep_rate_cov(ss,:,:,:));
%     EP_data(cc).ep_FF_est = squeeze(ep_FF_est(ss,:,:));
    
    EP_data(cc).simrpt_FF = squeeze(simRpt_FF_ests(ss,:,:));

    EP_data(cc).bin_PSTH_vars = squeeze(binned_PSTH_vars(ss,:,:));
    EP_data(cc).bin_tot_vars = squeeze(binned_tot_vars(ss,:,:));
    EP_data(cc).bin_tot_avgs = squeeze(binned_tot_avgs(ss,:,:));

    EP_data(cc).spk_vars = true_spk_vars(ss,:);
    EP_data(cc).noise_corr_spk = squeeze(ep_noise_corr_spk(ss,:,:,:));
end

%%
anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
if ~exist(anal_dir)
    mkdir(anal_dir)
end
cd(anal_dir);

sname = [sname sprintf('_ori%d',bar_ori)];
if rec_number > 1
    sname = strcat(sname,sprintf('_r%d',rec_number));
end

save(sname,'targs','EP_data','use_MUA','fit_unCor','ov_EP_data');
