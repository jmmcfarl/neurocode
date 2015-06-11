% clear all
% close all

addpath('~/other_code/fastBSpline/');

global Expt_name bar_ori monk_name rec_type rec_number

% Expt_name = 'M296';
% monk_name = 'lem';
% bar_ori = 45; %bar orientation to use (only for UA recs)
% rec_number = 1;

use_hres_ET = true; EP_params.use_hres_ET = use_hres_ET; %use high-res eye-tracking?
exclude_sacs = true; EP_params.exclude_sacs = exclude_sacs;
sub_trialavgs = true; EP_params.sub_trialavgs = sub_trialavgs; %subtract out trial avg spike counts

poss_us_change = [0.5 1 2 3 4]; %possible scaling factors to apply to the original stimulus bar size

use_LOOXV = 1; %[0 is no LOO; 1 is SUs only; 2 is SU + MU]

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

ov_RF_pos = Expts{expt_data.used_blocks(1)}.Stimvals.rf(1:2)/params.scale_fac;

sname = 'rpt_variability_barsize_sim';
%%
Expt_num = str2num(Expt_name(2:end));

cd(data_dir);
% load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
model_dir = ['~/Analysis/bruce/' Expt_name '/models'];
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';
mod_name = 'corrected_models_comp';

if rec_number > 1
    cluster_dir = [cluster_dir sprintf('/rec%d',rec_number)];
end

%if using coil info
if any(params.use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end
et_hres_anal_name = strcat(et_anal_name,'_hres');

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
et_hres_anal_name = [et_hres_anal_name sprintf('_ori%d',bar_ori)];
mod_name = [mod_name sprintf('_ori%d',bar_ori)];

if rec_number > 1
    mod_name = strcat(mod_name,sprintf('_r%d',rec_number));
    et_mod_data_name = strcat(et_mod_data_name,sprintf('r%d',rec_number));
    et_hres_anal_name = strcat(et_hres_anal_name,sprintf('r%d',rec_number));
    et_anal_name = strcat(et_anal_name,sprintf('r%d',rec_number));
end

% et_hres_anal_name = strcat(et_hres_anal_name,'_fullLOO');

%% LOAD EYE-TRACKING DATA
cd(et_dir)
load(et_mod_data_name,'all_mod*');
if use_hres_ET
    fprintf('Loading ET data %s\n',et_hres_anal_name);
    load(et_hres_anal_name)
else
    fprintf('Loading ET data %s\n',et_anal_name);
    load(et_anal_name);
end
tr_set = et_tr_set;

%for some recs the saccades for ET_data were detected using a slightly
%different algo. Eliminate the saccades causing the difference
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

%% LOAD IN MODEL FITS
load([model_dir '/' mod_name]);

%% BIN SPIKES FOR MU AND SU
all_binned_mua = spikes_int82double(spike_data.binned_mua);
all_binned_sua = spikes_int82double(spike_data.binned_sua);
Clust_data = spike_data.Clust_data;
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

n_units = size(all_binned_mua,2) + size(all_binned_sua,2);

%% extract time series to keep track of trial number and block number
orig_dt = params.dt;

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

saccades = ET_data.saccades(used_saccade_set);

sac_durs = [saccades(:).duration];
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

sacburst_set = find([saccades(:).isi] < params.sac_burst_isi | [saccades(:).next_isi] < params.sac_burst_isi);
micro_sacs = find([saccades(:).amplitude] < params.micro_thresh & ~used_is_blink');

msac_bursts = micro_sacs(ismember(micro_sacs,sacburst_set));
% micro_sacs(ismember(micro_sacs,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'

%guided saccades are those whose parallel component is large enough and
%that aren't blinks (and whose duration is not too long to be suspicious
big_sacs = find(abs(sac_deltaX) > params.gsac_thresh & ~used_is_blink' & sac_durs <= params.max_gsac_dur);
used_is_blink(sac_durs >= params.max_gsac_dur) = true;

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
%set of units where we have LOOXV on eye-tracking
if use_LOOXV == 2
    loo_set = tr_set;
elseif use_LOOXV == 1
    loo_set = tr_set(all_mod_SU(tr_set) > 0);
else
    loo_set = [];
end

%for model fitting
targs = setdiff(1:n_units,1:params.n_probes); %SU only
targs(targs > length(ModData)) = [];

sac_shift = et_params.sac_shift; %forward projection of saccade start times
if use_hres_ET %if using high-res ET
    [post_mean_EP,post_std_EP] = construct_eye_position(best_fix_cor,best_fix_std,...
        drift_post_mean,drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    sp_dx = et_params.sp_dx;
    post_mean_EP = post_mean_EP*sp_dx;
    
else %if using base resolution ET
    [post_mean_EP,post_std_EP] = construct_eye_position(it_fix_post_mean(end,:),it_fix_post_std(end,:),...
        drift_post_mean(end,:),drift_post_std(end,:),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    sp_dx = et_params.sp_dx;
    post_mean_EP = post_mean_EP*sp_dx;
end

%% IDENTIFY TIMES WITHIN SACCADES AND BLINKS
sac_buff = round(0.05/params.dt); EP_params.sac_buff = sac_buff; %window of data to exclude during saccades
sac_delay = round(0.03/params.dt); EP_params.sac_delay = sac_delay; %shift exclusion window to account for neural delay
blink_buff = round(0.1/params.dt); EP_params.blink_buff = blink_buff; %window for excluding blinks

nf = 400;
used_nf = nf-(params.beg_buffer + params.end_buffer)/params.dt;

in_sac_inds = false(NT,1);
nblink_start_inds = saccade_start_inds(~used_is_blink);
nblink_stop_inds = saccade_stop_inds(~used_is_blink);
for ii = 1:length(nblink_start_inds)
    cur_inds = (nblink_start_inds(ii):(nblink_stop_inds(ii) + sac_buff)) + sac_delay;
    cur_inds(cur_inds > NT) = [];
    in_sac_inds(cur_inds) = true;
end

in_blink_inds = false(NT,1);
blink_start_inds = saccade_start_inds(used_is_blink);
blink_stop_inds = saccade_stop_inds(used_is_blink);
for ii = 1:length(blink_start_inds)
    cur_inds = (blink_start_inds(ii):(blink_stop_inds(ii) + blink_buff)) + sac_delay;
    cur_inds(cur_inds > NT) = [];
    in_blink_inds(cur_inds) = true;
end

in_sac_inds(isnan(in_sac_inds)) = 0; in_blink_inds(isnan(in_blink_inds)) = 0;
in_sac_inds = logical(in_sac_inds); in_blink_inds = logical(in_blink_inds);


%% PROCESS MODEL FITS
has_stim_mod = false(length(targs),1);
for cc = 1:length(targs)
    if ~isempty(ModData(targs(cc)).bestGQM)
        cur_mod = ModData(targs(cc)).bestGQM;
        cur_block_filt = cur_mod.mods(1).filtK;
        cur_used_blocks = ModData(targs(cc)).unit_data.used_blocks;
        poss_used_blocks = ModData(targs(cc)).unit_data.poss_used_blocks;
        cur_used_blocks = find(ismember(cur_used_blocks,poss_used_blocks));
        cur_mod.spk_NL_params(1) = cur_mod.spk_NL_params(1) + mean(cur_block_filt(cur_used_blocks));
        cur_mod.mods(1) = []; %eliminate block filter
        stim_mod(cc) = cur_mod;
        has_stim_mod(cc) = true;
    end
end

%% use only completed trials
use_trials = unique(all_trialvec(used_inds)); %set of potentially usable trials

%use only completed trials
target_uf = 400 - (params.beg_buffer + params.end_buffer)/params.dt;
T = tabulate(all_trialvec(used_inds));
complete_trials = find(T(:,2) >= target_uf);
use_trials(~ismember(use_trials,complete_trials)) = [];

full_uinds = find(ismember(all_trialvec(used_inds),use_trials));
n_used_utrials = length(full_uinds)/target_uf;

%matrix of trial-by-trial EP estimates
fin_shift_cor = round(post_mean_EP(full_uinds)/modFitParams.sp_dx);
EP_tbt = reshape(fin_shift_cor,target_uf,n_used_utrials);

%TBT mats for sac and blink indicators
inblink_tbt = reshape(in_blink_inds(full_uinds),target_uf,n_used_utrials);
insac_tbt = reshape(in_sac_inds(full_uinds),target_uf,n_used_utrials);

%% generate stimulus

%select subset of pixels used for model fitting
buffer_pix = floor((params.full_nPix - modFitParams.use_nPix_us/modFitParams.spatial_usfac)/2);
[Xinds_up,~] = meshgrid(1/modFitParams.spatial_usfac:1/modFitParams.spatial_usfac:params.full_nPix,1:modFitParams.flen);
cur_use_pix = (1/modFitParams.spatial_usfac:1/modFitParams.spatial_usfac:(modFitParams.use_nPix_us/modFitParams.spatial_usfac)) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

%baselin stimulus properties
base_spatial_usfac = modFitParams.spatial_usfac;
base_nPix = params.full_nPix;
dds = mode(expt_data.expt_dds);

for uu = 1:length(poss_us_change)
    fprintf('US change %d of %d\n',uu,length(poss_us_change));

cur_us_change = poss_us_change(uu);

cur_usfac = cur_us_change*base_spatial_usfac; %current spatial up-sampling factor
cur_nPix = base_nPix/cur_us_change; %number of bars per frame

%make a random bar pattern stim
rand_Stim = zeros(stimComp.NT,cur_nPix);
nzero_vals = rand(stimComp.NT,cur_nPix) < dds/100;
rand_vals = rand(stimComp.NT,cur_nPix); rand_vals(rand_vals > 0.5) = 1; rand_vals(rand_vals < 0.5) = -1;
rand_Stim(nzero_vals) = rand_vals(nzero_vals);

%spatial up-sampling of the stimulus
full_nPix_us = cur_usfac*cur_nPix;
if cur_usfac > 1
    rand_Stim_up = zeros(size(rand_Stim,1),full_nPix_us);
    for ii = 1:size(rand_Stim,2)
        for jj = 1:cur_usfac
            rand_Stim_up(:,cur_usfac*(ii-1)+jj) = rand_Stim(:,ii);
        end
    end
elseif cur_usfac == 1
    rand_Stim_up = rand_Stim;
end

%%
stim_params_full = NMMcreate_stim_params([modFitParams.flen full_nPix_us]);
%COMPUTE XMAT
cur_shift_stimmat_up = rand_Stim_up;
for ii=1:length(full_uinds)
    cur_shift_stimmat_up(used_inds(full_uinds(ii)),:) = shift_matrix_Nd(cur_shift_stimmat_up(used_inds(full_uinds(ii)),:),-fin_shift_cor(ii),2);
end
all_Xmat_shift = create_time_embedding(cur_shift_stimmat_up,stim_params_full);
all_Xmat_shift = all_Xmat_shift(used_inds(full_uinds),use_kInds_up);


%% simulate neuron responses to this stimulus 
mod_rates1 = nan(length(full_uinds),length(targs));
for cc = 1:length(targs)
    if has_stim_mod(cc)
        [~,~,mod_rates1(:,cc)] = NMMmodel_eval(stim_mod(cc),[],all_Xmat_shift);
    end
end

%% permute the tbt eye position data, and compute a new version of the retinal stimulu
tperm = randperm(n_used_utrials);
cur_EP2 = EP_tbt(:,tperm);
cur_ep_cor = cur_EP2(:);

%COMPUTE XMAT
cur_shift_stimmat_up = rand_Stim_up;
for ii=1:length(full_uinds)
    cur_shift_stimmat_up(used_inds(full_uinds(ii)),:) = shift_matrix_Nd(cur_shift_stimmat_up(used_inds(full_uinds(ii)),:),-cur_ep_cor(ii),2);
end
all_Xmat_shift = create_time_embedding(cur_shift_stimmat_up,stim_params_full);
all_Xmat_shift = all_Xmat_shift(used_inds(full_uinds),use_kInds_up);

%% simulate neuron responses to this set of 'trials'
mod_rates2 = nan(length(full_uinds),length(targs));
for cc = 1:length(targs)
    if has_stim_mod(cc)
        [~,~,mod_rates2(:,cc)] = NMMmodel_eval(stim_mod(cc),[],all_Xmat_shift);
    end
end

%% compute psth covariance and overall rate covariance
mod_rates1(inblink_tbt(:),:) = nan;
mod_rates2(inblink_tbt(:),:) = nan;
if exclude_sacs
    mod_rates1(insac_tbt(:),:) = nan;
    mod_rates2(insac_tbt(:),:) = nan;
end
mod_rates1 = bsxfun(@minus,mod_rates1,nanmean(mod_rates1));
mod_rates2 = bsxfun(@minus,mod_rates2,nanmean(mod_rates2));

rate_covs(uu,:) = nanvar(mod_rates1);
psth_covs(uu,:) = nanmean(bsxfun(@times,mod_rates1,mod_rates2));
end

alphas = psth_covs./rate_covs; %alphas

%% store values for each cell
for cc = 1:length(targs)
    EP_sim_data(cc).rate_covs = rate_covs(:,cc);
    EP_sim_data(cc).psth_covs = psth_covs(:,cc);
    EP_sim_data(cc).alphas = alphas(:,cc);
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

save(sname,'targs','EP_sim_data','poss_us_change');

