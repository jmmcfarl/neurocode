% clear all
% close all

addpath('~/other_code/fastBSpline/');

global Expt_name bar_ori monk_name rec_type rec_number

% Expt_name = 'M011';
% monk_name = 'jbe';
% bar_ori = 160; %bar orientation to use (only for UA recs)
% rec_number = 1;
%
% [266-80 270-60 275-135 277-70 281-140 287-90 289-160 294-40 296-45 297-0/90 5-50 9-0 10-60 11-160 12-0 13-100 14-40 320-100]

sname = 'sim_variability_compact_FIN';

et_mod_data_name = 'full_eyetrack_initmods_FIN2_Rinit';
et_anal_name = 'full_eyetrack_FIN2_Rinit';
mod_name = 'corrected_models_comp_FIN2';

use_MUA = false; sim_params.use_MUA = use_MUA; %use MUA in model-fitting
use_hres_ET = true; sim_params.use_hres_ET = use_hres_ET; %use high-res eye-tracking?
exclude_sacs = false; sim_params.exclude_sacs = exclude_sacs; %exclude data surrounding microsaccades?
poss_SDs = [0.05:0.025:0.2];%range of possible EP SDs to test (last is the empirical)
poss_ubins = [1 2 5 10 20 50 100]; sim_params.poss_ubins = poss_ubins; %range of temporal downsampling factors to test

do_xcorr = true;
use_LOOXV = 1; %[0 is no LOO; 1 is SUs only; 2 is SU + MU]

%%
data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
if ~exist(anal_dir,'dir')
    mkdir(anal_dir)
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

%directories
et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
model_dir = ['~/Analysis/bruce/' Expt_name '/models'];

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
tr_set = et_tr_set; %set of units used in ET

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

%% EXTRACT RELEVANT SACCADE DATA
corrected_eye_vals_interp = ET_data.interp_eye_pos;
sac_start_times = [ET_data.saccades(:).start_time];
sac_stop_times = [ET_data.saccades(:).stop_time];

interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));

used_saccade_set = find(ismember(interp_sac_start_inds,used_inds)); %saccades occuring during used data

saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));

%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';
saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);

saccades = ET_data.saccades(used_saccade_set);
saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

used_is_blink = ET_data.is_blink(used_saccade_set); %which used events are blinks

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
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink)); %is this fixation following a blink
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
if use_MUA
    targs = 1:n_units; %SU and MU
else
    targs = setdiff(1:n_units,1:params.n_probes); %SU only
end
targs(targs > length(ModData)) = [];

sac_shift = et_params.sac_shift; %forward projection of saccade start times
if use_hres_ET %if using high-res ET
    [post_mean_EP,post_std_EP] = construct_eye_position(best_fix_cor,best_fix_std,...
        drift_post_mean,drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    sp_dx = et_params.sp_dx;
    post_mean_EP = post_mean_EP*sp_dx; %in deg
else %if using base resolution ET
    [post_mean_EP,post_std_EP] = construct_eye_position(it_fix_post_mean(end,:),it_fix_post_std(end,:),...
        drift_post_mean(end,:),drift_post_std(end,:),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    sp_dx = et_params.sp_dx;
    post_mean_EP = post_mean_EP*sp_dx;
end

%% get set of trials used for calcs, and construct TBT matrices
all_trial_Se = [trial_data(:).se];
n_rpt_seeds = length(params.rpt_seeds); EP_params.n_rpt_seeds = n_rpt_seeds; %seed values for repeat sequences

all_rpt_trials = find(ismember(all_trial_Se,params.rpt_seeds)); %trial indices of repeat trials

use_trials = unique(all_trialvec(used_inds)); %set of potentially usable trials
use_trials(ismember(use_trials,all_rpt_trials)) = []; %dont use repeat trials

%use only completed trials
target_uf = 400 - (params.beg_buffer + params.end_buffer)/params.dt;
T = tabulate(all_trialvec(used_inds));
complete_trials = find(T(:,2) >= target_uf);
use_trials(~ismember(use_trials,complete_trials)) = [];

full_uinds = find(ismember(all_trialvec(used_inds),use_trials));
n_utrials = length(full_uinds)/target_uf;
if mod(n_utrials,1) ~= 0
    error('non-integer number of trials');
end

%% IDENTIFY TIMES WITHIN SACCADES AND BLINKS
sac_buff = round(0.05/params.dt); EP_params.sac_buff = sac_buff; %window of data to exclude during saccades
sac_delay = round(0.03/params.dt); EP_params.sac_delay = sac_delay; %shift exclusion window to account for neural delay
blink_buff = round(0.1/params.dt); EP_params.blink_buff = blink_buff; %window for excluding blinks

%these are the non-blink saccades (real saccades)
nblink_start_inds = saccade_start_inds(~used_is_blink);
nblink_stop_inds = saccade_stop_inds(~used_is_blink);
in_sac_inds = false(NT,1); %store indices that count as 'affected' by a saccade
for ii = 1:length(nblink_start_inds)
    cur_inds = (nblink_start_inds(ii):(nblink_stop_inds(ii) + sac_buff)) + sac_delay;
    cur_inds(cur_inds > NT) = [];
    in_sac_inds(cur_inds) = true;
end

%find indices that count as affected by a blink
blink_start_inds = saccade_start_inds(used_is_blink);
blink_stop_inds = saccade_stop_inds(used_is_blink);
in_blink_inds = false(NT,1);
for ii = 1:length(blink_start_inds)
    cur_inds = (blink_start_inds(ii):(blink_stop_inds(ii) + blink_buff)) + sac_delay;
    cur_inds(cur_inds > NT) = [];
    in_blink_inds(cur_inds) = true;
end

%reshape these vectors into trial-by-trial matrices
in_sac_inds = in_sac_inds(full_uinds);
in_blink_inds = in_blink_inds(full_uinds);

%% PROCESS MODEL FITS

has_stim_mod = false(length(targs),1);
for cc = 1:length(targs) %loop over units used in analysis
    if ~isempty(ModData(targs(cc)).bestGQM) %if we have a model for this unit
        cur_mod = ModData(targs(cc)).bestGQM; %use the optimized GQM
        
        %remove the block-by-block variability in the model by
        %incorporating the avg output of the block-filter
        cur_block_filt = cur_mod.mods(1).filtK; %filter applied to the block index
        cur_used_blocks = ModData(targs(cc)).unit_data.used_blocks; %which blocks was this neuron isolated during
        poss_used_blocks = ModData(targs(cc)).unit_data.poss_used_blocks; %total set of used blocks
        cur_used_blocks = find(ismember(poss_used_blocks,cur_used_blocks)); %indices of blocks where this neuron was isolated
        cur_mod.spk_NL_params(1) = cur_mod.spk_NL_params(1) + mean(cur_block_filt(cur_used_blocks)); %add the average output of the block-filter to the spkNL offset
        cur_mod.mods(1) = []; %eliminate block filter
        
        %store model data
        stim_mod(cc) = cur_mod;
        EP_data(cc,1).unit_data = ModData(targs(cc)).unit_data;
        EP_data(cc,1).tune_props = ModData(targs(cc)).tune_props;
        EP_data(cc,1).bestGQM = ModData(targs(cc)).bestGQM;
        EP_data(cc,1).nullMod = ModData(targs(cc)).nullMod;
        EP_data(cc,1).useMod = cur_mod;
        has_stim_mod(cc) = true;
    end
end

%% get stimulus and eye position shifts
all_stim_mat = decompressTernNoise(stimComp);

%spatial up-sampling of the stimulus
full_nPix_us = modFitParams.spatial_usfac*params.full_nPix;
if modFitParams.spatial_usfac > 1
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:modFitParams.spatial_usfac
            all_stimmat_up(:,modFitParams.spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
        end
    end
elseif modFitParams.spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
end

%select subset of pixels used for model fitting
buffer_pix = floor((params.full_nPix - modFitParams.use_nPix_us/modFitParams.spatial_usfac)/2);
[Xinds_up,~] = meshgrid(1/modFitParams.spatial_usfac:1/modFitParams.spatial_usfac:params.full_nPix,1:modFitParams.flen);
cur_use_pix = (1/modFitParams.spatial_usfac:1/modFitParams.spatial_usfac:(modFitParams.use_nPix_us/modFitParams.spatial_usfac)) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

stim_params_full = NMMcreate_stim_params([modFitParams.flen full_nPix_us]);

EP_tbt = reshape(post_mean_EP(full_uinds),[],n_utrials);
base_EP_SD = robust_std_dev(EP_tbt(:));

inblink_tbt = reshape(in_blink_inds,[],n_utrials);
insac_tbt = reshape(in_sac_inds,[],n_utrials);
%%
if do_xcorr
    [CI,CJ] = meshgrid(1:length(targs));
    un_pairs = CJ >= CI; %all unique pairs of neurons (including self-pairs)
    Cpairs = [CI(un_pairs) CJ(un_pairs)]; %I and J indices for each pair
    n_cell_pairs = size(Cpairs,1);    
end
%%
poss_SDs = [poss_SDs base_EP_SD]; %make the last EP SD to test the empirical one
sim_params.poss_SDs = poss_SDs; 
max_shift = round(full_nPix_us*0.8); %maximum shift size (to avoid going trying to shift more than the n

sim_psth_vars = nan(length(poss_SDs),length(poss_ubins),length(targs));
sim_tot_vars = nan(length(poss_SDs),length(poss_ubins),length(targs));
sim_mean_rates = nan(length(poss_SDs),length(poss_ubins),length(targs));

sim_psth_covars = nan(length(poss_SDs),length(poss_ubins),length(targs),length(targs));
sim_tot_covars = nan(length(poss_SDs),length(poss_ubins),length(targs),length(targs));

ep_rates = nan(length(full_uinds),length(targs));
ep_rates2 = nan(length(full_uinds),length(targs));
for sd = 1:length(poss_SDs)
    fprintf('SD %d of %d\n',sd,length(poss_SDs));
    
    %rescale EP data to have specified (robust) SD
    cur_EP = EP_tbt./base_EP_SD*poss_SDs(sd);
    fin_shift_cor = round(cur_EP(:)/modFitParams.sp_dx);
    fin_shift_cor(fin_shift_cor > max_shift) = max_shift;
    fin_shift_cor(fin_shift_cor < -max_shift) = -max_shift;
    
    %RECOMPUTE XMAT
    cur_shift_stimmat_up = all_stimmat_up;
    for i=1:length(full_uinds)
        cur_shift_stimmat_up(used_inds(full_uinds(i)),:) = shift_matrix_Nd(cur_shift_stimmat_up(used_inds(full_uinds(i)),:),-fin_shift_cor(i),2);
    end
    all_Xmat_shift = create_time_embedding(cur_shift_stimmat_up,stim_params_full);
    all_Xmat_shift = all_Xmat_shift(used_inds(full_uinds),use_kInds_up);
    
    %if using tent-basis spatial-upsampling
    if modFitParams.add_usfac > 1
        all_Xmat_shift = tb_proc_stim(all_Xmat_shift,modFitParams.add_usfac,modFitParams.flen);
    end
    
    %compute model-predicted rates given this retinal stim
    for cc = 1:length(targs)
        if has_stim_mod(cc)
            [~,~,ep_rates(:,cc)] = NMMmodel_eval(stim_mod(cc),[],all_Xmat_shift);
        end
    end
    uset1 = ~inblink_tbt;
    if exclude_sacs
       uset1 = uset1 & ~insac_tbt; 
    end
    %subtract off avg rate
    cur_avg_rates = nanmean(ep_rates(uset1(:),:));
    ep_rates = bsxfun(@minus,ep_rates,cur_avg_rates);
    
    %sample again from the same ET distribution but randomly permute trial
    %assignments
    tperm = randperm(n_utrials);
    cur_EP2 = EP_tbt(:,tperm);
    cur_EP2 = cur_EP2(:)./base_EP_SD*poss_SDs(sd); %rescale to specified SD
    
    cur_inblink = inblink_tbt(:,tperm);
    cur_insac = insac_tbt(:,tperm);
    uset2 = ~cur_inblink;
    if exclude_sacs
       uset2 = uset2 & ~inblink_tbt; 
    end
    fin_shift_cor = round(cur_EP2(:)/modFitParams.sp_dx);
    fin_shift_cor(fin_shift_cor > max_shift) = max_shift;
    fin_shift_cor(fin_shift_cor < -max_shift) = -max_shift;
    
    cur_shift_stimmat_up = all_stimmat_up;
    for i=1:length(full_uinds)
        cur_shift_stimmat_up(used_inds(full_uinds(i)),:) = shift_matrix_Nd(all_stimmat_up(used_inds(full_uinds(i)),:),-fin_shift_cor(i),2);
    end
    all_Xmat_shift = create_time_embedding(cur_shift_stimmat_up,stim_params_full);
    all_Xmat_shift = all_Xmat_shift(used_inds(full_uinds),use_kInds_up);
    
    %if using tent-basis spatial-upsampling
    if modFitParams.add_usfac > 1
        all_Xmat_shift = tb_proc_stim(all_Xmat_shift,modFitParams.add_usfac,modFitParams.flen);
    end

    %compute model-rates for the second version of the retinal stim
    for cc = 1:length(targs)
        if has_stim_mod(cc)
            [~,~,ep_rates2(:,cc)] = NMMmodel_eval(stim_mod(cc),[],all_Xmat_shift);
        end
    end
    ep_rates2 = bsxfun(@minus,ep_rates2,nanmean(ep_rates2(uset2(:),:))); %subtract off avg rates

    both_uset = uset1 & uset2; %set of indices where both simulated rates are usable
    ep_rates(~both_uset(:),:) = nan;
    ep_rates2(~both_uset(:),:) = nan;
    tbt_ep_rates = reshape(ep_rates,target_uf,n_utrials,length(targs));
    tbt_ep_rates2 = reshape(ep_rates2,target_uf,n_utrials,length(targs));
    for pp = 1:length(poss_ubins) %loop over possible temporal down-sampling factors
        bin_dsfac = poss_ubins(pp);
        n_newbins = floor(target_uf/bin_dsfac);
        new_ep_rates = zeros(n_newbins,n_utrials,length(targs),bin_dsfac);
        new_ep_rates2 = zeros(n_newbins,n_utrials,length(targs),bin_dsfac);
        for ii = 1:bin_dsfac
            new_ep_rates(:,:,:,ii) = tbt_ep_rates(ii:bin_dsfac:(ii+bin_dsfac*(n_newbins-1)),:,:,:);
            new_ep_rates2(:,:,:,ii) = tbt_ep_rates2(ii:bin_dsfac:(ii+bin_dsfac*(n_newbins-1)),:,:,:);
        end
        new_ep_rates = squeeze(sum(new_ep_rates,4)); %now sum over the coarser time bins. This makes the whole bin a NAN if any component bins are NAN
        new_ep_rates2 = squeeze(sum(new_ep_rates2,4)); 
        
        new_ep_rates = reshape(new_ep_rates,[],length(targs));
        new_ep_rates2 = reshape(new_ep_rates2,[],length(targs));
        new_ep_rates = bsxfun(@minus,new_ep_rates,nanmean(new_ep_rates)); %make sure mean rates are still 0
        new_ep_rates2 = bsxfun(@minus,new_ep_rates2,nanmean(new_ep_rates2));
        
        sim_psth_vars(sd,pp,:) = nanmean(new_ep_rates.*new_ep_rates2); %this gives the PSTH variance
        sim_tot_vars(sd,pp,:) = 0.5*nanvar(new_ep_rates) + 0.5*nanvar(new_ep_rates2); %total rate variance
        sim_mean_rates(sd,pp,:) = cur_avg_rates*bin_dsfac; %mean spike prob (per bin)
        
        %get PSTH and total covariance matrices
        if do_xcorr
            sim_psth_covars(sd,pp,:,:) = nanmean(bsxfun(@times,new_ep_rates,reshape(new_ep_rates2,[],1,length(targs))));
            sim_tot_covars(sd,pp,:,:) = 0.5*nanmean(bsxfun(@times,new_ep_rates,reshape(new_ep_rates,[],1,length(targs)))) + ...
                0.5*nanmean(bsxfun(@times,new_ep_rates2,reshape(new_ep_rates2,[],1,length(targs))));
        end
    end
end

%% store results
for cc = 1:length(targs)
    if has_stim_mod(cc)
        sim_data(cc).psth_vars = squeeze(sim_psth_vars(:,:,cc));
        sim_data(cc).tot_vars = squeeze(sim_tot_vars(:,:,cc));
        sim_data(cc).mean_rates = squeeze(sim_mean_rates(:,:,cc));
        sim_data(cc).across_trial_vars = sim_data(cc).tot_vars - sim_data(cc).psth_vars;
        sim_data(cc).FF_bias = sim_data(cc).across_trial_vars./sim_data(cc).mean_rates;
        sim_data(cc).alphas = sim_data(cc).across_trial_vars./sim_data(cc).tot_vars;
    end
end

if do_xcorr
    for cc = 1:n_cell_pairs
        sim_pairs(cc).ids = Cpairs(cc,:); %store index values of the neurons in this pair
        sim_pairs(cc).tot_xcovar = squeeze(sim_tot_covars(:,:,Cpairs(cc,1),Cpairs(cc,2))); %raw spk count covariances (this matrix is already symmetrized so we only need to grab one value)
        sim_pairs(cc).psth_xcovar = squeeze(0.5*sim_psth_covars(:,:,Cpairs(cc,1),Cpairs(cc,2)) + 0.5*sim_psth_covars(:,:,Cpairs(cc,2),Cpairs(cc,1))); %avg over two elements of psth covariance matrix to get best estimate
        cell1_noisevars = sim_data(Cpairs(cc,1)).across_trial_vars + sim_data(Cpairs(cc,1)).mean_rates;
        cell2_noisevars = sim_data(Cpairs(cc,2)).across_trial_vars + sim_data(Cpairs(cc,2)).mean_rates;
        if ~isempty(cell1_noisevars) & ~isempty(cell2_noisevars)
            sim_pairs(cc).covar_norm = sqrt(cell1_noisevars.*cell2_noisevars);
        else
            sim_pairs(cc).covar_norm = nan(length(poss_SDs),length(poss_ubins));
        end
    end
end
%%
cd(anal_dir);

sname = [sname sprintf('_ori%d',bar_ori)];
if rec_number > 1
    sname = strcat(sname,sprintf('_r%d',rec_number));
end
if do_xcorr
    save(sname,'targs','sim_data','sim_pairs','sim_params');
else
    save(sname,'targs','sim_data','sim_params');
end