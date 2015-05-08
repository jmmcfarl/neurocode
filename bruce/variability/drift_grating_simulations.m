% clear all
% close all

global Expt_name bar_ori monk_name rec_type

% Expt_name = 'M012';
% monk_name = 'jbe';
% bar_ori = 0; %bar orientation to use (only for UA recs)
% rec_number = 1;

use_MUA = false;
fit_unCor = false;
use_hres_ET = true; EP_params.use_hres_ET = use_hres_ET; %use high-res eye-tracking?

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

save_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
save_name = 'grating_sim';
save_name = strcat(save_name,sprintf('_ori%d',bar_ori));
if rec_number > 1
    save_name = strcat(save_name,sprintf('_r%d',rec_number));
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

load([model_dir '/' mod_name]);

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

%% Recon retinal stim
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
sac_buff = round(0.05/dt);
sac_delay = round(0.03/dt);
blink_buff = round(0.1/dt);

% in_sac_inds = zeros(NT,1);
% nblink_start_inds = saccade_start_inds(~used_is_blink);
% for ii = 1:(sac_buff+1)
%     cur_inds = nblink_start_inds + sac_delay + (ii-1);
%     uu = find(cur_inds <= NT);
%     uu(all_trialvec(used_inds(cur_inds(uu))) ~= all_trialvec(used_inds(nblink_start_inds(uu)))) = [];
%     in_sac_inds(cur_inds(uu)) = 1;
% end
% in_sac_inds = logical(in_sac_inds);
% 
% blink_start_inds = saccade_start_inds(used_is_blink);
% in_blink_inds = zeros(NT,1);
% for ii = 1:(blink_buff+1)
%     cur_inds = blink_start_inds + sac_delay + (ii-1);
%     uu = find(cur_inds <= NT);
%     uu(all_trialvec(used_inds(cur_inds(uu))) ~= all_trialvec(used_inds(blink_start_inds(uu)))) = [];
%     in_blink_inds(cur_inds(uu)) = 1;
% end
% in_blink_inds = logical(in_blink_inds);

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

%% absorb block filter into spkNL offset parameter
for cc = targs
    if ~isempty(ModData(cc).bestGQM)
        cur_mod = ModData(cc).bestGQM;
        %absorb block-by-block offsets into overall spkNL offset param
        cur_block_filt = cur_mod.mods(1).filtK;
        cur_used_blocks = ModData(cc).unit_data.used_blocks;
        poss_used_blocks = ModData(cc).unit_data.poss_used_blocks;
        cur_used_blocks = find(ismember(cur_used_blocks,poss_used_blocks));
        cur_mod.spk_NL_params(1) = cur_mod.spk_NL_params(1) + mean(cur_block_filt(cur_used_blocks));
        cur_mod.mods(1) = [];
        GQM_mod{cc} = cur_mod;
    else
        GQM_mod{cc} = [];
    end
end

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
base_EP_est = post_mean_EP(full_uinds) - nanmedian(post_mean_EP(full_uinds));
EP_tbt = reshape(base_EP_est,target_uf,n_utrials);

%TBT mats for sac and blink indicators
inblink_tbt = reshape(in_blink_inds(full_uinds),target_uf,n_utrials);
insac_tbt = reshape(in_sac_inds(full_uinds),target_uf,n_utrials);


%% compute stats of model output on the stim ensemble it was trained on

flen = modFitParams.flen;
nPix_us = modFitParams.use_nPix_us;
stim_params = NMMcreate_stim_params([flen nPix_us]);
mod_usfac = modFitParams.spatial_usfac;
nPix = nPix_us/mod_usfac;

test_NT = 1e4;
test_dd = mode(expt_data.expt_dds);

test_stim = randi(2,test_NT,nPix);
test_stim(test_stim == 2) = -1;
is_zero = rand(test_NT,nPix) > test_dd/100;
test_stim(is_zero) = 0;

if mod_usfac > 1
    test_stim_up = zeros(size(test_stim,1),nPix_us);
    for ii = 1:size(test_stim,2)
        for jj = 1:mod_usfac
            test_stim_up(:,mod_usfac*(ii-1)+jj) = test_stim(:,ii);
        end
    end
elseif mod_usfac == 1
    test_stim_up = all_stim_mat;
end

test_stim_SD = nanstd(test_stim_up(:));

Xmat = create_time_embedding(test_stim_up,stim_params);

[test_gSD,test_gMean] = deal(cell(length(targs),1));
test_meanrate = nan(length(targs),1);
for ii = 1:length(targs)
    if ~isempty(GQM_mod{targs(ii)})
    [~,~,test_prate,~,gint] = NMMeval_model(GQM_mod{targs(ii)},[],Xmat);
    test_gSD{ii} = std(gint);
    test_gMean{ii} = mean(gint);
    test_meanrate(ii) = mean(test_prate);
    end
end

%%
poss_ubins = [1 2 5 10 20 40 80 100];
poss_grate_sf = [1 2 4];
poss_grate_tf = [2 4 8];

clear epscale*
target_EP_SD = 0.1;
EP_scale = target_EP_SD/robust_std_dev(EP_tbt(:));
full_EP_orth = EP_scale*EP_tbt(:);

xax = (1:nPix_us)*sp_dx;
nf = 375;
tax  = (1:nf)*params.dt;
throwout_win = 0.15/params.dt; %exclude this amount of data from the beginning of each trial to handle edge effects

for sf = 1:length(poss_grate_sf)
    for tf = 1:length(poss_grate_tf)
        fprintf('Simulating grating repeats for SF %d/%d and TF %d/%d\n',sf,length(poss_grate_sf),tf,length(poss_grate_tf));
        
        gr_sf = poss_grate_sf(sf);
        gr_tf = poss_grate_tf(tf);
        
        [XX,TT] = meshgrid(xax,tax);
        base_dg = sin(2*pi*(gr_sf*XX + gr_tf*TT));
        base_dg_rev = sin(2*pi*(gr_sf*XX - gr_tf*TT));
        
        
        [XX,TT,RR] = meshgrid(xax,tax,1:n_utrials); 
        flip_trials = rand(n_utrials,1) > 0.5; %randomly set half of trials to have opposite grating direction
        TT(:,:,flip_trials) = -TT(:,:,flip_trials);
        TT = reshape(permute(TT,[1 3 2]),[],length(xax));
        XX = reshape(permute(XX,[1 3 2]),[],length(xax));
        
        drift_grating = sin(2*pi*(gr_sf*bsxfun(@plus,XX,full_EP_orth) + gr_tf*TT));
        
        %% compute firing rate output of each model neuron
        
        Xmat = create_time_embedding(drift_grating,stim_params);
        base_Xmat = create_time_embedding(base_dg,stim_params);
        base_Xmat_rev = create_time_embedding(base_dg_rev,stim_params);
        
        base_prates = nan(nf,length(targs));
        base_prates_R = nan(nf,length(targs));
        rpt_prates = nan(nf,n_utrials,length(targs));
        for cc = 1:length(targs)
            if ~isempty(GQM_mod{targs(cc)})
            [~,~,cur_prate,~,gint] = NMMeval_model(GQM_mod{targs(cc)},[],Xmat);
            gint(inblink_tbt(:),:) = nan;
            gint(insac_tbt(:),:) = nan;
%             contrast_scale = mean(test_gSD{cc})/mean(nanstd(gint)); %adjust the contrast of the grating to match the mean SD of filter outputs 
            contrast_scale = 1;

            [~,~,pred_rate,~,gint] = NMMeval_model(GQM_mod{targs(cc)},[],Xmat*contrast_scale);
            rpt_prates(:,:,cc) = reshape(pred_rate,(nf),n_utrials);
            
            %these are the firing rates in response to the unperturbed
            %(left and rightward) gratings
            [~,~,base_prates(:,cc)] = NMMeval_model(GQM_mod{targs(cc)},[],base_Xmat*contrast_scale);
            [~,~,base_prates_R(:,cc)] = NMMeval_model(GQM_mod{targs(cc)},[],base_Xmat_rev*contrast_scale);
            end
        end
        
        base_prates(1:throwout_win,:) = [];
        base_prates_R(1:throwout_win,:) = [];
        
        rpt_prates = reshape(rpt_prates,[],length(targs));
        rpt_prates(inblink_tbt(:),:) = nan;
        rpt_prates_noSac = rpt_prates; %make a copy of the tbt firing rates with intrasac data nanned out
        rpt_prates_noSac(insac_tbt(:),:) = nan;
        rpt_prates = reshape(rpt_prates,[],n_utrials,length(targs));
        rpt_prates_noSac = reshape(rpt_prates_noSac,[],n_utrials,length(targs));
        
        rpt_prates(1:throwout_win,:,:) = [];
        rpt_prates_noSac(1:throwout_win,:,:) = [];
        
        pref_revdir = nanmean(base_prates_R) > nanmean(base_prates); %these units prefer the opposite motion direction
        dir_selectivity = abs(nanmean(base_prates_R) - nanmean(base_prates))./(nanmean(base_prates_R) + nanmean(base_prates));
        base_prates(:,pref_revdir) = base_prates_R(:,pref_revdir);
        
        amp_spectra = fft(base_prates)/size(base_prates,1);
        amp_spectra = abs(amp_spectra(1:size(base_prates,1)/2+1,:));
        f = 1/dt/2*linspace(0,1,size(base_prates,1)/2+1);
        first_harmonic = interp1(f,amp_spectra,gr_tf);
        second_harmonic = interp1(f,amp_spectra,2*gr_tf);
        F1F0 = 2*first_harmonic./amp_spectra(1,:);
        F2F1 = second_harmonic./first_harmonic;
 
        %% compute across-trial stats at different time resolutions
        
        [PSTH_vars,tot_vars,tot_means,FF_ests] = deal(nan(length(poss_ubins),length(targs)));
        [PSTH_vars_NS,tot_vars_NS,tot_means_NS,FF_ests_NS] = deal(nan(length(poss_ubins),length(targs)));
        for pp = 1:length(poss_ubins)
            bin_usfac = poss_ubins(pp);
            n_newbins = floor((nf-throwout_win)/bin_usfac);
            
            %'rebin' firnig rate data
            new_ep_rates = zeros(n_newbins,n_utrials,length(targs),bin_usfac);
            new_ep_rates_NS = zeros(n_newbins,n_utrials,length(targs),bin_usfac);
            for ii = 1:bin_usfac
                new_ep_rates(:,:,:,ii) = rpt_prates(ii:bin_usfac:(ii+bin_usfac*(n_newbins-1)),:,:);
                new_ep_rates_NS(:,:,:,ii) = rpt_prates_noSac(ii:bin_usfac:(ii+bin_usfac*(n_newbins-1)),:,:);
            end
            new_ep_rates = squeeze(nanmean(new_ep_rates,4))*bin_usfac;
            new_ep_rates_NS = squeeze(nanmean(new_ep_rates_NS,4))*bin_usfac;
            
            tot_vars(pp,:) = nanvar(reshape(new_ep_rates(:,~flip_trials,:),[],length(targs)));
            tot_means(pp,:) = nanmean(reshape(new_ep_rates(:,~flip_trials,:),[],length(targs)));
            tot_vars(pp,pref_revdir) = nanvar(reshape(new_ep_rates(:,flip_trials,pref_revdir),[],sum(pref_revdir)));
            tot_means(pp,pref_revdir) = nanmean(reshape(new_ep_rates(:,flip_trials,pref_revdir),[],sum(pref_revdir)));
            tot_vars_NS(pp,:) = nanvar(reshape(new_ep_rates_NS(:,~flip_trials,:),[],length(targs)));
            tot_means_NS(pp,:) = nanmean(reshape(new_ep_rates_NS(:,~flip_trials,:),[],length(targs)));
            tot_vars_NS(pp,pref_revdir) = nanvar(reshape(new_ep_rates_NS(:,flip_trials,pref_revdir),[],sum(pref_revdir)));
            tot_means_NS(pp,pref_revdir) = nanmean(reshape(new_ep_rates_NS(:,flip_trials,pref_revdir),[],sum(pref_revdir)));

            at_avgs = squeeze(nanmean(new_ep_rates(:,~flip_trials,:),2));
            at_vars = squeeze(nanvar(new_ep_rates(:,~flip_trials,:),[],2));
            at_avgs(:,pref_revdir) = squeeze(nanmean(new_ep_rates(:,flip_trials,pref_revdir),2));
            at_vars(:,pref_revdir) = squeeze(nanvar(new_ep_rates(:,flip_trials,pref_revdir),[],2));

            PSTH_vars(pp,:) = nanvar(at_avgs);
            FF_ests(pp,:) = nanmean((at_avgs + at_vars)./at_avgs);
            
            at_avgs_NS = squeeze(nanmean(new_ep_rates_NS(:,~flip_trials,:),2));
            at_vars_NS = squeeze(nanvar(new_ep_rates_NS(:,~flip_trials,:),[],2));
            at_avgs_NS(:,pref_revdir) = squeeze(nanmean(new_ep_rates_NS(:,flip_trials,pref_revdir),2));
            at_vars_NS(:,pref_revdir) = squeeze(nanvar(new_ep_rates_NS(:,flip_trials,pref_revdir),[],2));
           
            PSTH_vars_NS(pp,:) = nanvar(at_avgs_NS);
            FF_ests_NS(pp,:) = nanmean((at_avgs_NS + at_vars_NS)./at_avgs_NS);
        end
        
        %%
        for cc = 1:length(targs)
            grate_Cdata(cc).PSTH_vars(sf,tf,:) = PSTH_vars(:,cc);
            grate_Cdata(cc).PSTH_vars_NS(sf,tf,:) = PSTH_vars_NS(:,cc);
            
            grate_Cdata(cc).tot_vars(sf,tf,:) = tot_vars(:,cc);
            grate_Cdata(cc).tot_vars_NS(sf,tf,:) = tot_vars_NS(:,cc);
            
            grate_Cdata(cc).tot_means(sf,tf,:) = tot_means(:,cc);
            grate_Cdata(cc).tot_means_NS(sf,tf,:) = tot_means_NS(:,cc);
            
            grate_Cdata(cc).FF_ests(sf,tf,:) = FF_ests(:,cc);
            grate_Cdata(cc).FF_ests_NS(sf,tf,:) = FF_ests_NS(:,cc);
            
            grate_Cdata(cc).F1F0(sf,tf) = F1F0(cc);
            grate_Cdata(cc).F2F1(sf,tf) = F2F1(cc);
            grate_Cdata(cc).pref_revdir(sf,tf) = pref_revdir(cc);
            grate_Cdata(cc).dir_selectivity(sf,tf) = dir_selectivity(cc);
        end
        
    end
end

%%
cd(save_dir)
save(save_name,'grate_Cdata','poss_ubins','poss_grate_sf','poss_grate_tf');