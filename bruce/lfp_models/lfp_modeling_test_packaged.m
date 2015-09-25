clear all
close all
%
global Expt_name bar_ori monk_name rec_type rec_number

Expt_name = 'M012';
monk_name = 'jbe';
bar_ori = 0; %bar orientation to use (only for UA recs)
rec_number = 1;

fit_unCor = false; %also fit models without eye corrections?
use_MUA = false; %use MUA in model-fitting
fit_rect = false; %split quad linear filter into two rectified

save_name = 'corrected_models_comp_FIN2';
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
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
fprintf('Loading %s\n',data_name);
load(data_name);

%%
Expt_num = str2num(Expt_name(2:end));

cd(data_dir);
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
save_dir = ['~/Analysis/bruce/' Expt_name '/models'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
et_mod_data_name = 'full_eyetrack_initmods_FIN2_Rinit';
et_anal_name = 'full_eyetrack_FIN2_Rinit';

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

%% general parameters
has_repeats = any(~isnan(params.rpt_seeds));
xv_frac = 0.2; %fraction of trials to use for XVAL

flen = 15; %time lags for ST filters
spatial_usfac = et_params.spatial_usfac;%spatial up-sampling factor
add_usfac = et_params.add_usfac; %additional (tent-basis) spatial up-sampling
base_usfac = spatial_usfac/add_usfac; %number of pixels used in model fits

full_nPix = params.full_nPix; %total number of bars to keep track of in stimulus
use_nPix = et_params.use_nPix; %number of (bars) used in models
% use_nPix = full_nPix;
dt = params.dt;

use_nPix_us = use_nPix*spatial_usfac; %number of pixels in stimulus filters
klen_us = use_nPix_us*flen; %number of parameters in stim-filters
sp_dx = et_params.sp_dx; %pixel size (deg)
mod_dx = sp_dx*add_usfac; %model pixel size

use_LOOXV = 1; %[0 is no LOO; 1 is SUs only; 2 is SU + MU]
all_stim_mat = decompressTernNoise(stimComp);

%RF center position on screen
ov_RF_pos = Expts{expt_data.used_blocks(1)}.Stimvals.rf(1:2)/params.scale_fac;
fix_point = [Expts{expt_data.used_blocks(1)}.Stimvals.fx Expts{expt_data.used_blocks(1)}.Stimvals.fy];
ov_RF_pos = ov_RF_pos - fix_point; %account for non-zero fixation point location


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
if has_repeats%don't use repeat trials for training or cross-validation here.
    rpt_trials = find(ismember([trial_data(:).se],params.rpt_seeds));
    rpt_inds = find(ismember(all_trialvec(used_inds),rpt_trials));
else
    rpt_trials = [];
    rpt_inds = [];
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

sac_shift = et_params.sac_shift; %forward projection of saccade start times

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

%% LOAD LFP DATA
% cd(data_dir)
% if strcmp(rec_type,'LP') %if its a laminar probe rec
%     Fs = 1000; %LFP sampling freq
%     dsf = 5; %integer down-sampling factor
%     Fsd = Fs/dsf; %resulting Fs
%     niqf = Fs/2; %nyquist freq
%     new_niqf = Fsd/2; %nyquist after downsampling
%     [bb,aa] = butter(2,0.8*new_niqf/niqf); %anti-aliasing filter
%     use_lfps = 1:1:params.n_probes;
% elseif strcmp(rec_type,'UA') %if its a Utah array rec
%     Fs = 400.0032;
%     dsf = 2;
%     Fsd = Fs/dsf;
%     niqf = Fs/2;
%     new_niqf = Fsd/2;
%     [bb,aa] = butter(4,0.8*new_niqf/niqf);
%     use_lfps = 1:1:params.n_probes;
% end
% 
% thresh_match_diff = .01; %maximum temporal difference (in sec) for finding trial matches in LFP data
% 
% %wavelet parameters
% nwfreqs = 10; %number of frequency channels
% min_freq = 2; max_freq = 30; %freq range
% min_scale = 1/max_freq*Fsd;
% max_scale = 1/min_freq*Fsd;
% wavetype = 'cmor1-1'; %mother wavelet
% scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
% wfreqs = scal2frq(scales,wavetype,1/Fsd);
% 
% [full_lfps,full_lfp_taxis,full_cwts] = deal([]);
% cur_toffset = 0; %intialize time-offset tracker
% ublock_set = 1:n_blocks;
% for ee = ublock_set
%     
%     fprintf('Loading LFPs, Expt %d of %d\n',ee,n_blocks);
%     cur_trials = find([trial_data(:).block_nums] == ee);
%     cur_trial_start_times = [trial_data(cur_trials).start_times];
%     if strcmp(rec_type,'LP')
%         fname = sprintf('%s%sA.%d.lfp.mat',monk_name,Expt_name,expt_data.used_blocks(ee));
%         load(fname);
%         
%         tlens = arrayfun(@(X) length(X.ftime),LFP.Trials); %number of LFP chunks per trial
%         bad_trials = find(tlens == 0); %get rid of trials with no recorded LFP
%         LFP.Trials(bad_trials) = [];
%         
%         %get trial timestamps
%         lfp_rec_starts = [LFP.Trials(:).ftime]/1e4 + cur_toffset; %set of timestamps for start of LFP data (often before start of trial)
%         lfp_trial_starts = [LFP.Trials(:).Start]/1e4 + cur_toffset; %set of timestamps for start of trial
%         lfp_trial_ends = [LFP.Trials(:).End]/1e4 + cur_toffset; %timestamps of end of LFP data
%         
%         [expt_lfp_t_axis,expt_lfps,expt_cwt] = deal([]);
%         for tt = 1:length(cur_trials) %loop over used trials
%             [match_diff,matching_LFP_trial] = min(abs(cur_trial_start_times(tt) - lfp_trial_starts));%find matching trial in LFP data
%             if match_diff > thresh_match_diff
%                 fprintf('Match error of %.4f detected\n',match_diff);
%             end
%             cur_LFP = double(LFP.Trials(matching_LFP_trial).LFP(:,use_lfps));
%             %build time-axis for LFP trial data
%             cur_npts = size(cur_LFP,1);
%             cur_t_end = lfp_rec_starts(matching_LFP_trial)+(cur_npts-1)/Fs;
%             cur_t_axis = (lfp_rec_starts(matching_LFP_trial):1/Fs:cur_t_end);
%             if dsf > 1 %anti-alias filter and down-sample if needed
%                 cur_LFP = filtfilt(bb,aa,cur_LFP);
%                 cur_LFP = downsample(cur_LFP,dsf);
%                 cur_t_axis = downsample(cur_t_axis,dsf);
%             end
%             %compute CWT on full-trial LFP data
%             cur_cwt = nan(length(cur_t_axis),length(wfreqs),length(use_lfps));
%             for cc = 1:length(use_lfps)
%                 cur_cwt(:,:,cc) = cwt(cur_LFP(:,cc),scales,'cmor1-1')';
%             end
%             if ~isempty(expt_lfp_t_axis)
%                 cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
%                 if cur_sp > 1 %if there is overlap with existing LFP data, just take the non-overlapping portion
%                     cur_LFP = cur_LFP(cur_sp:end,:);
%                     cur_t_axis = cur_t_axis(cur_sp:end);
%                     cur_cwt = cur_cwt(cur_sp:end,:,:);
%                 end
%             end
%             expt_lfp_t_axis = cat(1,expt_lfp_t_axis, cur_t_axis(:));
%             expt_lfps = cat(1,expt_lfps,cur_LFP);
%             expt_cwt = cat(1,expt_cwt,cur_cwt);
%         end
%     else
%         lfp_fname = sprintf('Expt%d_LFP.mat',cur_block_set(ee));
%         load(lfp_fname);
%         
%         cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
%         if dsf > 1
%             cur_lfps = filtfilt(bb,aa,cur_lfps);
%             expt_lfps = downsample(cur_lfps,dsf);
%             expt_lfp_t_axis = downsample(lfp_t_ax',dsf);
%         end
%         
%         %apply CWT to trial-concatenated data. Need to be careful about edge artifacts
%         expt_cwt = nan(size(expt_lfps,1),length(wfreqs),length(use_lfps));
%         for cc = 1:length(use_lfps)
%             expt_cwt(:,:,cc) = cwt(expt_lfps(:,cc),scales,'cmor1-1')';
%         end
%     end
%     
%     cur_uset = find(all_blockvec == ee); %data during this block
%     if ~isempty(cur_uset) %keep only LFP data within the range of the used data
%         uinds = find(expt_lfp_t_axis >= all_t_axis(cur_uset(1)) & expt_lfp_t_axis <= all_t_axis(cur_uset(end)));
%         full_lfps = cat(1,full_lfps,expt_lfps(uinds,:));
%         full_cwts = cat(1,full_cwts,expt_cwt(uinds,:,:));
%         full_lfp_taxis = cat(1,full_lfp_taxis,expt_lfp_t_axis(uinds));
%     end
%     cur_toffset = time_data.trial_toffset(ee);
% end

%% interpolate LFP data onto modeling time axis
% interp_lfps = interp1(full_lfp_taxis,full_lfps,all_t_axis);
% interp_lfps = interp_lfps(used_inds,:);
% interp_lfps = nanzscore(interp_lfps); %zscore normalize LFP amplitudes
% 
% %interpolate real and imaginary components of the CWT separately
% interp_cwt_real = interp1(full_lfp_taxis,real(full_cwts),all_t_axis);
% interp_cwt_imag = interp1(full_lfp_taxis,imag(full_cwts),all_t_axis);
% interp_cwt_mag = sqrt(interp_cwt_real.^2 + interp_cwt_imag.^2);
% 
% %normalize CWT components by the SD of the amplitude of each CWT channel
% interp_cwt_real = bsxfun(@rdivide,interp_cwt_real,nanstd(interp_cwt_mag));
% interp_cwt_imag = bsxfun(@rdivide,interp_cwt_imag,nanstd(interp_cwt_mag));

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

%% model-fitting parameters
% poss_smoothreg_scalefacs = logspace(-2,2,10); %possible scale factors to apply to smoothness reg strength
poss_smoothreg_scalefacs = 1; %possible scale factors to apply to smoothness reg strength

base_lambda_d2XT = 100;
base_lambda_L1 = 5;
init_lambda_d2XT = 100;

%max number of quad filter to use for fitting SUA
max_Emods = 4; %max excitatory quad filters
max_Imods = 4; %max inhibitory quad filters

mod_stim_params(1) = NIM.create_stim_params([flen use_nPix_us/add_usfac],'stim_dt',dt);
mod_stim_params(2) = NIM.create_stim_params([n_blocks],'stim_dt',dt);
Xmat{2} = Xblock(used_inds,:);

% fit_lfp_chs = 1:2:params.n_probes; %use every other probe for fitting
% Xreal = interp_cwt_real(used_inds,:,fit_lfp_chs);
% Ximag = interp_cwt_imag(used_inds,:,fit_lfp_chs);
% XLFP = cat(3,Xreal,Ximag); %concatenate real and imag along 'channel' axis (dim 2)
% mod_stim_params(3) = NIM.create_stim_params([length(wfreqs) 2*length(fit_lfp_chs) 1],...
%     'boundary_conds',[0 0 0],'split_pts',[2 length(fit_lfp_chs) 0]);
% Xmat{3} = [reshape(Xreal,NT,[]) reshape(Ximag,NT,[])];
% clear interp_cwt* XLFP

%%
silent = 1;

% for cc = targs
cc = 36;
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
    
    %create X matrix
    if add_usfac > 1
        %if doing additional spatial up-sampling use tent basis functions
        Xmat{1} = tb_proc_stim(all_Xmat_shift(used_inds,use_kInds_up),add_usfac,flen);
    else
        Xmat{1} = all_Xmat_shift(used_inds,use_kInds_up);
    end
    
    %% compute null model
    
    nullMod = NIM(mod_stim_params,'lin',1,'Xtargets',2);
    nullMod = nullMod.fit_filters(cur_Robs,Xmat,cur_tr_inds,'silent',silent);
    nullMod_full = nullMod.fit_filters(cur_Robs,Xmat,cur_full_inds,'silent',silent);
    nullMod = nullMod.fit_spkNL(cur_Robs,Xmat,cur_tr_inds,'silent',silent);
    nullMod_full = nullMod_full.fit_spkNL(cur_Robs,Xmat,cur_full_inds,'silent',silent);
    nullMod_xvLL = nullMod.eval_model(cur_Robs,Xmat,cur_xv_inds);
    
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
    
    %start out with just a linear filter
    nEfilts = 3;
    nIfilts = 3;
    optim_params.optTol = 1e-4;
    optim_params.progTol = 1e-8;
    optim_params.silent = 0;

    stim_mod_signs = [1 1 ones(1,nEfilts) -1*ones(1,nIfilts)];
    stim_NL_types = [{'lin'} {'lin'} repmat({'quad'},1,nEfilts) repmat({'quad'},1,nIfilts)];
    stim_Xtargs = [2 1 ones(1,nEfilts+nIfilts)];
    gqm = NIM(mod_stim_params,stim_NL_types,stim_mod_signs,'Xtargets',stim_Xtargs,...
        'd2xt',[0 ones(1,sum(stim_Xtargs==1))]*cur_init_lambda_d2XT);
    gqm = gqm.fit_filters(cur_Robs,Xmat,cur_tr_inds,'optim_params',optim_params);
    
%     [~,~,mod_internals] = gqm.eval_model(cur_Robs,Xmat,cur_tr_inds);
%     vgint = var(mod_internals.gint(:,stim_Xtargs==1));
%     sgint = std(mod_internals.gint(:,stim_Xtargs==1));
%     new_d2xt = cur_base_lambda_d2XT./vgint;
%     new_l1 = cur_base_lambda_L1./sgint;
%     gqm = gqm.set_reg_params('sub_inds',find(stim_Xtargs==1),'d2xt',new_d2xt','l1',new_l1);
%     
%     gqm = gqm.fit_filters(cur_Robs,Xmat,cur_tr_inds,'optim_params',optim_params);
%     gqm = gqm.fit_spkNL(cur_Robs,Xmat,cur_tr_inds,'silent',silent);
    if isempty(cur_xv_inds)
        error('Need xval trials');
    end
    gqm_xvLL = gqm.eval_model(cur_Robs,Xmat,cur_xv_inds);
    gqm_xvLLimp = (gqm_xvLL - nullMod_xvLL)/log(2);
    
    %%
    cur_gam = 1.5;
    optim_params.MaxIter = 1000;
    filt_UF = 5;
    stim_NL_types = [{'lin'} repmat({'rectpow'},1,nEfilts*filt_UF) repmat({'rectpow'},1,nIfilts*filt_UF)];
    stim_mod_signs = [1 ones(1,nEfilts*filt_UF) -1*ones(1,nIfilts*filt_UF)];
    stim_Xtargs = [2 ones(1,nEfilts*filt_UF+nIfilts*filt_UF)];
    NLparams = cell(1,length(stim_mod_signs));
    [NLparams{3:end}] = deal(cur_gam);
    
    gqm2 = NIM(mod_stim_params,stim_NL_types,stim_mod_signs,'Xtargets',stim_Xtargs,'nlparams',NLparams,...
        'd2xt',[0 ones(1,sum(stim_Xtargs==1))]*cur_init_lambda_d2XT*4);
%     gqm2 = gqm2.fit_filters(cur_Robs,Xmat,cur_tr_inds,'fit_offsets',0,'optim_params',optim_params);
    gqm2 = gqm2.fit_filters(cur_Robs,Xmat,cur_tr_inds,'fit_offsets',1,'optim_params',optim_params);
%     gqm2 = gqm2.fit_spkNL(cur_Robs,Xmat,cur_tr_inds,'optim_params',optim_params);
    gqm2_xvLL = gqm2.eval_model(cur_Robs,Xmat,cur_xv_inds);

    %%
    cur_gam = 1.5;
    optim_params.MaxIter = 500;
    optim_params.optTol = 1e-4;
    optim_params.progTol = 1e-8;
    optim_params.silent = 0;
    nE = 4;
    nI = 4;
%     stim_NL_types = [{'lin'} repmat({'rectlin'},1,nE) repmat({'rectlin'},1,nI)];
    stim_NL_types = [{'lin'} repmat({'rectpow'},1,nE) repmat({'rectpow'},1,nI)];
    stim_mod_signs = [1 ones(1,nE) -1*ones(1,nI)];
    stim_Xtargs = [2 ones(1,nE+nI)];
    NLparams = cell(1,length(stim_mod_signs));
    [NLparams{2:end}] = deal(cur_gam);
    
    grm = NIM(mod_stim_params,stim_NL_types,stim_mod_signs,'Xtargets',stim_Xtargs,'nlparams',NLparams,...
        'd2xt',[0 ones(1,sum(stim_Xtargs==1))]*cur_init_lambda_d2XT*1);
    grm = grm.fit_filters(cur_Robs,Xmat,cur_tr_inds,'fit_offsets',1,'optim_params',optim_params);
    grm_xvLL = grm.eval_model(cur_Robs,Xmat,cur_xv_inds);
    
    %%
    grm2 = grm;
    %%
    grm2 = grm2.fit_NLparams(cur_Robs,Xmat,cur_tr_inds,'optim_params',optim_params);
    grm2 = grm2.fit_filters(cur_Robs,Xmat,cur_tr_inds,'fit_offsets',1,'optim_params',optim_params);

    %%
    optim_params.MaxIter = 1000;
    nE = 3;
    nI = 4;
    stim_NL_types = [{'lin'} {'lin'} repmat({'quad'},1,nE) repmat({'quad'},1,nI)];
    stim_mod_signs = [1 1 ones(1,nE) -1*ones(1,nI)];
    stim_Xtargs = [2 1 ones(1,nE+nI)];    
    gqm2 = NIM(mod_stim_params,stim_NL_types,stim_mod_signs,'Xtargets',stim_Xtargs,...
        'd2xt',[0 ones(1,sum(stim_Xtargs==1))]*cur_init_lambda_d2XT);
    gqm2 = gqm2.fit_filters(cur_Robs,Xmat,cur_tr_inds,'fit_offsets',0,'optim_params',optim_params);
    gqm2_xvLL = gqm2.eval_model(cur_Robs,Xmat,cur_xv_inds);
   
    %%
    [~,~,gqm_internals] = gqm2.eval_model(cur_Robs,Xmat);
    
    SL_Xmat{1} = gqm_internals.fgint(:,2:end);
    SL_Xmat{2} = Xmat{2};
    
    SL_stim_params(1) = NIM.create_stim_params(sum(stim_Xtargs == 1));
    SL_stim_params(2) = mod_stim_params(2);
    SL_nEfilts = 6;
    SL_nIfilts = 6;
    SL_NL_types = [{'lin'} {'lin'} repmat({'quad'},1,SL_nEfilts) repmat({'quad'},1,SL_nIfilts)];
    SL_mod_signs = [1 1 ones(1,SL_nEfilts) -1*ones(1,SL_nIfilts)];
    SL_Xtargs = [2 1 ones(1,SL_nEfilts+SL_nIfilts)];
    NLparams = cell(1,length(SL_mod_signs));
%     [NLparams{3:end}] = deal(cur_gam);
    
    slm = NIM(SL_stim_params,SL_NL_types,SL_mod_signs,'Xtargets',SL_Xtargs,'nlparams',NLparams,'l2',[0 500*ones(1,sum(SL_Xtargs == 1))*1]);
%     slm.subunits(2).filtK(find(stim_mod_signs(2:end)==1)) = 1;
%     slm.subunits(2).filtK(find(stim_mod_signs(2:end)==-1)) = -1;
    slm = slm.fit_filters(cur_Robs,SL_Xmat,cur_tr_inds,'fit_offsets',1,'optim_params',optim_params,'fit_subs',1:2);
    slm_xvLL = slm.eval_model(cur_Robs,SL_Xmat,cur_xv_inds);
    
    %%
    min_gam = 1;
    max_gam = 2.5;
    optim_params.silent = 1;
    for ii = 1:100
        init_gam(ii) = rand*(max_gam- min_gam) + min_gam;
        [NLparams{3:end}] = deal(init_gam(ii));
        
        gqm2 = NIM(mod_stim_params,stim_NL_types,stim_mod_signs,'Xtargets',stim_Xtargs,'nlparams',NLparams,...
            'd2xt',[0 ones(1,sum(stim_Xtargs==1))]*cur_init_lambda_d2XT);
        gqm2 = gqm2.fit_filters(cur_Robs,Xmat,cur_tr_inds,'fit_offsets',1,'optim_params',optim_params);
        all_mods(ii) = gqm2;
        all_LLs(ii) = gqm2.fit_props.LL;
    end
    
    %%
    lambda_d2t = 100;
    lambda_d2x = 100;
    lambda_l2 = 0;
    optim_params.optTol = 1e-4;
    optim_params.progTol = 1e-8;
    optim_params.silent = 0;

    [~,~,mod_internals] = gqm.eval_model(cur_Robs,Xmat);
    gqm_stim_subs = find(stim_Xtargs == 1);
    mod_outs = bsxfun(@times,mod_internals.fgint(:,gqm_stim_subs),stim_mod_signs(gqm_stim_subs));
    Efilts = find(stim_mod_signs(gqm_stim_subs) == 1); Efilts(1) = [];
    exc_out = sum(mod_outs(:,Efilts),2);
    inh_out = sum(mod_outs(:,stim_mod_signs(gqm_stim_subs) == -1),2);
    
    lfp_mod_signs = [1 1 1];
    lfp_NL_types = {'lin','lin','lin'};
    lfp_Xtargs = [2 3 4];
    lfp_add = NIM(mod_stim_params,lfp_NL_types,lfp_mod_signs,'Xtargets',lfp_Xtargs,...
        'd2t',[0 lambda_d2t 0],'d2x',[0 lambda_d2x 0],'l2',[0 lambda_l2 0]);
    
    lfp_add = lfp_add.fit_filters(cur_Robs,Xmat,cur_tr_inds,'optim_params',optim_params);
    lfp_add_xvLL = lfp_add.eval_model(cur_Robs,Xmat,cur_xv_inds);  
    
    %%
    mod_stim_params(4) = NIM.create_stim_params(length(gqm_stim_subs));
    Xmat{4} = mod_outs;
    
    cur_mod_signs = cat(2,lfp_mod_signs,ones(1,length(gqm_stim_subs)));
    cur_NL_types = cat(2,lfp_NL_types,repmat({'lin'},1,length(gqm_stim_subs)));
    cur_Xtargs = cat(2,lfp_Xtargs,ones(1,length(gqm_stim_subs))*3);
    
    gain_funs = cat(2,ones(NT,length(lfp_mod_signs)),Xmat{4});
    lfp_mult = NIM(mod_stim_params,cur_NL_types,cur_mod_signs,'Xtargets',cur_Xtargs,...
        'd2t',[0 lambda_d2t 0 lambda_d2t*ones(1,length(gqm_stim_subs))],...
        'd2x',[0 lambda_d2x 0 lambda_d2x*ones(1,length(gqm_stim_subs))],...
        'l2',[0 lambda_l2 0 lambda_l2*ones(1,length(gqm_stim_subs))]);
    lfp_mult.subunits(3).filtK(:) = 1;
    fit_subs = setdiff(1:length(cur_mod_signs),3);
    lfp_mult = lfp_mult.fit_filters(cur_Robs,Xmat,cur_tr_inds,'fit_subs',fit_subs,'gain_funs',gain_funs,'optim_params',optim_params);
    lfp_mult_xvLL = lfp_mult.eval_model(cur_Robs,Xmat,cur_xv_inds,'gain_funs',gain_funs);
    
    %%
    lfp_gain_filts = cell2mat(lfp_mult.get_filtKs(4:length(cur_mod_signs))');
    lfp_gain_filts = reshape(lfp_gain_filts,length(wfreqs),length(fit_lfp_chs),2,[]);
    lfp_gain_filt_amps = squeeze(sqrt(lfp_gain_filts(:,:,1,:).^2 + lfp_gain_filts(:,:,2,:).^2));
    lfp_gain_filt_phases = squeeze(atan2(lfp_gain_filts(:,:,1,:),lfp_gain_filts(:,:,2,:)));

    f1 = figure();
    for ii = 1:length(gqm_stim_subs)
        subplot(2,4,ii);
        pcolor(wfreqs,1:length(fit_lfp_chs),squeeze(lfp_gain_filt_amps(:,:,ii))');shading flat
    end
    f2 = figure();
    for ii = 1:length(gqm_stim_subs)
        subplot(2,4,ii);
        pcolor(wfreqs,1:length(fit_lfp_chs),squeeze(lfp_gain_filt_phases(:,:,ii))');shading flat
    end
    
    %%
    mod_stim_params(4) = NIM.create_stim_params(2);
    Xmat{4} = [exc_out inh_out];
    lambda_d2t = 500;
    lambda_d2x = 500;
    lambda_l2 = 0;
    cur_mod_signs = cat(2,lfp_mod_signs,ones(1,2));
    cur_NL_types = cat(2,lfp_NL_types,repmat({'lin'},1,2));
    cur_Xtargs = cat(2,lfp_Xtargs,ones(1,2)*3);
    
    gain_funs = cat(2,ones(NT,length(lfp_mod_signs)),Xmat{4});
    lfp_EI = NIM(mod_stim_params,cur_NL_types,cur_mod_signs,'Xtargets',cur_Xtargs,...
        'd2t',[0 lambda_d2t 0 lambda_d2t*ones(1,2)],...
        'd2x',[0 lambda_d2x 0 lambda_d2x*ones(1,2)],...
        'l2',[0 lambda_l2 0 lambda_l2*ones(1,2)]);
    lfp_EI.subunits(3).filtK(:) = 1;
    fit_subs = setdiff(1:length(cur_mod_signs),3);
    lfp_EI = lfp_EI.fit_filters(cur_Robs,Xmat,cur_tr_inds,'fit_subs',fit_subs,'gain_funs',gain_funs,'optim_params',optim_params);
    lfp_EI_xvLL = lfp_EI.eval_model(cur_Robs,Xmat,cur_xv_inds,'gain_funs',gain_funs);
    
    %%
    lfp_gain_filts = cell2mat(lfp_EI.get_filtKs(4:length(cur_mod_signs))');
    lfp_gain_filts = reshape(lfp_gain_filts,length(wfreqs),length(fit_lfp_chs),2,[]);
    lfp_gain_filt_amps = squeeze(sqrt(lfp_gain_filts(:,:,1,:).^2 + lfp_gain_filts(:,:,2,:).^2));
    lfp_gain_filt_phases = squeeze(atan2(lfp_gain_filts(:,:,1,:),lfp_gain_filts(:,:,2,:)));
    f1 = figure();
    for ii = 1:2
        subplot(2,1,ii);
        pcolor(wfreqs,1:length(fit_lfp_chs),squeeze(lfp_gain_filt_amps(:,:,ii))');shading flat
    end
    f2 = figure();
    for ii = 1:2
        subplot(2,1,ii);
        pcolor(wfreqs,1:length(fit_lfp_chs),squeeze(lfp_gain_filt_phases(:,:,ii))');shading flat
    end

    %%
    ch = 12;
    freq_ind = 5;
    cur_phase_signal = atan2(Ximag(:,freq_ind,ch),Xreal(:,freq_ind,ch));
    n_Gbins = 15;
    n_phase_bins = 20;
    Xtick = linspace(-pi,pi,n_phase_bins); 

    neural_signal = exc_out;
    
    exc_Ytick = my_prctile(neural_signal,linspace(0,99.5,n_Gbins)); %equipopulated binning
    TB_stim = [cur_phase_signal neural_signal];
    TB_used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
        TB_stim(:,2) >= exc_Ytick(1) & TB_stim(:,2) <= exc_Ytick(end));

    %initialize TBs
    TB = TentBasis2D(Xtick, exc_Ytick);
    [TB_Xmat{1},TB_counts] = TB.InputNL2D(TB_stim(TB_used_data,:));
    
    %phase is first dim, neural signal is second
    TB_stim_params = NIM.create_stim_params([n_phase_bins n_Gbins 1],'boundary_conds',[-1 Inf 0]);
    TB_mod = NIM(TB_stim_params,'lin',1,'d2t',10,'d2x',10,'d2xt',0);
    TB_mod  = TB_mod.fit_filters(cur_Robs(TB_used_data),TB_Xmat,'optim_params',optim_params);
     
    TB_K = reshape(TB_mod.subunits(1).filtK,n_phase_bins,n_Gbins);
    bin_areas = TB.GetBinAreas();
    TB_dist = TB_counts./bin_areas;
    exc_TB_rate = log(1 + exp(TB_K + TB_mod.spkNL.theta));

    
    neural_signal = inh_out;
    inh_Ytick = my_prctile(neural_signal,linspace(0,99.5,n_Gbins)); %equipopulated binning
    TB_stim = [cur_phase_signal neural_signal];
    TB_used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
        TB_stim(:,2) >= inh_Ytick(1) & TB_stim(:,2) <= inh_Ytick(end));

    %initialize TBs
    TB = TentBasis2D(Xtick, inh_Ytick);
    [TB_Xmat{1},TB_counts] = TB.InputNL2D(TB_stim(TB_used_data,:));
    
    %phase is first dim, neural signal is second
    TB_stim_params = NIM.create_stim_params([n_phase_bins n_Gbins 1],'boundary_conds',[-1 Inf 0]);
    TB_mod = NIM(TB_stim_params,'lin',1,'d2t',10,'d2x',10,'d2xt',0);
    TB_mod  = TB_mod.fit_filters(cur_Robs(TB_used_data),TB_Xmat,'optim_params',optim_params);
     
    TB_K = reshape(TB_mod.subunits(1).filtK,n_phase_bins,n_Gbins);
    bin_areas = TB.GetBinAreas();
    TB_dist = TB_counts./bin_areas;
    inh_TB_rate = log(1 + exp(TB_K + TB_mod.spkNL.theta));

    figure
    subplot(2,2,1)
    pcolor(Xtick,exc_Ytick,exc_TB_rate');shading flat;
    subplot(2,2,2)
    pcolor(Xtick,exc_Ytick,bsxfun(@rdivide,exc_TB_rate,mean(exc_TB_rate,1))');shading flat;
    subplot(2,2,3)
    pcolor(Xtick,inh_Ytick,inh_TB_rate');shading flat;set(gca,'yscale','log');
    subplot(2,2,4)
    pcolor(Xtick,inh_Ytick,bsxfun(@rdivide,inh_TB_rate,mean(inh_TB_rate,1))');shading flat;set(gca,'yscale','log');
 
    %%
    n_Ebins = 20;
    n_Ibins = 15;

    exc_tick = my_prctile(exc_out,linspace(0,99.5,n_Ebins)); %equipopulated binning
    inh_tick = my_prctile(inh_out,linspace(0,99.5,n_Ibins)); %equipopulated binning
    TB_stim = [exc_out inh_out];
    TB_used_data = find(TB_stim(:,1) >= exc_tick(1) & TB_stim(:,1) <= exc_tick(end) & ...
        TB_stim(:,2) >= inh_tick(1) & TB_stim(:,2) <= inh_tick(end));

    %initialize TBs
    TB = TentBasis2D(exc_tick, inh_tick);
    [TB_Xmat{1},TB_counts] = TB.InputNL2D(TB_stim(TB_used_data,:));
    
    %phase is first dim, neural signal is second
    TB_stim_params = NIM.create_stim_params([n_Ebins n_Ibins 1],'boundary_conds',[Inf Inf 0]);
    TB_mod = NIM(TB_stim_params,'lin',1,'d2t',10,'d2x',10,'d2xt',0);
    TB_mod  = TB_mod.fit_filters(cur_Robs(TB_used_data),TB_Xmat,'optim_params',optim_params);
     
    TB_K = reshape(TB_mod.subunits(1).filtK,n_Ebins,n_Ibins);
    bin_areas = TB.GetBinAreas();
    TB_dist = TB_counts./bin_areas;
    EI_TB_rate = log(1 + exp(TB_K + TB_mod.spkNL.theta));

    
end

