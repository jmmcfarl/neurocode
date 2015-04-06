clear all
close all

global Expt_name bar_ori monk_name rec_type

Expt_name = 'M296';
monk_name = 'lem';
bar_ori = 45; %bar orientation to use (only for UA recs)

% [266-80 270-60 275-135 277-70 281-140 287-90 289-160 294-40 296-45 297-0/90]

fit_unCor = false; %use eye correction
use_MUA = false; %use MUA in model-fitting
fit_rect = false; %split quad linear filter into two rectified
use_hres_ET = true; %use high-res eye-tracking?

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
fprintf('Loading %s\n',data_name);
load(data_name);

ov_RF_pos = Expts{expt_data.used_blocks(1)}.Stimvals.rf(1:2)/params.scale_fac;

%%
Expt_num = str2num(Expt_name(2:end));

cd(data_dir);
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';

%if using coil info
if any(params.use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end
et_hres_anal_name = strcat(et_anal_name,'_hres');

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
et_hres_anal_name = [et_hres_anal_name sprintf('_ori%d',bar_ori)];

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

%% BIN SPIKES FOR MU AND SU
all_binned_mua = spikes_int82double(spike_data.binned_mua);
all_binned_sua = spikes_int82double(spike_data.binned_sua);
Clust_data = spike_data.Clust_data;
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

n_units = size(all_binned_mua,2) + size(all_binned_sua,2);

%% extract time series to keep track of trial number and block number
NT = length(used_inds);
fullNT = size(all_binned_mua,1);
n_trials = length(time_data.trial_flip_ids);
% n_blocks = length(expt_data.used_blocks);
n_blocks = length(time_data.block_flip_ids);

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

%% get binned spike data on finer time scale
orig_dt = params.dt;
base_dt = 0.005;

if n_trials ~= length(trial_data)
    error('Trial data mismatch');
end
up_nf = unique(expt_data.expt_nf)/base_dt*orig_dt;
tbt_binned_spikes = nan(up_nf,n_trials,length(SU_numbers));
tbt_t_axis = nan(up_nf,n_trials);
for ii = 1:n_trials
    cur_bin_edges = [trial_data(ii).start_times:base_dt:(trial_data(ii).start_times + base_dt*(up_nf))];
    cur_bin_cents = 0.5*cur_bin_edges(1:end-1) + 0.5*cur_bin_edges(2:end);
    for cc = 1:length(SU_numbers)
       cur_hist = histc(spike_data.SU_spk_times{cc},cur_bin_edges);
       tbt_binned_spikes(:,ii,cc) = cur_hist(1:end-1);
    end
    tbt_t_axis(:,ii) = cur_bin_cents;
end

orig_t_ind = round(interp1(time_data.t_axis,1:fullNT,tbt_t_axis(:)));
orig_t_ind(1) = 1;
up_used_inds = find(ismember(orig_t_ind,used_inds));

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

sac_shift = et_params.sac_shift; %forward projection of saccade start times
if use_hres_ET %if using high-res ET
    [post_mean_EP,post_std_EP] = construct_eye_position(best_fix_cor,best_fix_std,...
        drift_post_mean,drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    sp_dx = et_params.sp_dx;
    post_mean_EP = post_mean_EP*sp_dx;
    
    post_mean_EP_LOO = nan(length(loo_set),NT);
    for ss = 1:length(loo_set)
        [post_mean_EP_LOO(ss,:)] = construct_eye_position(best_fix_cor,best_fix_std,...
            drift_post_mean_LOO(ss,:),drift_post_std_LOO(ss,:),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    end
    post_mean_EP_LOO = post_mean_EP_LOO*sp_dx;
    
else %if using base resolution ET
    [post_mean_EP,post_std_EP] = construct_eye_position(it_fix_post_mean(end,:),it_fix_post_std(end,:),...
        drift_post_mean(end,:),drift_post_std(end,:),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    sp_dx = et_params.sp_dx;
    post_mean_EP = post_mean_EP*sp_dx;
    
    post_mean_EP_LOO = nan(length(loo_set),NT);
    for ss = 1:length(loo_set)
        [post_mean_EP_LOO(ss,:)] = construct_eye_position(squeeze(it_fix_post_mean_LOO(ss,end,:)),squeeze(it_fix_post_std_LOO(ss,end,:)),...
            squeeze(drift_post_mean_LOO(ss,end,:)),squeeze(drift_post_std_LOO(ss,end,:)),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    end
    post_mean_EP_LOO = post_mean_EP_LOO*sp_dx;
end

%%
all_trial_Se = [trial_data(:).se];
rpt_trials = find(ismember(all_trial_Se,params.rpt_seeds));
n_rpt_trials = length(rpt_trials);

rpt_taxis = (1:round(params.trial_dur)/base_dt)*base_dt-base_dt/2;
rpt_taxis(rpt_taxis < params.beg_buffer) = [];
rpt_taxis(params.trial_dur - rpt_taxis < params.end_buffer) = [];

bad_rpt_trials = find([trial_data(rpt_trials).rpt_frames] > 0); %get rid of any repeat trials where there were repeat frames
rpt_trials(bad_rpt_trials) = [];
fprintf('Using %d of %d repeat trials\n',length(rpt_trials),length(rpt_trials)+length(bad_rpt_trials));

n_rpts = length(rpt_trials);
all_rpt_inds = find(ismember(all_trialvec(used_inds),rpt_trials));
all_nonrpt_inds = find(~ismember(all_trialvec(used_inds),rpt_trials));

%% DEFINE IN-SAC AND IN-BUFF INDS
postsac_buff = round(0.03/base_dt);
rpt_saccade_set = find(ismember(saccade_trial_inds,rpt_trials));

in_sac_inds = zeros(up_nf*n_trials,1);
for ii = 1:sac_buff+1
    cur_inds = saccade_stop_inds + (ii-1);
    cur_inds(cur_inds > NT) = [];
    in_sac_inds(cur_inds) = 1;
end
for ii = 1:length(saccade_start_inds)
    cur_inds = saccade_start_inds(ii):saccade_stop_inds(ii);
    in_sac_inds(cur_inds) = 1;
end

in_sac_inds = logical(in_sac_inds);

blink_inds = find(used_is_blink);
in_blink_inds = zeros(NT,1);
for ii = 1:length(blink_inds)
    cur_inds = saccade_start_inds(blink_inds(ii)):saccade_stop_inds(blink_inds(ii));
    in_blink_inds(cur_inds) = 1;
end
in_blink_inds = logical(in_blink_inds);

%%
interp_post_mean_EP = interp1(time_data.t_axis(used_inds),post_mean_EP,tbt_t_axis(up_used_inds));
tbt_EP = nan(up_nf,n_trials);
tbt_EP(up_used_inds) = interp_post_mean_EP;


