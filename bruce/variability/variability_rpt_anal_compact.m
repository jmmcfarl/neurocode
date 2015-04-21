clear all
close all

addpath('~/other_code/fastBSpline/');

global Expt_name bar_ori monk_name rec_type

Expt_name = 'M275';
monk_name = 'lem';
bar_ori = 135; %bar orientation to use (only for UA recs)

% [266-80 270-60 275-135 277-70 281-140 287-90 289-160 294-40 296-45 297-0/90]

use_MUA = false; %use MUA in model-fitting
use_hres_ET = true; EP_params.use_hres_ET = use_hres_ET; %use high-res eye-tracking?
exclude_sacs = true; EP_params.exclude_sacs = exclude_sacs;
base_dt = 0.01; EP_params.base_dt = base_dt;

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

sname = 'rpt_variability_compact';

%%
Expt_num = str2num(Expt_name(2:end));

cd(data_dir);
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
model_dir = ['~/Analysis/bruce/' Expt_name '/models'];
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';
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
    post_mean_EP = post_mean_EP*sp_dx;
    
    if exist('drift_post_mean_LOO','var')
        post_mean_EP_LOO = nan(length(loo_set),NT);
        for ss = 1:length(loo_set)
            [post_mean_EP_LOO(ss,:)] = construct_eye_position(best_fix_cor,best_fix_std,...
                drift_post_mean_LOO(ss,:),drift_post_std_LOO(ss,:),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
        end
        post_mean_EP_LOO = post_mean_EP_LOO*sp_dx;
    else
        warning('No LOO variables detected');
    end
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

%% IDENTIFY REPEAT TRIALS
all_trial_Se = [trial_data(:).se];
n_rpt_seeds = length(params.rpt_seeds); EP_params.n_rpt_seeds = n_rpt_seeds;

all_rpt_trials = find(ismember(all_trial_Se,params.rpt_seeds));
bad_rpt_trials = find([trial_data(all_rpt_trials).rpt_frames] > 0); %get rid of any repeat trials where there were repeat frames
all_rpt_trials(bad_rpt_trials) = [];

nf = 400; %number of frames
used_Tinds = (params.beg_buffer/params.dt + 1):(nf - params.end_buffer/params.dt); %set of time points within each trial for analysis (excluding buffer windows)
used_nf = length(used_Tinds); %number of used time points per trial at this resolution

T = tabulate(all_trialvec(used_inds));
rpt_tdurs = T(all_rpt_trials,2);
bad_rpt_trials = find(rpt_tdurs ~= used_nf);
if ~isempty(bad_rpt_trials)
fprintf('Eliminating %d/%d incomplete repeats\n',length(bad_rpt_trials),length(all_rpt_trials));
all_rpt_trials(bad_rpt_trials) = [];
end

all_rpt_seqnum = nan(size(all_rpt_trials));
for ii = 1:n_rpt_seeds
    cur_trials = find(all_trial_Se(all_rpt_trials) == params.rpt_seeds(ii));
    all_rpt_seqnum(cur_trials) = ii;
end

rpt_taxis = (1:round(params.trial_dur)/base_dt)*base_dt-base_dt/2;
rpt_taxis(rpt_taxis < params.beg_buffer) = [];
rpt_taxis(params.trial_dur - rpt_taxis < params.end_buffer) = [];

tot_nrpts = length(all_rpt_seqnum);
fprintf('Using %d repeat trials, %d sequences\n',tot_nrpts,length(params.rpt_seeds));

rpt_trial_block = [trial_data(all_rpt_trials).block_nums];
%% get binned spike data on finer time scale

if n_trials ~= length(trial_data)
    error('Trial data mismatch');
end
up_nf = nf/base_dt*orig_dt;
tbt_binned_spikes = nan(up_nf,tot_nrpts,length(SU_numbers));
tbt_t_axis = nan(up_nf,tot_nrpts);
for ii = 1:tot_nrpts
    cur_bin_edges = [trial_data(all_rpt_trials(ii)).start_times:base_dt:(trial_data(all_rpt_trials(ii)).start_times + base_dt*(up_nf))];
    cur_bin_cents = 0.5*cur_bin_edges(1:end-1) + 0.5*cur_bin_edges(2:end);
    for cc = 1:length(SU_numbers)
        if ~isnan(Clust_data.SU_block_probes(cc,rpt_trial_block(ii)))
            cur_hist = histc(spike_data.SU_spk_times{cc},cur_bin_edges);
            tbt_binned_spikes(:,ii,cc) = cur_hist(1:end-1);
        end
    end
    tbt_t_axis(:,ii) = cur_bin_cents;
end

orig_t_ind = round(interp1(time_data.t_axis,1:fullNT,tbt_t_axis(:)));
orig_t_ind(1) = 1;
up_used_inds = find(ismember(orig_t_ind,used_inds));

tbt_binned_spikes = reshape(tbt_binned_spikes,[],length(SU_numbers));
for ii = 1:length(SU_numbers)
    tbt_binned_spikes(isnan(all_binned_sua(orig_t_ind,ii)),ii) = nan;
end
tbt_binned_spikes = reshape(tbt_binned_spikes,up_nf,tot_nrpts,length(SU_numbers));

used_Tinds = (params.beg_buffer/base_dt + 1):(up_nf - params.end_buffer/base_dt); %set of time points within each trial for analysis (excluding buffer windows)
used_up_nf = length(used_Tinds); %number of used time points per trial at this resolution

%% IDENTIFY TIMES WITHIN SACCADES AND BLINKS
sac_buff = round(0.06/params.dt); EP_params.sac_buff = sac_buff; %window of data to exclude during saccades
sac_delay = round(0.04/params.dt); EP_params.sac_delay = sac_delay; %shift exclusion window to account for neural delay
blink_buff = round(0.1/params.dt); EP_params.blink_buff = blink_buff; %window for excluding blinks

in_sac_inds = zeros(NT,1);
nblink_start_inds = saccade_start_inds(~used_is_blink);
for ii = 1:(sac_buff+1)
    cur_inds = nblink_start_inds + sac_delay + (ii-1);
    uu = find(cur_inds <= NT);
    uu(all_trialvec(used_inds(cur_inds(uu))) ~= all_trialvec(used_inds(nblink_start_inds(uu)))) = [];
    in_sac_inds(cur_inds(uu)) = 1;
end

blink_start_inds = saccade_start_inds(used_is_blink);
in_blink_inds = zeros(NT,1);
for ii = 1:(blink_buff+1)
    cur_inds = blink_start_inds + sac_delay + (ii-1);
    uu = find(cur_inds <= NT);
    uu(all_trialvec(used_inds(cur_inds(uu))) ~= all_trialvec(used_inds(blink_start_inds(uu)))) = [];
    in_blink_inds(cur_inds(uu)) = 1;
end

in_sac_inds = reshape(logical(ceil(interp1(all_t_axis(used_inds),in_sac_inds,tbt_t_axis(up_used_inds)))),used_up_nf,tot_nrpts);
in_blink_inds = reshape(logical(ceil(interp1(all_t_axis(used_inds),in_blink_inds,tbt_t_axis(up_used_inds)))),used_up_nf,tot_nrpts);

%% process trial-by-trial binned spike data (exclude blinks, sacs, subtract trial avgs)

tbt_BS_ms = reshape(tbt_binned_spikes(used_Tinds,:,:),[],length(SU_numbers));
tbt_BS_ms(in_blink_inds(:),:) = nan;

if exclude_sacs
    tbt_BS_ms(in_sac_inds(:),:) = nan;
end

tbt_BS_ms = reshape(tbt_BS_ms,used_up_nf,tot_nrpts,length(SU_numbers));

%subtract out trial-avg spike count
ov_avg_BS = nanmean(reshape(tbt_binned_spikes(used_Tinds,:,:),[],length(SU_numbers)));
trial_avg_BS = nanmean(tbt_binned_spikes(used_Tinds,:,:));
tbt_BS_ms = bsxfun(@minus,tbt_BS_ms,trial_avg_BS);
trial_avg_BS = squeeze(trial_avg_BS);

%now subtract out overall avg spike count
tbt_BS_ms = bsxfun(@minus,tbt_BS_ms,reshape(nanmean(reshape(tbt_BS_ms,[],length(SU_numbers))),[1 1 length(SU_numbers)]));

for cc = 1:length(targs)
    EP_data(cc).ov_avg_BS = ov_avg_BS(cc);
    EP_data(cc).trial_avg_BS = trial_avg_BS(:,cc);
end

%% COMPUTE BASIC STATS FOR REPEAT TRIALS
for rr = 1:n_rpt_seeds
    cur_trial_set = find(all_rpt_seqnum == rr);
    cur_nrpts = length(cur_trial_set);
    %% GET BASIC STATS OF TBT DATA
    
    psths = squeeze(nanmean(tbt_BS_ms(:,cur_trial_set,:),2));
    psth_var = nanvar(psths);
    tot_resp_var = nanvar(reshape(tbt_BS_ms(:,cur_trial_set,:),[],length(SU_numbers)));
    
    n_utrials = squeeze(mean(sum(~isnan(tbt_BS_ms(:,cur_trial_set,:)),2)));
    n_spikes = squeeze(nansum(reshape(tbt_binned_spikes(used_Tinds,:,:),[],length(SU_numbers))));
    
    avg_temp_var = squeeze(nanmean(nanvar(tbt_BS_ms(:,cur_trial_set,:)))); %avg (across trials) of across-time variance
    psth_var_cor = psth_var.*(n_utrials'./(n_utrials'-1)) - avg_temp_var'./n_utrials'; %sahani linden correction for PSTH sampling noise
    
    trial_avg_var = squeeze(nanvar(trial_avg_BS)); %variance of trial-avg rates
    
    for cc = 1:length(targs)
        EP_data(cc).psths(rr,:) = psths(:,cc);
        EP_data(cc).psth_var(rr) = psth_var(cc);
        EP_data(cc).psth_var_cor(rr) = psth_var_cor(cc);
        EP_data(cc).tot_var(rr) = tot_resp_var(cc);
        EP_data(cc).n_utrials(rr) = n_utrials(cc);
        EP_data(cc).n_spikes(rr) = n_spikes(cc);
        EP_data(cc).trial_avg_var(rr) = trial_avg_var(cc);
    end
end

%% ESTIMATE OVERALL DISTRIBUTION OF DELTA_X
emb_win = round(0.06/base_dt); EP_params.emb_win = emb_win; %look back this many time steps to parse EP trajectories
emb_shift = round(0.04/base_dt); EP_params.emb_shift = emb_shift;

n_eval_pts = 100; EP_params.n_eval_pts = n_eval_pts; %number of points to evaluate spline fit

interp_post_mean_EP = interp1(time_data.t_axis(used_inds),post_mean_EP,tbt_t_axis(up_used_inds)); %interpolate EP data

tbt_EP = reshape(interp_post_mean_EP,used_up_nf,tot_nrpts);

%initialize a time-embedded version of the EPs
sp = NMMcreate_stim_params(emb_win + emb_shift);
tbt_EP_emb = create_time_embedding(tbt_EP(:),sp);
tbt_EP_emb(in_blink_inds(:),:) = nan;
if exclude_sacs
    tbt_EP_emb(in_sac_inds(:),:) = nan;
end
tbt_EP_emb = reshape(tbt_EP_emb(:,(emb_shift+1):end),used_up_nf,tot_nrpts,[]);

% compute the distribution of delta_X
rand_T = [];
for rr = 1:n_rpt_seeds
    cur_trial_set = find(all_rpt_seqnum == rr);
    
    %% estimate quantiles of the distribution of pairwise EP similarities by this metric
    rset = randi(used_up_nf,100,1);
    for ii = 1:length(rset)
        cur_Dmat = abs(squareform(pdist(squeeze(tbt_EP_emb(rset(ii),cur_trial_set,:)))))/sqrt(emb_win);
        cur_Dmat(logical(eye(length(cur_trial_set)))) = nan;
        cur_Dmat = cur_Dmat(~isnan(cur_Dmat));
        rand_T = cat(1,rand_T,cur_Dmat);
    end
end

eval_xx = unique([0 prctile(rand_T,linspace(0,100,n_eval_pts))]); %x-axis for evaluating spline models
EP_params.eval_xx = eval_xx;

%%
xvfold = 10; EP_params.xvfold = xvfold; %cross-val fold for estimating optimal number of splines
poss_n_splines = [2 3 4 6 8 10 14]; EP_params.poss_N_splines = poss_n_splines; %range of possible values for number of splines.
n_boot_samps = 20; EP_params.n_boot_samps = n_boot_samps; %number of bootstrap samples for estimating spline uncertainty

for cc = 1:length(targs)
    fprintf('SU %d/%d\n',cc,length(targs));
    loo_ind = find(loo_set == targs(cc));
    if ~isempty(loo_ind)
        interp_post_mean_EP = interp1(time_data.t_axis(used_inds),post_mean_EP_LOO(loo_ind,:),tbt_t_axis(up_used_inds));
        tbt_EP = reshape(interp_post_mean_EP,used_up_nf,tot_nrpts);
        
        %initialize a time-embedded version of the EPs
        sp = NMMcreate_stim_params(emb_win+emb_shift);
        loo_tbt_EP_emb = create_time_embedding(tbt_EP(:),sp);
        loo_tbt_EP_emb(in_blink_inds(:),:) = nan;
        if exclude_sacs
            loo_tbt_EP_emb(in_sac_inds(:),:) = nan;
        end
        loo_tbt_EP_emb = reshape(loo_tbt_EP_emb(:,(emb_shift+1):end),used_up_nf,tot_nrpts,[]);
        
        all_D = [];
        all_base_D = [];
        all_X = [];
        for rr = 1:n_rpt_seeds %loop over unique repeat seeds
            cur_trial_set = find(all_rpt_seqnum == rr);
            cur_nrpts = length(cur_trial_set);
            [II,JJ] = meshgrid(1:cur_nrpts);
            uset = JJ > II; %only need to count each unique trial pair once
            n_unique_pairs = sum(uset(:));
            
            cur_D = nan(n_unique_pairs*used_up_nf,1);
            cur_base_D = nan(n_unique_pairs*used_up_nf,1);
            cur_X = nan(n_unique_pairs*used_up_nf,1);
            for tt = 1:used_up_nf
                cur_inds = (tt-1)*n_unique_pairs + (1:n_unique_pairs);
                Y1 = squeeze(tbt_BS_ms(tt,cur_trial_set,cc));
                
                cur_Dmat = abs(squareform(pdist(squeeze(loo_tbt_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                cur_Dmat(logical(eye(cur_nrpts))) = nan;
                cur_D(cur_inds) = cur_Dmat(uset);
                cur_Dmat = abs(squareform(pdist(squeeze(tbt_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                cur_Dmat(logical(eye(cur_nrpts))) = nan;
                cur_base_D(cur_inds) = cur_Dmat(uset);
                
                cur_Xmat = bsxfun(@times,Y1(II),Y1(JJ));
                cur_X(cur_inds) = cur_Xmat(uset);
            end
            cur_upts = find(~isnan(cur_D) & ~isnan(cur_X));
            all_D = cat(1,all_D,cur_D(cur_upts));
            all_base_D = cat(1,all_base_D,cur_base_D(cur_upts));
            all_X = cat(1,all_X,cur_X(cur_upts));
        end
        n_data_points = length(all_X);
        
        %assign each data point (trial pair) to one of the cross-validation
        %sets
        n_xvpts = round(n_data_points/xvfold);
        pt_xvsets = ceil((1:n_data_points)/n_xvpts);
        rperm = randperm(n_data_points);
        pt_xvsets = pt_xvsets(rperm);
        
        xv_err = nan(xvfold,length(poss_n_splines));
        for pp = 1:length(poss_n_splines)
            fprintf('Computing xval for spline set %d/%d\n',pp,length(poss_n_splines));
            n_splines = poss_n_splines(pp);
            spline_DS = prctile(rand_T,100/n_splines:100/n_splines:(100-100/n_splines));
            knot_pts = [0 0 0 0 spline_DS spline_DS(end) spline_DS(end) spline_DS(end)];
            
            %get cross-validated MSE for each value of n_knots
            for ii = 1:xvfold
                cur_tr_sets = setdiff(1:xvfold,ii);
                cur_tr_pts = ismember(pt_xvsets,cur_tr_sets);
                cur_xv_pts = ~ismember(pt_xvsets,cur_tr_sets);
                sp = fastBSpline.lsqspline(knot_pts,3,all_D(cur_tr_pts),all_X(cur_tr_pts));
                sp_pred = sp.evalAt(all_D(cur_xv_pts));
                xv_err(ii,pp) = sqrt(mean((sp_pred - all_X(cur_xv_pts)).^2));
            end
        end
        avg_xv_err = mean(xv_err);
        [~,best_xval] = min(avg_xv_err);
        best_n_knots = poss_n_splines(best_xval);
        
        spline_DS = prctile(rand_T,100/best_n_knots:100/best_n_knots:(100-100/best_n_knots));
        knot_pts = [0 0 0 0 spline_DS spline_DS(end) spline_DS(end) spline_DS(end)];
        boot_spline_pred = nan(n_boot_samps,length(eval_xx));
        for nn =1:n_boot_samps
            fprintf('Computing bootstrap sample %d/%d\n',nn,n_boot_samps);
            cur_upts = randi(n_data_points,n_data_points,1); %sampling from data points with replacement
            sp = fastBSpline.lsqspline(knot_pts,3,all_D(cur_upts),all_X(cur_upts));
            boot_spline_pred(nn,:) = sp.evalAt(eval_xx);
        end
        
        sp = fastBSpline.lsqspline(knot_pts,3,all_base_D,all_X);
        EP_data(cc).spline_pred_baseEP = sp.evalAt(eval_xx);
        
        sp = fastBSpline.lsqspline(knot_pts,3,all_D,all_X);
        EP_data(cc).spline_pred_looEP = sp.evalAt(eval_xx);
        
        EP_data(cc).pair_psth_var = mean(all_X);
        EP_data(cc).spline_boot_avg = mean(boot_spline_pred);
        EP_data(cc).spline_boot_sd = std(boot_spline_pred);
        EP_data(cc).avg_xv_err = avg_xv_err;
        EP_data(cc).n_knots = best_n_knots;
    end
end
%% PROCESS MODEL FITS
has_stim_mod = false(length(targs),1);
for cc = 1:length(targs)
    if ~isempty(ModData(targs(cc)).bestGQM)
        cur_mod = ModData(targs(cc)).bestGQM;
        cur_block_filt = cur_mod.mods(1).filtK;
        cur_mod.spk_NL_params(1) = cur_mod.spk_NL_params(1) + mean(cur_block_filt); %absorb block filter into offset parameter
        cur_mod.mods(1) = []; %eliminate block filter
        stim_mod(cc) = cur_mod;
        EP_data(cc).unit_data = ModData(targs(cc)).unit_data;
        EP_data(cc).tune_props = ModData(targs(cc)).tune_props;
        has_stim_mod(cc) = true;
    end
end

%% get stimulus xmat during repeat trials

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

%indices of repeat trials
full_rpt_inds = find(ismember(all_trialvec,all_rpt_trials));
used_rpt_inds = find(ismember(all_trialvec(used_inds),all_rpt_trials));

fin_shift_cor = round(post_mean_EP/modFitParams.sp_dx); %use overall EP estimate

%RECOMPUTE XMAT
best_shift_stimmat_up = all_stimmat_up;
for i=1:length(used_rpt_inds)
    best_shift_stimmat_up(used_inds(used_rpt_inds(i)),:) = shift_matrix_Nd(all_stimmat_up(used_inds(used_rpt_inds(i)),:),-fin_shift_cor(used_rpt_inds(i)),2);
end
all_Xmat_shift = create_time_embedding(best_shift_stimmat_up(full_rpt_inds,:),stim_params_full);
all_Xmat_shift = all_Xmat_shift(ismember(full_rpt_inds,used_inds),use_kInds_up);

%%

all_mod_emp_prates = nan(length(used_rpt_inds),length(targs));
for cc = 1:length(targs)
    if has_stim_mod(cc)
        [~,~,all_mod_emp_prates(:,cc)] = NMMmodel_eval(stim_mod(cc),[],all_Xmat_shift);
    end
end
all_mod_emp_prates(in_blink_inds(:),:) = nan;
if exclude_sacs
    all_mod_emp_prates(in_sac_inds(:),:) = nan;
end
all_mod_emp_prates = reshape(all_mod_emp_prates,used_nf,tot_nrpts,length(targs));

for rr = 1:n_rpt_seeds
    cur_trial_set = find(all_rpt_seqnum == rr);
    mod_psths = squeeze(nanmean(all_mod_emp_prates(:,cur_trial_set,:),2));
    mod_cond_vars = squeeze(nanvar(all_mod_emp_prates(:,cur_trial_set,:),[],2));
    mod_tot_vars = squeeze(nanvar(reshape(all_mod_emp_prates(:,cur_trial_set,:),[],length(targs))));
    
    for cc = 1:length(targs)
        if has_stim_mod(cc)
            EP_data(cc).mod_psths(rr,:) = mod_psths(:,cc);
            EP_data(cc).mod_cond_vars(rr,:) = mod_cond_vars(:,cc);
            EP_data(cc).mod_tot_vars(rr) = mod_tot_vars(cc);
            EP_data(cc).mod_psth_vars(rr) = var(mod_psths(:,cc));
            EP_data(cc).mod_ep_vars(rr) = mean(mod_cond_vars(:,cc));
            EP_data(cc).mod_alphas(rr) = EP_data(cc).mod_psth_vars(rr)/EP_data(cc).mod_tot_vars(rr);
        end
    end
end
%%
anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
if ~exist(anal_dir,'dir')
    mkdir(anal_dir)
end
cd(anal_dir);

sname = [sname sprintf('_ori%d',bar_ori)];

save(sname,'targs','EP_data','EP_params','use_MUA');