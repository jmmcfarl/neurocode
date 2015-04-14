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
exclude_sacs = true;
base_dt = 0.01;

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

%% get binned spike data on finer time scale

if n_trials ~= length(trial_data)
    error('Trial data mismatch');
end
up_nf = unique(expt_data.expt_nf)/base_dt*orig_dt;
tbt_binned_spikes = nan(up_nf,n_rpts,length(SU_numbers));
tbt_t_axis = nan(up_nf,n_rpts);
for ii = 1:n_rpts
    cur_bin_edges = [trial_data(rpt_trials(ii)).start_times:base_dt:(trial_data(rpt_trials(ii)).start_times + base_dt*(up_nf))];
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

tbt_binned_spikes = reshape(tbt_binned_spikes,[],length(SU_numbers));
for ii = 1:length(SU_numbers)
   tbt_binned_spikes(isnan(all_binned_sua(orig_t_ind,ii)),ii) = nan; 
end
tbt_binned_spikes = reshape(tbt_binned_spikes,up_nf,n_rpts,length(SU_numbers));

%% DEFINE IN-SAC AND IN-BUFF INDS
postsac_buff = round(0.03/base_dt); %window after saccade stop time to exclude data
rpt_saccade_set = find(ismember(saccade_trial_inds,rpt_trials));

rpt_sac_start_inds = round(interp1(tbt_t_axis(:),1:length(tbt_t_axis(:)),sac_start_times(used_saccade_set(rpt_saccade_set))));
rpt_sac_stop_inds = round(interp1(tbt_t_axis(:),1:length(tbt_t_axis(:)),sac_stop_times(used_saccade_set(rpt_saccade_set))));

in_sac_inds = zeros(up_nf*n_rpts,1);
in_blink_inds = zeros(up_nf*n_rpts,1);
for ii = 1:length(rpt_saccade_set)
    cur_ind_range = rpt_sac_start_inds(ii):(rpt_sac_stop_inds(ii) + postsac_buff);
    cur_ind_range(all_trialvec(orig_t_ind(cur_ind_range)) ~= all_trialvec(orig_t_ind(cur_ind_range(1)))) = [];
    
    if ~used_is_blink(rpt_saccade_set(ii)) 
    in_sac_inds(cur_ind_range) = 1;
    else
    in_blink_inds(cur_ind_range) = 1;
    end
end

in_sac_inds = logical(reshape(in_sac_inds,up_nf,n_rpts));
in_blink_inds = logical(reshape(in_blink_inds,up_nf,n_rpts));

%% process trial-by-trial binned spike data (exclude blinks, sacs, subtract trial avgs)
tbt_BS_ms = reshape(tbt_binned_spikes,[],length(SU_numbers));
tbt_BS_ms(in_blink_inds(:),:) = nan;

if exclude_sacs
    tbt_BS_ms(in_sac_inds(:),:) = nan;
end

tbt_BS_ms = reshape(tbt_BS_ms,up_nf,n_rpts,length(SU_numbers));
%subtract out trial-avg spike count
trial_avg_BS = nanmean(tbt_binned_spikes);
tbt_BS_ms = bsxfun(@minus,tbt_binned_spikes,nanmean(trial_avg_BS));
trial_avg_BS = squeeze(trial_avg_BS);
%%
psths = squeeze(nanmean(tbt_BS_ms,2));
psth_var = nanvar(psths);
tot_resp_var = nanvar(reshape(tbt_BS_ms,[],length(SU_numbers)));

n_utrials = squeeze(mean(sum(~isnan(tbt_BS_ms),2)));

avg_temp_var = squeeze(nanmean(nanvar(tbt_BS_ms))); %avg (across trials) of across-time variance
psth_var_cor = psth_var.*(n_utrials'./(n_utrials'-1)) - avg_temp_var'./n_utrials'; %sahani linden correction for PSTH sampling noise

trial_avg_var = squeeze(nanvar(trial_avg_BS)); %variance of trial-avg rates


%%
interp_post_mean_EP = interp1(time_data.t_axis(used_inds),post_mean_EP,tbt_t_axis(up_used_inds));
tbt_EP = nan(up_nf,n_rpts);
tbt_EP(up_used_inds) = interp_post_mean_EP;

tbt_EP(in_blink_inds) = nan;
tbt_EP(in_sac_inds) = nan;
uinds = (params.beg_buffer/base_dt + 1):(up_nf - params.end_buffer/base_dt);

%%
Xtick = rpt_taxis(1):0.01:rpt_taxis(end);
Ytick = -0.5:0.025:0.5;
TB = TentBasis2D(Xtick, Ytick);
[RR,TT] = meshgrid(1:n_rpts,rpt_taxis);
TB_stim = [TT(:) reshape(tbt_EP(uinds,:),[],1)];

%process data with TBs
[TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim);

%%
su_num = 1;

robs = reshape(tbt_binned_spikes(uinds,:,su_num),[],1);
uset = find(~isnan(robs));

temp_smooth = 1;
eye_smooth = 1;
TB_stim_params = NMMcreate_stim_params([length(Xtick) length(Ytick)]);
TB_reg_params = NMMcreate_reg_params('lambda_d2T',temp_smooth,'lambda_d2X',eye_smooth,'boundary_conds',[Inf Inf]);
TB_mod = NMMinitialize_model(TB_stim_params,1,{'lin'},TB_reg_params);
TB_mod = NMMfit_filters(TB_mod,robs,TB_Xmat,[],uset,0);
TB_filt = reshape(TB_mod.mods(1).filtK,[length(Xtick) length(Ytick)]);
TB_rate = log(1+exp(TB_filt + TB_mod.spk_NL_params(1)))/base_dt;

% %%
% xv_frac = 0.2;
% 
% poss_trials = unique(RR(uset));
% nuse_trials = length(poss_trials);
% n_xvtrials = round(xv_frac*(length(poss_trials)));
% xv_trials = randperm(nuse_trials);
% xv_trials(n_xvtrials+1:end) = [];
% xv_trials = poss_trials(xv_trials);
% tr_trials = setdiff(poss_trials,xv_trials);
% tr_inds = uset(ismember(RR(uset),tr_trials));
% xv_inds = uset(ismember(RR(uset),xv_trials));
% 
% temp_smooth = 10;
% poss_eye_smooth = [0.1 1 10 100 1000];
% optim_params.maxIter = 200;
% clear mod_xv TB_rate
% for pp = 1:length(poss_eye_smooth)
%     tr_mod = TB_mod;
%     tr_mod.mods(1).reg_params.lambda_d2X = poss_eye_smooth(pp);
%     tr_mod = NMMfit_filters(tr_mod,robs,TB_Xmat,[],tr_inds,0,optim_params);
%     TB_filt = reshape(tr_mod.mods(1).filtK,[length(Xtick) length(Ytick)]);
%     TB_rate{pp} = log(1+exp(TB_filt + tr_mod.spk_NL_params(1)))/base_dt;
%     mod_xv(pp) = NMMeval_model(tr_mod,robs,TB_Xmat,[],xv_inds)
% end
% 
%%
TB_dist = bsxfun(@rdivide,TB_counts,sum(TB_counts,2));
avg_dist = mean(TB_dist);
% [~,best_xv] = max(mod_xv);
% best_TB_rate = TB_rate{best_xv};
sp_mean = sum(bsxfun(@times,TB_rate,avg_dist),2);
sp_var = sum(bsxfun(@times,bsxfun(@minus,TB_rate,sp_mean).^2,avg_dist),2);

%%
back_look = 7; %look back this many time steps to parse EP trajectories

% %make an anticausal filter for processing EP history
% back_kern = zeros(back_look*2+1,1);
% back_kern(1:back_look+1) = 1;
% back_kern = flipud(back_kern/sum(back_kern));

%initialize a time-embedded version of the EPs
sp = NMMcreate_stim_params(back_look);
tbt_EP_emb = create_time_embedding(tbt_EP(:),sp);

tbt_EP_emb = tbt_EP_emb(:,4); back_look = 1;

tbt_EP_emb = reshape(tbt_EP_emb,up_nf,n_rpts,[]);

%%
rand_T = [];
tsamps = uinds(1:5:end);
for ii = 1:length(tsamps)
    cur_Dmat = abs(squareform(pdist(squeeze(tbt_EP_emb(tsamps(ii),:,:)))))/sqrt(back_look);
    cur_Dmat(logical(eye(n_rpts))) = nan;
    cur_Dmat = cur_Dmat(~isnan(cur_Dmat));
    rand_T = cat(1,rand_T,cur_Dmat);
end

n_EP_bins = 100;
EP_bin_edges = prctile(rand_T,linspace(0,100,n_EP_bins+1));
EP_bin_centers = (EP_bin_edges(1:end-1)+EP_bin_edges(2:end))/2;
%% COMPUTE VARIANCES FOR BINNED DELTA ES FOR ALL CELLS (using full recon et)
[II,JJ] = meshgrid(1:n_rpts);

cur_XC = nan(up_nf,n_EP_bins,length(targs));
cur_cnt = zeros(n_EP_bins,length(targs));
rand_XC = nan(up_nf,length(targs));
for tt = uinds
    Y1 = squeeze(tbt_binned_spikes(tt,:,:));

    cur_Dmat = abs(squareform(pdist(squeeze(tbt_EP_emb(tt,:,:)))))/sqrt(back_look);
    cur_Dmat(logical(eye(n_rpts))) = nan;
    for jj = 1:n_EP_bins
        curset = find(cur_Dmat > EP_bin_edges(jj) & cur_Dmat <= EP_bin_edges(jj+1));
        cur_XC(tt,jj,:) = squeeze(nanmean(bsxfun(@times,Y1(II(curset),:),Y1(JJ(curset),:)),1));
        cur_cnt(jj,:) = cur_cnt(jj,:) + sum(~isnan(Y1(II(curset),:)));
    end
    curset = ~isnan(cur_Dmat);
    rand_XC(tt,:) = squeeze(nanmean(bsxfun(@times,Y1(II(curset),:),Y1(JJ(curset),:)),1));
end
new_psth_var = nanmean(rand_XC);

var_ep_binned = squeeze(nanmean(cur_XC));
all_relprobs = bsxfun(@rdivide,cur_cnt,sum(cur_cnt));

%%
n_splines = 10;
spline_bounds = [0 max(EP_bin_centers)];
spline_DS = prctile(rand_T,100/n_splines:100/n_splines:(100-100/n_splines));
knot_pts = [0 0 0 0 spline_DS spline_DS(end) spline_DS(end) spline_DS(end)];

sp1 = fastBSpline.lsqspline(knot_pts,3,EP_bin_centers,squeeze(var_ep_binned(:,3)))

%% FIT SPLINES TO VAR VS ET FUNCTIONS (again using full et recon)
spline_spacing = 0.015;
spline_knots = spline_spacing:spline_spacing:maxlag_ED;
spline_eval = 0:0.005:maxlag_ED;

nboots = 0;

var_spline_ZPT = nan(n_chs,1);
var_spline_funs = nan(n_chs,length(spline_eval));
var_spline_ZPT_boot = nan(n_chs,nboots);
for ii = 1:n_chs
    x = ED_bin_centers;
    y = squeeze(var_ep_binned(:,ii));
    bad = find(isnan(y));  x(bad) = []; y(bad) = [];
    ss = fnxtr(csape(spline_knots,y(:).'/fnval(fnxtr(csape(spline_knots,eye(length(spline_knots)),'var')),x(:).'),'var'));
    var_spline_ZPT(ii) = fnval(ss,0);
    var_spline_funs(ii,:) = fnval(ss,spline_eval);
    
    npts = length(x);
    for nn = 1:nboots
        rsmp = randi(npts,npts,1);
        ss = fnxtr(csape(spline_knots,y(rsmp).'/fnval(fnxtr(csape(spline_knots,eye(length(spline_knots)),'var')),x(rsmp).'),'var'));
        var_spline_ZPT_boot(ii,nn) = fnval(ss,0);
    end
end

% psth_var_frac = ms_psth_var'./var_spline_ZPT;
new_psth_var_frac = new_psth_var'./var_spline_ZPT;
psth_var_frac_cor = ms_psth_var_cor'./var_spline_ZPT;

new_psth_var_frac_boot = bsxfun(@times,(var_spline_ZPT_boot.^(-1)),new_psth_var');
new_psth_var_frac_CIs = prctile(new_psth_var_frac_boot',[5 95]);


%%
load([model_dir '/' mod_name]);

for tt = 1:length(targs)
    cur_mod = ModData(targs(tt)).bestGQM;
    cur_block_filt = cur_mod.mods(1).filtK;
    cur_mod.spk_NL_params(1) = cur_mod.spk_NL_params(1) + mean(cur_block_filt);
    cur_mod.mods(1) = [];
    reduced_mods(tt) = cur_mod;
end

%%
full_rpt_inds = find(ismember(all_trialvec,rpt_trials));

all_stim_mat = decompressTernNoise(stimComp);

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

%%
fin_shift_cor = round(post_mean_EP/modFitParams.sp_dx);

%RECOMPUTE XMAT
best_shift_stimmat_up = all_stimmat_up;
for i=1:length(all_rpt_inds)
    best_shift_stimmat_up(used_inds(all_rpt_inds(i)),:) = shift_matrix_Nd(all_stimmat_up(used_inds(all_rpt_inds(i)),:),-fin_shift_cor(all_rpt_inds(i)),2);
end
all_Xmat_shift = create_time_embedding(best_shift_stimmat_up(full_rpt_inds,:),stim_params_full);
all_Xmat_shift = all_Xmat_shift(ismember(full_rpt_inds,used_inds),use_kInds_up);

%%
first_rpt_inds = find(all_trialvec == rpt_trials(1));
first_rpt_uinds = find(ismember(first_rpt_inds,used_inds));
base_rpt_stim = all_stimmat_up(first_rpt_inds,:);
% base_rpt_stim = create_time_embedding(all_stimmat_up(first_rpt_inds,:),stim_params_full);
% base_rpt_stim = base_rpt_stim(ismember(first_rpt_inds,used_inds),:);
%%
max_shift = round(0.5/modFitParams.sp_dx);
poss_EP_shifts = -max_shift:max_shift;

all_mod_prates = nan(length(first_rpt_uinds),length(poss_EP_shifts),length(targs));
for pp = 1:length(poss_EP_shifts)
    shifted_stim = shift_matrix_Nd(base_rpt_stim,poss_EP_shifts(pp),2);
    shifted_stim = create_time_embedding(shifted_stim,stim_params_full);
    shifted_stim = shifted_stim(first_rpt_uinds,use_kInds_up);
    for tt = 1:length(targs)
        [~,~,all_mod_prates(:,pp,tt)] = NMMmodel_eval(reduced_mods(tt),[],shifted_stim);
    end
end
all_mod_prates = all_mod_prates/modFitParams.dt;

ueps = fin_shift_cor(all_rpt_inds(abs(fin_shift_cor(all_rpt_inds)) <= max_shift));
avg_dist = hist(ueps,poss_EP_shifts); 
avg_dist = avg_dist/sum(avg_dist);
mod_sp_mean = sum(bsxfun(@times,all_mod_prates,avg_dist),2);
mod_sp_var = sum(bsxfun(@times,bsxfun(@minus,all_mod_prates,mod_sp_mean).^2,avg_dist),2);
mod_sp_mean = squeeze(mod_sp_mean);
mod_sp_var = squeeze(mod_sp_var);
%%
all_mod_emp_prates = nan(length(all_rpt_inds),length(targs));
for tt = 1:length(targs)
    [~,~,all_mod_emp_prates(:,tt)] = NMMmodel_eval(reduced_mods(tt),[],all_Xmat_shift);
end
all_mod_emp_prates = reshape(all_mod_emp_prates,length(first_rpt_uinds),n_rpts,length(targs))/modFitParams.dt;

mod_psths = squeeze(nanmean(all_mod_emp_prates,2));
mod_vars = squeeze(nanvar(all_mod_emp_prates,[],2));

uinds = (params.beg_buffer/base_dt + 1):(up_nf - params.end_buffer/base_dt);
all_mod_emp_prates = reshape(all_mod_emp_prates,[],length(targs));
cur_in_blink = in_blink_inds(uinds,:);
cur_in_sac = in_sac_inds(uinds,:);
all_mod_emp_prates(cur_in_blink(:),:) = nan;
all_mod_emp_prates(cur_in_sac(:),:) = nan;
all_mod_emp_prates = reshape(all_mod_emp_prates,length(uinds),n_rpts,length(targs));

ns_mod_psths = squeeze(nanmean(all_mod_emp_prates,2));
ns_mod_vars = squeeze(nanvar(all_mod_emp_prates,[],2));
