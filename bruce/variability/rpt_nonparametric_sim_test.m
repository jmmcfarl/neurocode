clear all
close all

addpath('~/other_code/fastBSpline/');

global Expt_name bar_ori monk_name rec_type rec_number

Expt_name = 'M012';
monk_name = 'jbe';
bar_ori = 0; %bar orientation to use (only for UA recs)
rec_number = 1;

use_hres_ET = true; EP_params.use_hres_ET = use_hres_ET; %use high-res eye-tracking?
exclude_sacs = true; EP_params.exclude_sacs = exclude_sacs;
sub_trialavgs = true; EP_params.sub_trialavgs = sub_trialavgs; %subtract out trial avg spike counts

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

sname = 'rpt_variability_compact_multDT';
if sub_trialavgs
    sname = strcat(sname,'_subTrial');
end
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


in_sac_inds = reshape(in_sac_inds,used_nf,n_trials);
in_blink_inds = reshape(in_blink_inds,used_nf,n_trials);
in_sac_inds(isnan(in_sac_inds)) = 0; in_blink_inds(isnan(in_blink_inds)) = 0;
in_sac_inds = logical(in_sac_inds); in_blink_inds = logical(in_blink_inds);

tbt_EP = reshape(post_mean_EP,used_nf,n_trials);
used_NF = (params.beg_buffer/params.dt+1):(params.trial_dur - params.end_buffer)/params.dt;

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

%%
poss_ntrials = [25 50 100 150 200];
max_sim_trials = max(poss_ntrials);
n_sim_rpts = 25;

ex_mod = 3;
cur_NMM = stim_mod(ex_mod);
for rr = 1:n_sim_rpts
    fprintf('Sim rpt %d/%d\n',rr,n_sim_rpts);
    
    sim_trial_set = randperm(n_trials); %pick a subset of total trials for EP data
    sim_trial_set = sim_trial_set(1:max_sim_trials);
    sim_tbt_EP = tbt_EP(:,sim_trial_set);
    
    fin_shift_cor = round(sim_tbt_EP/modFitParams.sp_dx); %use overall EP estimate
    fin_shift_cor(isnan(fin_shift_cor)) = 0;
    
    rand_trial_stim = randi(n_trials,1,1);
    cur_trial_inds = find(all_trialvec == rand_trial_stim);
    
    base_rpt_stim = all_stimmat_up(cur_trial_inds,:);
    
    %RECOMPUTE XMAT
    best_shift_stimmat_up = repmat(base_rpt_stim,[1 1 max_sim_trials]);
    for ii = 1:max_sim_trials
        for jj = 1:length(used_NF)
            best_shift_stimmat_up(used_NF(jj),:,ii) = shift_matrix_Nd(base_rpt_stim(used_NF(jj),:),-fin_shift_cor(jj,ii),2);
        end
    end
    best_shift_stimmat_up = permute(best_shift_stimmat_up,[1 3 2]);
    best_shift_stimmat_up = reshape(best_shift_stimmat_up,[],full_nPix_us);
    
    all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_full);
    all_Xmat_shift = reshape(all_Xmat_shift,[],max_sim_trials,full_nPix_us*modFitParams.flen);
    all_Xmat_shift = reshape(all_Xmat_shift(used_NF,:,use_kInds_up),[],length(use_kInds_up));
    
    %%
    all_mod_emp_prates = nan(size(all_Xmat_shift,1),length(targs));
    for cc = 1:length(targs)
        if has_stim_mod(cc)
            [~,~,all_mod_emp_prates(:,cc)] = NMMmodel_eval(stim_mod(cc),[],all_Xmat_shift);
        end
    end
    if exclude_blinks
    all_mod_emp_prates(reshape(in_blink_inds(:,sim_trial_set),[],1),:) = nan;
    end
    if exclude_sacs
        all_mod_emp_prates(reshape(in_sac_inds(:,sim_trial_set),[],1),:) = nan;
    end
    all_mod_emp_prates = reshape(all_mod_emp_prates,used_nf,max_sim_trials,length(targs));
    
    %% Construct Time embedded eye position sequences for repeat trials
    emb_win = (bin_dt + 0.05); EP_params.emb_win = emb_win; %look back this many time steps to parse EP trajectories
    emb_shift = 0.03; EP_params.emb_shift = emb_shift;
    emb_win = round(emb_win/orig_dt);
    emb_shift = round(emb_shift/orig_dt);
    
    %initialize a time-embedded version of the EPs
    sp = NMMcreate_stim_params(emb_win + emb_shift);
    tbt_EP_emb = create_time_embedding(sim_tbt_EP(:),sp);
    if exclude_blinks
    tbt_EP_emb(reshape(in_blink_inds(:,sim_trial_set),[],1),:) = nan;
    end
    if exclude_sacs
        tbt_EP_emb(reshape(in_sac_inds(:,sim_trial_set),[],1),:) = nan;
    end
    tbt_EP_emb = reshape(tbt_EP_emb(:,(emb_shift+1):end),used_nf,max_sim_trials,[]);
    
    %% ESTIMATE OVERALL DISTRIBUTION OF DELTA_X to determine quantiles
    
    % compute the distribution of delta_X
    rand_T = [];
    
    % estimate quantiles of the distribution of pairwise EP similarities by this metric
    rset = randi(length(used_NF),100,1);
    for ii = 1:length(rset)
        cur_Dmat = abs(squareform(pdist(squeeze(tbt_EP_emb(rset(ii),:,:)))))/sqrt(emb_win);
        cur_Dmat(logical(eye(length(cur_trial_set)))) = nan;
        cur_Dmat = cur_Dmat(~isnan(cur_Dmat));
        rand_T = cat(1,rand_T,cur_Dmat);
    end
    
    %%
    
    for pp = 1:length(poss_ntrials)
        use_ntrials = poss_ntrials(pp);
        
        cur_trial_set = randperm(max_sim_trials); cur_trial_set = cur_trial_set(1:use_ntrials);
        cur_tbt_prates = squeeze(all_mod_emp_prates(:,cur_trial_set,ex_mod));
        
        simStats(rr,pp).acrossTrialVar = mean(nanvar(cur_tbt_prates,[],2));
        simStats(rr,pp).psthVar = var(nanmean(cur_tbt_prates,2));
        simStats(rr,pp).totVar = nanvar(cur_tbt_prates(:));
        
        simSpikes = poissrnd(cur_tbt_prates);
        
        %%
        %subtract out trial-avg spk counts if desired
        trial_avg_spikes = nanmean(simSpikes);
        if sub_trialavgs
            simSpikes = bsxfun(@minus,simSpikes,trial_avg_spikes);
        end
        
        %now subtract out overall avg spike count
        simSpikes = simSpikes - nanmean(simSpikes(:));
        
        %% MAIN WITHIN-CELL ANALYSIS LOOP
        maxD_prc = 50;
        
        poss_n_splines = [3:8]; EP_params.poss_N_splines = poss_n_splines; %range of possible values for number of splines.
        best_n_knots = 4; EP_params.best_n_knots = 4;
        n_EP_bins = 50; EP_params.n_EP_bins = n_EP_bins;
        EP_bin_edges = prctile(rand_T,linspace(0,maxD_prc,n_EP_bins+1));
        EP_bin_centers = (EP_bin_edges(1:end-1)+EP_bin_edges(2:end))/2;  EP_params.EP_bin_centers = EP_bin_centers;
        maxD = prctile(rand_T,maxD_prc);
        
        poss_eps_sizes = [.005 .01 .02]; EP_params.poss_eps_sizes = poss_eps_sizes;
        
        n_eval_pts = 100; EP_params.n_eval_pts = n_eval_pts; %number of points to evaluate spline fit
        eval_xx = unique([0 prctile(rand_T,linspace(0,maxD_prc,n_eval_pts))]); %x-axis for evaluating spline models
        EP_params.eval_xx = eval_xx;
        
        all_D = [];
        all_X = [];
        [II,JJ] = meshgrid(1:use_ntrials);
        uset = JJ > II; %only need to count each unique trial pair once
        n_unique_pairs = sum(uset(:));
        
        cur_D = nan(n_unique_pairs*used_nf,1);
        cur_X = nan(n_unique_pairs*used_nf,1);
        for tt = 1:used_nf
            cur_inds = (tt-1)*n_unique_pairs + (1:n_unique_pairs);
            Y1 = squeeze(simSpikes(tt,:));
            
            cur_Dmat = abs(squareform(pdist(squeeze(tbt_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
            cur_Dmat(logical(eye(use_ntrials))) = nan;
            cur_D(cur_inds) = cur_Dmat(uset);
            
            cur_Xmat = bsxfun(@times,Y1(II),Y1(JJ));
            cur_X(cur_inds) = cur_Xmat(uset);
        end
        
        cur_upts = find(~isnan(cur_D) & ~isnan(cur_X));
        all_D = cat(1,all_D,cur_D(cur_upts));
        all_X = cat(1,all_X,cur_X(cur_upts));
        
        n_data_points = length(all_X);
        
        [bincnts,binids] = histc(all_D,EP_bin_edges);
        var_ep_binned = nan(n_EP_bins,1);
        for bb = 1:n_EP_bins
            var_ep_binned(bb) = mean(all_X(binids == bb));
        end
        
        spline_DS = prctile(all_D,maxD_prc/(best_n_knots-1):maxD_prc/(best_n_knots-1):(maxD_prc-maxD_prc/(best_n_knots-1)));
        knot_pts = [0 0 0 0 spline_DS maxD maxD maxD];
        upts = find(all_D <= knot_pts(end));
        
        eps_ball_var = nan(length(poss_eps_sizes),1);
        eps_ball_cnts = nan(length(poss_eps_sizes),1);
        for bb = 1:length(poss_eps_sizes)
            curset = find(all_D < poss_eps_sizes(bb));
            eps_ball_var(bb) = mean(all_X(curset));
            eps_ball_cnts(bb) = length(curset);
        end
        
        sp = fastBSpline.lsqspline(knot_pts,3,all_D(upts),all_X(upts));
        spline_pred = sp.evalAt(eval_xx);
        simStats(rr,pp).spline_var = spline_pred(1);
        simStats(rr,pp).eps_vars = eps_ball_var;
        
        simStats(rr,pp).pair_psth_var = mean(all_X);
        
    end
end

%%
atVars = arrayfun(@(x) x.acrossTrialVar,simStats);
totVars = arrayfun(@(x) x.totVar,simStats);
psthVars = arrayfun(@(x) x.psthVar,simStats);
noiseVars = bsxfun(@rdivide,atVars,poss_ntrials);

pairVar = arrayfun(@(x) x.pair_psth_var,simStats);
splineVar = arrayfun(@(x) x.spline_var,simStats);

totVarErr = (splineVar-totVars)./totVars*100;
psthVarErr = (pairVar-psthVars)./psthVars*100;
psthVarErr_cor = (pairVar-(psthVars-noiseVars))./(psthVars-noiseVars)*100;

%look at bias and variance of psth and total variance estimates using
%trial-pairs
f1 = figure();
subplot(2,1,1); hold on
errorbar(poss_ntrials,nanmean(psthVarErr),nanstd(psthVarErr));
errorbar(poss_ntrials,nanmean(psthVarErr_cor),nanstd(psthVarErr_cor),'r');
subplot(2,1,2); 
errorbar(poss_ntrials,nanmean(totVarErr),nanstd(totVarErr));

%look at bias and variance of total variance decomposition
totErrs = (totVars-(atVars+psthVars))./totVars*100;
totErrs_cor = (totVars-(atVars+psthVars-noiseVars))./totVars*100;
f2 = figure(); hold on
errorbar(poss_ntrials,nanmean(totErrs),nanstd(totErrs));
errorbar(poss_ntrials,nanmean(totErrs_cor),nanstd(totErrs_cor),'r');