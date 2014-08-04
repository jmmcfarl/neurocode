clear alel
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

fig_dir = '/home/james/Analysis/bruce/saccade_modulation/ori_tuning/';

%run on 86 [0, 90];
Expt_num = 277;
Expt_name = sprintf('M%.3d',Expt_num);
if Expt_num > 280 
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end

cd(data_dir);

% load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
% load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

good_coils = [1 1];

has_data = cellfun(@(x) length(x),Expts) > 0;
expt_names(has_data) = cellfun(@(x) x.Header.expname,Expts(has_data),'UniformOutput',false);

cur_block_set = find(strcmp('rls.orXFaRC',expt_names));
% % cur_block_set = find(strcmp('rls.orRC',expt_names));
% cur_block_set = find(strcmp('rls.or',expt_names));

if Expt_num == 277
cur_block_set(cur_block_set == 26) = []; %different set of orientations used in this block
end

n_blocks = length(cur_block_set);
%%
flen = 15;
% full_nPix = 18;
% full_nPix = 9;
stim_fs = 100; %in Hz
dt = 0.01;
Fr = 1;

min_trial_dur = 0.75;

beg_buffer = 0.15;
end_buffer = 0.05;
if Expt_num < 266
    trial_dur = 2;
else
trial_dur = 4;
end
n_probes = 24;

use_right_eye = false;

n_use_blocks = Inf;

sac_backlag = round(0.1/dt);
sac_forlag = round(0.25/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_stim_or = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
for ee = 1:n_blocks;
    cur_block = cur_block_set(ee);
    
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
    if length(un_ids) < length(trial_ids)
        rpt_trials = true;
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        cur_stim_or = Expts{cur_block}.Trials(use_trials(tt)).or;
        cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        n_frames = length(cur_stim_or);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times+cur_toffset < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times(:) + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis(:) + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_stim_or = [all_stim_or; cur_stim_or(:)];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
    trial_toffset(ee) = all_t_bin_edges(end);
    cur_toffset = trial_toffset(ee);
    cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
end

orig_stim_or = all_stim_or;
%%
poss_oris = unique(all_stim_or);
poss_oris(poss_oris < -90) = [];

bin_edges = linspace(-80,90,19);
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
coarse_stim_or = all_stim_or;
for ii = 1:length(bin_cents)
    cur_set = find(all_stim_or >= bin_edges(ii) & all_stim_or < bin_edges(ii+1));
    coarse_stim_or(cur_set) = bin_cents(ii);
end
coarse_stim_or(all_stim_or >= bin_edges(end)) = bin_cents(end);

all_stim_or = coarse_stim_or;
poss_oris = bin_cents;

full_nPix = length(poss_oris);
all_stim_mat = zeros(length(all_t_axis),full_nPix);
for ii = 1:full_nPix
    cur_set = find(all_stim_or == poss_oris(ii));
    all_stim_mat(cur_set,ii) = 1;
end

all_blank_mat = zeros(length(all_t_axis),1);
all_blank_mat(orig_stim_or == -1009) = 1;
all_junk_mat = zeros(length(all_t_axis),1);
all_junk_mat(orig_stim_or == -1005) = 1;

stim_params = NMMcreate_stim_params([flen full_nPix],dt);
uflen = 1;
ustim_params = NMMcreate_stim_params([uflen full_nPix],dt);
bstim_params = NMMcreate_stim_params([flen 1],dt);

all_Xmat = create_time_embedding(all_stim_mat,stim_params);
all_Xmat_blank = create_time_embedding(all_blank_mat,bstim_params);
all_Xmat_junk = create_time_embedding(all_junk_mat,bstim_params);

% [Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
% u_ks = find(Tinds == 6);
% all_UXmat = all_Xmat(:,u_ks);

%% BIN SPIKES FOR MU AND SU
rec_type = 'LP';
clust_params.n_probes = n_probes;
if strcmp(rec_type,'LP')
    clust_params.exclude_adjacent = true;
else
    clust_params.exclude_adjacent = false;
end
[all_binned_mua,all_binned_sua,Clust_data] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;


%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_stim_times),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% PROCESS EYE TRACKING DATA
expt_bar_ori = zeros(size(cur_block_set));

% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,good_coils);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori,used_inds);

[saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params);

is_blink = detect_blinks(all_eye_ts,all_eye_vals,saccades,et_params);

[saccades,is_blink] = merge_blinks(saccades,is_blink);

sac_start_times = [saccades(:).start_time];
sac_stop_times = [saccades(:).stop_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
sac_stop_times(bad_sacs) = [];
saccades(bad_sacs) = [];
is_blink(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro);
micro_sacs = find(is_micro);

sac_durs = [saccades(:).duration];
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

sac_burst_isi = 0.15;
sacburst_set = find([saccades(used_saccade_set).isi] < sac_burst_isi | [saccades(used_saccade_set).next_isi] < sac_burst_isi);
micro_sacs(ismember(micro_sacs,sacburst_set)) = [];

%%
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

if Expt_num == 277 %problem with fixation point being in RFs with certain saccade directions
    trial_sacpos = nan(length(trial_start_inds),4);
    for ii = 1:length(trial_start_inds)
        cur_set = find(all_eye_ts >= all_trial_start_times(ii) & all_eye_ts <= all_trial_end_times(ii));
        if ~isempty(cur_set)
            trial_sacpos(ii,:) = nanmean(all_eye_vals(cur_set,:),1);
        end
    end
    trial_sac_angle = atan2(trial_sacpos(:,2),trial_sacpos(:,1))*180/pi;
    g1 = find(trial_sac_angle < -160 | trial_sac_angle > 160);
    g2 = find(trial_sac_angle > 100 & trial_sac_angle < 160);
    g3 = find(trial_sac_angle < 100 & trial_sac_angle > 40);
    g4 = find(trial_sac_angle > -160 & trial_sac_angle < -100);
    
    temp = all_tsince_start(used_inds(saccade_start_inds(big_sacs)));
    outward_sacs = big_sacs((temp >= 0.5 & temp < 1.3) | ...
        (temp >= 2.0 & temp < 2.6 | temp >= 3.5));
    inward_sacs = big_sacs((temp >= 1.3 & temp < 2.0) | ...
        (temp >= 2.6 & temp < 3.5));
    
    % exclude_ts = [g2 g4];
    exclude_ts = [g4];
    big_sacs = big_sacs(~ismember(all_trialvec(used_inds(saccade_start_inds(big_sacs))),exclude_ts));    
end
%%
saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

Xsac = zeros(NT,n_sac_bins);
Xmsac = zeros(NT,n_sac_bins);
for ii = 1:n_sac_bins
    cur_sac_target = saccade_start_inds(big_sacs) + sac_bincents(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(big_sacs(uu))) = [];
    Xsac(cur_sac_target,ii) = 1;
    
    cur_sac_target = saccade_start_inds(micro_sacs) + sac_bincents(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(micro_sacs(uu))) = [];
    Xmsac(cur_sac_target,ii) = 1;
end


%%
xv_frac = 0;
% tr_inds = 1:length(used_inds);
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

% if strcmp(xv_type,'rpt')
%     xv_trials = find(all_trial_Se==rpt_seed);
%     n_xv_trials = length(xv_trials);
% else
    n_xv_trials = round(xv_frac*nuse_trials);
    xv_trials = randperm(nuse_trials);
    xv_trials(n_xv_trials+1:end) = [];
    xv_trials = use_trials(xv_trials);
% end
tr_trials = setdiff(use_trials,xv_trials);
n_tr_trials = length(tr_trials);
fprintf('Initializing models with %d training trials and %d xval trials\n',n_tr_trials,n_xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

%%
% backlag = round(0.1/dt);
% forwardlag = round(4/dt);
% all_norm_binned_mua = bsxfun(@rdivide,all_binned_mua,nanmean(all_binned_mua));
% all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
% all_trial_start_inds(isnan(all_trial_start_inds)) = 1;
% [ov_gsac_trig_avg,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds,backlag,forwardlag);
% [ov_gsac_trig_g4,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g4),backlag,forwardlag);
% [ov_gsac_trig_g2,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g2),backlag,forwardlag);
% [ov_gsac_trig_g3,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g3),backlag,forwardlag);
% [ov_gsac_trig_g1,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g1),backlag,forwardlag);


%%
all_norm_binned_mua = bsxfun(@rdivide,all_binned_mua,nanmean(all_binned_mua));
backlag = round(0.1/dt);
forwardlag = round(0.3/dt);
[all_sac_avgs,lags] = get_event_trig_avg(all_norm_binned_mua(used_inds,:),saccade_start_inds(big_sacs),backlag,forwardlag);
[all_msac_avgs,lags] = get_event_trig_avg(all_norm_binned_mua(used_inds,:),saccade_start_inds(micro_sacs),backlag,forwardlag);

blank_locs = find(all_stim_or(used_inds) == -1009);
[all_blank_avgs,lags] = get_event_trig_avg(all_norm_binned_mua(used_inds,:),blank_locs,backlag,forwardlag);

%%
X{1} = all_Xmat(used_inds,:);
X{2} = all_Xmat_blank(used_inds,:);
X{3} = all_Xmat_junk(used_inds,:);

mod_stim_params(1) = stim_params;
mod_stim_params(2) = bstim_params;
mod_stim_params(3) = bstim_params;

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS

mod_signs = [1 1];
NL_types = {'lin','lin','lin'};
init_Xtargs = [1 2];

lambda_custom = 0;
lambda_L1 = 0;

%create L2 mat
L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,flen*full_nPix);

init_reg_params = NMMcreate_reg_params('lambda_custom',[lambda_custom; 0; 0]);

cd(anal_dir);

tot_nUnits = length(su_probes) + n_probes;
all_mod_SU = zeros(tot_nUnits,1);
all_mos_SUnum = zeros(tot_nUnits,1);
silent = 1;
for ss = 1:n_probes;
    fprintf('Computing base LLs for MU %d of %d\n',ss,n_probes);
    cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
        
%         gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
%         gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);
% 
%         [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,all_Xmat(used_inds(cur_tr_inds),:));
%         
%         gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_custom',base_lambda_d2X./var(gint)');
%         gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
%         gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);  

        gqm1 = NMMinitialize_model(mod_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds),[],[],silent,[],L2_mat);
        gqm1 = NMMfit_logexp_spkNL_withoffset(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds));

all_mu_mods(ss) = gqm1;
    end
end

%%

mod_signs = [1 1];
NL_types = {'lin','lin','lin'};
init_Xtargs = [1 2];

lambda_custom = 0;
lambda_L1 = 0;

%create L2 mat
L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,flen*full_nPix);


init_reg_params = NMMcreate_reg_params('lambda_custom',[lambda_custom; 0; 0]);
for ss = 1:length(su_probes)
    fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
    cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
        
        gqm1 = NMMinitialize_model(mod_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds),[],[],silent,[],L2_mat);
        gqm1 = NMMfit_logexp_spkNL_withoffset(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds));

        [~,~,cprate] = NMMmodel_eval(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds));
        all_su_mods(ss) = gqm1;
    end
end

%%
% new_mod = NMMfit_logexp_spkNL_withoffset(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds),[],1);
% [~,~,cprate] = NMMmodel_eval(new_mod,Robs,get_Xcell_tInds(X,cur_tr_inds));

%%
%create L2 mat
L2_params = create_L2_params([],[1 uflen*full_nPix],[uflen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,uflen*full_nPix);
% L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
% L2_mat = generate_L2_mat(L2_params,flen*full_nPix);

for ss = 1:length(su_probes)
    ss
    %%
    cur_Robs = all_binned_sua(used_inds,ss);
    cc_uinds = find(~isnan(cur_Robs));
    if ~isempty(cc_uinds)
    cur_Robs = cur_Robs(cc_uinds);
    
    cur_tr_inds = find(ismember(cc_uinds,tr_inds));
    cur_xv_inds = find(ismember(cc_uinds,xv_inds));
    
    stim_filt = reshape([all_su_mods(ss).mods(1).filtK],[flen full_nPix]);
    temp_kern = std(stim_filt,[],2);
    [~,best_lag] = max(temp_kern);
    [Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
    u_ks = find(Tinds == best_lag);
    [Xinds,Tinds] = meshgrid(1,1:flen);
    u_ks_blank = find(Tinds == best_lag);
    all_UXmat = all_Xmat(used_inds(cc_uinds),u_ks);
    all_UXmat_blank = all_Xmat_blank(used_inds(cc_uinds),u_ks_blank);
    
    gqm1 = NMMinitialize_model(ustim_params,1,NL_types,init_reg_params,init_Xtargs);
%     gqm1 = NMMinitialize_model(ustim_params,[1 -1],NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,cur_Robs,all_UXmat,[],[],silent,[],L2_mat);
%   gqm1 = all_su_mods(ss);
  
%     gqm1 = NMMinitialize_model(stim_params,1,NL_types,init_reg_params,init_Xtargs);
%     gqm1 = NMMfit_filters(gqm1,cur_Robs,all_Xmat(used_inds(cc_uinds),:),[],[],silent,[],L2_mat);

    %% FIT UPSTREAM STIM-MODULATION
    cur_Xsac = Xsac(cc_uinds,:);
%     cur_Xsac = Xmsac(cc_uinds,:);
    any_sac_inds = find(any(cur_Xsac > 0,2));
    tr_sac_inds = any_sac_inds(ismember(any_sac_inds,cur_tr_inds));
    xv_sac_inds = any_sac_inds(ismember(any_sac_inds,cur_xv_inds));
  
    %%
%     Xsac_mat = cur_Xsac(any_sac_inds,:);
%     sacGainMod = fit_sacgain_model(gqm1,cur_Robs(any_sac_inds),all_UXmat(any_sac_inds,:),Xsac_mat,lambda_d2T,lambda_L2);

    %%
    has_bar = max(all_UXmat,[],2) > 0;
    for ii = 1:length(sac_bincents)
        cur_set = find(cur_Xsac(:,ii) == 1 & has_bar);
%         cur_set = find(cur_Xsac(:,ii) == 1);
        sac_avg = sum(bsxfun(@times,all_UXmat(cur_set,:),cur_Robs(cur_set)))/sum(cur_Robs(cur_set));
        ov_avg = mean(all_UXmat(cur_set,:));
        ori_sta(ii,:) = sac_avg-ov_avg;
        ori_avg(ii,:) = mean(bsxfun(@times,all_UXmat(cur_set,:),cur_Robs(cur_set)))./mean(all_UXmat(cur_set,:));
        ori_arate(ii) = mean(cur_Robs(cur_set));

    end
    [~,mmloc] = minmax(ori_arate);
    mmloc = sort(mmloc);
    mmloc(1) = mmloc(1)-2;
    
    mlrange = diff(mmloc)+1;
    cmap = jet(mlrange);
        
    sta_ssinfo = nanmean(ori_avg.*log2(bsxfun(@rdivide,ori_avg,ori_arate')),2)./ori_arate';

    ov_avg = sum(bsxfun(@times,all_UXmat,cur_Robs))/sum(cur_Robs);
    ov_sta = ov_avg - mean(all_UXmat);
   
    all_su_stainfo(ss,:) = sta_ssinfo;
    %% FIT POST-INTEGRATION GAIN
    
    lambda_d2T = 1;
    lambda_L2 = 0.1;
    sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
    
    optim_params.optTol = 1e-6;
    optim_params.progTol = 1e-10;
    optim_params.Method = 'lbfgs';
    optim_params.verbose = 1;
    
%     [LL, penLL, pred_rate, G, gint, fgint] = NMMmodel_eval(gqm1,cur_Robs,get_Xcell_tInds(X,cc_uinds));
%     g_tot = sum(fgint,2);
    
    g_tot = all_UXmat*gqm1.mods(1).filtK;
%     g_tot = all_Xmat(used_inds(cc_uinds),:)*gqm1.mods(1).filtK;
    Xsac_tot = bsxfun(@times,cur_Xsac,g_tot);
    clear tr_stim
    tr_stim{1} = [g_tot];
    tr_stim{2} = cur_Xsac;
    tr_stim{3} = Xsac_tot;
    clear sac_stim_params
    sac_stim_params(1) = NMMcreate_stim_params(1);
    sac_stim_params(2) = NMMcreate_stim_params([size(Xsac_tot,2)]);
    sac_stim_params(3) = NMMcreate_stim_params([size(Xsac_tot,2)]);
    
    mod_signs = [1 1 1];
    Xtargets = [1 2 3];
    NL_types = {'lin','lin','lin'};
    spost_gsac_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    spost_gsac_mod = NMMfit_filters(spost_gsac_mod,cur_Robs,tr_stim,[],[],silent);
    spost_gsac_mod = NMMfit_logexp_spkNL_withoffset(spost_gsac_mod,cur_Robs,tr_stim);

    %%
%     sc_Xmat = reshape(all_Xmat(used_inds(cc_uinds),:),length(cc_uinds),flen,[]);
    sc_Xmat = reshape(all_UXmat,length(cc_uinds),[]);
%     resh_filt = reshape(gqm1.mods(1).filtK,flen,[]);
    resh_filt = gqm1.mods(1).filtK;
    g_tot = nan(length(cc_uinds),length(poss_oris));
    g_tot = bsxfun(@times,sc_Xmat,resh_filt');
%     for ii = 1:length(poss_oris)
%         g_tot(:,ii) = sc_Xmat(:,:,ii)*resh_filt(:,ii);
%     end
    
    Xsac_tot = bsxfun(@times,cur_Xsac,reshape(g_tot,[],1,length(poss_oris)));
    clear tr_stim
    tr_stim{1} = [g_tot];
    tr_stim{2} = cur_Xsac;
    tr_stim{3} = reshape(Xsac_tot,length(cc_uinds),[]);
    clear sac_stim_params
    sac_stim_params(1) = NMMcreate_stim_params(length(poss_oris));
    sac_stim_params(2) = NMMcreate_stim_params(n_sac_bins);
    sac_stim_params(3) = NMMcreate_stim_params([n_sac_bins length(poss_oris)]);
    
     sac_reg_params = NMMcreate_reg_params('lambda_d2T',50,'boundary_conds',[0 0 0]);

    mod_signs = [1 1 1];
    Xtargets = [1 2 3];
    NL_types = {'lin','lin','lin'};
    spost_gsac_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    spost_gsac_mod.mods(1).reg_params = NMMcreate_reg_params();
    spost_gsac_mod.mods(3).reg_params = NMMcreate_reg_params('lambda_d2T',5,'lambda_L2',5);
    spost_gsac_mod = NMMfit_filters(spost_gsac_mod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds),[],[],silent);

    %%
%                 [gainLL,gain_pred_rate] = eval_sacgain_mod( temp_gain_mod, cur_Robs, all_Xmat(used_inds(cc_uinds),:), cur_Xsac);

    [spost_LL,~,spost_pred_rate,~,~,~,nullLL] = NMMmodel_eval(spost_gsac_mod, cur_Robs,tr_stim);
    [sac_spost_info,sac_info,sac_arate,sac_spost_LL,sac_nullLL,sac_Nspks] = deal(nan(length(sac_bincents),1));
    for ii = 1:length(sac_bincents)
        tempr = find(cur_Xsac(:,ii) == 1);
        sac_arate(ii) = mean(cur_Robs(tempr));
        cur_avg_rate = mean(cur_Robs(tempr))*ones(size(tempr));
        sac_nullLL(ii) = nansum(cur_Robs(tempr).*log2(cur_avg_rate) - cur_avg_rate);
        sac_Nspks(ii) = sum(cur_Robs(tempr));
        
        %                         sac_LL(ii) = nansum(cur_Robs(tempr).*log2(gain_pred_rate(tempr)) - gain_pred_rate(tempr));
        %                 sac_info(ii) = nanmean(gain_pred_rate(tempr).*log2(gain_pred_rate(tempr)/mean(gain_pred_rate(tempr))))/mean(gain_pred_rate(tempr));

        sac_spost_LL(ii) = nansum(cur_Robs(tempr).*log2(spost_pred_rate(tempr)) - spost_pred_rate(tempr));
        sac_spost_info(ii) = nanmean(spost_pred_rate(tempr).*log2(spost_pred_rate(tempr)/mean(spost_pred_rate(tempr))))/mean(spost_pred_rate(tempr));
    end
    
    gsac_spostLL_info = (sac_spost_LL - sac_nullLL)./sac_Nspks;
    gsac_ov_LLinfo = (spost_LL - nullLL)/log(2);
    
    all_su_modinfo(ss,:) = sac_spost_info;
    all_su_arate(ss,:) = sac_arate;
    
    %%
%     cur_filt = reshape(gqm1.mods(1).filtK,flen,[]);
%     cur_X = reshape(all_Xmat(used_inds(cc_uinds),:),length(cc_uinds),flen,[]);
%     n_ochans = size(cur_filt,2);
%     for ii = 1:n_ochans
%        cur_chan_out(ii,:) = cur_X(:,:,ii)*cur_filt(:,ii); 
%     end
%     clear tr_stim
%     tr_stim{1} = cur_chan_out';
%     tr_stim{2} = cur_Xsac;
%     tr_stim{3} = reshape(bsxfun(@times,cur_Xsac,reshape(cur_chan_out',length(cc_uinds),1,[])),length(cc_uinds),[]);
%     cur_stim_params(1) = NMMcreate_stim_params([n_ochans]);
%     cur_stim_params(2) = NMMcreate_stim_params([n_sac_bins]);
%     cur_stim_params(3) = NMMcreate_stim_params([n_sac_bins n_ochans]);
%     cur_reg_params = NMMcreate_reg_params('lambda_d2T',2);
%     mod_signs = [1 1];
%     NL_types = {'lin','lin','lin'};
%     cur_Xtargs = [2 3];
%     cur_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,cur_reg_params,cur_Xtargs);
%     cur_mod = NMMfit_filters(cur_mod,cur_Robs,tr_stim);
    %%
n_Gbins = 15;
    
    Xtick = -(sac_backlag-1/2):(1):(sac_forlag+1/2);
            n_sbins = length(Xtick);
            
            cur_sac_starts = saccade_start_inds(big_sacs);
            cur_sac_stops = saccade_stop_inds(big_sacs);
            t_since_sac_start = nan(NT,1);
            for ii = 1:length(cur_sac_starts)
                prev_tstart = find(trial_start_inds <= cur_sac_starts(ii),1,'last');
                next_tstop = find(trial_end_inds >= cur_sac_starts(ii),1,'first');
                cur_inds = (cur_sac_starts(ii) - sac_backlag):(cur_sac_starts(ii) + sac_forlag);
                cur_uset = find(cur_inds > trial_start_inds(prev_tstart) & cur_inds < trial_end_inds(next_tstop));
                t_since_sac_start(cur_inds(cur_uset)) = sac_bincents(cur_uset);
            end
            
            TB_stim = [t_since_sac_start(cc_uinds) g_tot];
            Ytick = linspace(my_prctile(TB_stim(any_sac_inds,2),0.1),my_prctile(TB_stim(any_sac_inds,2),100-1),n_Gbins);
            TB = TentBasis2D(Xtick, Ytick);
            
            TB_lambda = 10;
            used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
                TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
            [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
            L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.025 1]);
            TB_fitmod = regGLM_fit(TB_Xmat,cur_Robs(used_data),L2_params,TB_lambda,[],[],silent);
            [LL, penLL, pred_rate, G] = regGLM_eval(TB_fitmod,cur_Robs(used_data),TB_Xmat);
            TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
            bin_areas = TB.GetBinAreas();
            gsac_TB_dist = TB_counts./bin_areas;
            gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
            gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
            sacStimProc(cc).gsac_TB_rate = gsac_TB_rate;
            
            %INFO CALS
            cur_avg_rate = mean(cur_Robs(used_data));
            marg_gdist = sum(gsac_TB_dist,2);
            marg_sdist = sum(gsac_TB_dist);
            marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
            marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
            gsacdep_info = nan(1,n_sac_bins);
            for tt = 1:n_sbins
                gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
            end
            gcumdist = cumsum(marg_gdist)/sum(marg_gdist);
            
            sacStimProc(cc).gsac_ov_TB_info = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate))/cur_avg_rate;
            
            sacStimProc(cc).gsac_TB_avg_rate = marg_gsacrate;
            sacStimProc(cc).gsac_TB_info = gsacdep_info./marg_gsacrate;
            sacStimProc(cc).gsac_TB_gdist = marg_gdist;
            sacStimProc(cc).gsac_TB_grate = marg_grate;
            
%             
%             
%    
    %%
    f1 = figure;
    subplot(3,2,1); hold on
    plot(sac_bincents*dt,ori_arate/dt,'k');
    for ii = 1:mlrange
        plot(sac_bincents(mmloc(1)-1+ii)*dt,ori_arate(mmloc(1)-1+ii)/dt,'o','color',cmap(ii,:),'linewidth',2);
    end
    xlabel('Time since sac (s)');
    ylabel('Firing rate (Hz)');
    xlim(sac_bincents([1 end])*dt);
    
    subplot(3,2,3);
    hold on
    plot(sac_bincents*dt,spost_gsac_mod.mods(2).filtK);
    hold on
    plot(sac_bincents*dt,spost_gsac_mod.mods(3).filtK,'r');
    xlabel('Time since sac (s)');
    ylabel('Filter amp');
    legend('Offset','Gain');
    xlim(sac_bincents([1 end])*dt);
    
    subplot(3,2,5); hold on
    plot(sac_bincents*dt,sac_spost_info,'r','linewidth',2);
    plot(sac_bincents*dt,sta_ssinfo,'b','linewidth',2);
    line(sac_bincents([1 end])*dt,[gsac_ov_LLinfo gsac_ov_LLinfo],'color','k')
    xlabel('Time since sac (s)');
    ylabel('SS info (bits/spk)');
    legend('Model-based','STA-based');
    xlim(sac_bincents([1 end])*dt);
    
    subplot(3,2,2)
    plot(mod(poss_oris,180),ov_sta,'k','linewidth',2);
    hold on
    for ii = 1:mlrange
        plot(mod(poss_oris,180),ori_sta(mmloc(1)-1+ii,:),'color',cmap(ii,:),'linewidth',1);
    end
    xlim([0 180]);
    xlabel('Orientation (deg)');
    ylabel('STA');
    
    subplot(3,2,4); hold on
    plot(mod(poss_oris,180),mean(ori_avg)/dt,'k','linewidth',2);
    for ii = 1:mlrange
        plot(mod(poss_oris,180),ori_avg(mmloc(1)-1+ii,:)/dt,'color',cmap(ii,:),'linewidth',1);
    end
    xlim([0 180]);
    xlabel('Orientation (deg)');
    ylabel('Firing rate (Hz)');
    
    
    
    fig_width = 8; rel_height = 1.1;
    figufy(f1);
    fname = [fig_dir sprintf('Gsac_orimod_E%d_SU%d.pdf',Expt_num,ss)];
    exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(f1);
% %         
% % %     pause
% % %     close all
    end
end

%%
%%
%create L2 mat
% L2_params = create_L2_params([],[1 uflen*full_nPix],[uflen full_nPix],2,3,[Inf -1]);
% L2_mat = generate_L2_mat(L2_params,uflen*full_nPix);
L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,flen*full_nPix);

for ss = 1:24
    ss
    %%
    cur_Robs = all_binned_mua(used_inds,ss);
    cc_uinds = find(~isnan(cur_Robs));
    if ~isempty(cc_uinds)
    cur_Robs = cur_Robs(cc_uinds);
    
    cur_tr_inds = find(ismember(cc_uinds,tr_inds));
    cur_xv_inds = find(ismember(cc_uinds,xv_inds));
    
    stim_filt = reshape([all_mu_mods(ss).mods(1).filtK],[flen full_nPix]);
    temp_kern = std(stim_filt,[],2);
    [~,best_lag] = max(temp_kern);
    [Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
    u_ks = find(Tinds == best_lag);
    [Xinds,Tinds] = meshgrid(1,1:flen);
    u_ks_blank = find(Tinds == best_lag);
    all_UXmat = all_Xmat(used_inds(cc_uinds),u_ks);
    all_UXmat_blank = all_Xmat_blank(used_inds(cc_uinds),u_ks_blank);
    
%     gqm1 = NMMinitialize_model(ustim_params,1,NL_types,init_reg_params,init_Xtargs);
%     gqm1 = NMMfit_filters(gqm1,cur_Robs,all_UXmat,[],[],silent,[],L2_mat);
%     gqm1 = NMMinitialize_model(stim_params,1,NL_types,init_reg_params,init_Xtargs);
%     gqm1 = NMMfit_filters(gqm1,cur_Robs,all_UXmat,[],[],silent,[],L2_mat);
    
gqm1 = all_mu_mods(ss);

    %% FIT UPSTREAM STIM-MODULATION
    cur_Xsac = Xsac(cc_uinds,:);
%     cur_Xsac = Xmsac(cc_uinds,:);
    any_sac_inds = find(any(cur_Xsac > 0,2));
    tr_sac_inds = any_sac_inds(ismember(any_sac_inds,cur_tr_inds));
    xv_sac_inds = any_sac_inds(ismember(any_sac_inds,cur_xv_inds));
        
    %%
%     cur_backlag = 0; cur_forlag = flen-1;
%     clear or_ta
%     for ii = 1:length(poss_oris)
%         cur_set = find(all_stim_or(used_inds(cc_uinds)) == poss_oris(ii));
%         [or_ta(ii,:),or_lags] = get_event_trig_avg_v3(cur_Robs,cur_set,cur_backlag,cur_forlag);
%     end
%     cur_set = find(all_stim_or(used_inds(cc_uinds)) == -1005);
%     blank_ta = get_event_trig_avg_v3(cur_Robs,cur_set,cur_backlag,cur_forlag);
%     
%     ffr = reshape([gqm1.mods(1).filtK],flen,[]);
%     ffr_rate = log(1+exp(ffr+gqm1.spk_NL_params(1)));
    %%
    has_bar = max(all_UXmat,[],2) > 0;
    for ii = 1:length(sac_bincents)
        cur_set = find(cur_Xsac(:,ii) == 1 & has_bar);
%         cur_set = find(cur_Xsac(:,ii) == 1);
        sac_avg = sum(bsxfun(@times,all_UXmat(cur_set,:),cur_Robs(cur_set)))/sum(cur_Robs(cur_set));
        ov_avg = mean(all_UXmat(cur_set,:));
        ori_sta(ii,:) = sac_avg-ov_avg;
        ori_avg(ii,:) = mean(bsxfun(@times,all_UXmat(cur_set,:),cur_Robs(cur_set)))./mean(all_UXmat(cur_set,:));
        ori_arate(ii) = mean(cur_Robs(cur_set));
    end
    [~,mmloc] = minmax(ori_arate);
    mlrange = diff(mmloc)+1;
    cmap = jet(mlrange);
        
    sta_ssinfo = nanmean(ori_avg.*log2(bsxfun(@rdivide,ori_avg,ori_arate')),2)./ori_arate';

    ov_avg = sum(bsxfun(@times,all_UXmat,cur_Robs))/sum(cur_Robs);
    ov_sta = ov_avg - mean(all_UXmat);
    
    all_mu_ss_info(ss,:) = sta_ssinfo;
    %% FIT POST-INTEGRATION GAIN
    
    lambda_d2T = 10;
    lambda_L2 = 0.1;
    sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
    
    optim_params.optTol = 1e-6;
    optim_params.progTol = 1e-10;
    optim_params.Method = 'lbfgs';
    optim_params.verbose = 1;
    
    [LL, penLL, pred_rate, G, gint, fgint] = NMMmodel_eval(gqm1,cur_Robs,get_Xcell_tInds(X,cc_uinds));
    g_tot = sum(fgint,2);

%     g_tot = all_UXmat*gqm1.mods(1).filtK;
    Xsac_tot = bsxfun(@times,cur_Xsac,g_tot);
    clear tr_stim
    tr_stim{1} = [g_tot];
    tr_stim{2} = cur_Xsac;
    tr_stim{3} = Xsac_tot;
    clear sac_stim_params
    sac_stim_params(1) = NMMcreate_stim_params(1);
    sac_stim_params(2) = NMMcreate_stim_params([size(Xsac_tot,2)]);
    sac_stim_params(3) = NMMcreate_stim_params([size(Xsac_tot,2)]);
    
    mod_signs = [1 1 1];
    Xtargets = [1 2 3];
    NL_types = {'lin','lin','lin'};
    spost_gsac_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    spost_gsac_mod = NMMfit_filters(spost_gsac_mod,cur_Robs,tr_stim,[],[],silent);
    spost_gsac_mod = NMMfit_logexp_spkNL_withoffset(spost_gsac_mod,cur_Robs,tr_stim);
    
    %%
    [spost_LL,~,spost_pred_rate,~,~,~,nullLL] = NMMmodel_eval(spost_gsac_mod, cur_Robs,tr_stim);
    [sac_spost_info,sac_info,sac_arate,sac_spost_LL,sac_nullLL,sac_Nspks] = deal(nan(length(sac_bincents),1));
    for ii = 1:length(sac_bincents)
        tempr = find(cur_Xsac(:,ii) == 1);
        sac_arate(ii) = mean(cur_Robs(tempr));
        cur_avg_rate = mean(cur_Robs(tempr))*ones(size(tempr));
        sac_nullLL(ii) = nansum(cur_Robs(tempr).*log2(cur_avg_rate) - cur_avg_rate);
        sac_Nspks(ii) = sum(cur_Robs(tempr));
        
        sac_spost_LL(ii) = nansum(cur_Robs(tempr).*log2(spost_pred_rate(tempr)) - spost_pred_rate(tempr));
        sac_spost_info(ii) = nanmean(spost_pred_rate(tempr).*log2(spost_pred_rate(tempr)/mean(spost_pred_rate(tempr))))/mean(spost_pred_rate(tempr));
    end
    
    gsac_spostLL_info = (sac_spost_LL - sac_nullLL)./sac_Nspks;
    gsac_ov_LLinfo = (spost_LL - nullLL)/log(2);
    
    all_mu_info(ss,:) = sac_spost_info;
    all_mu_arate(ss,:) = sac_arate;
    
    %%
    f1 = figure;
    subplot(3,2,1); hold on
    plot(sac_bincents*dt,ori_arate/dt,'k');
    for ii = 1:mlrange
        plot(sac_bincents(mmloc(1)-1+ii)*dt,ori_arate(mmloc(1)-1+ii)/dt,'o','color',cmap(ii,:),'linewidth',2);
    end
    xlabel('Time since sac (s)');
    ylabel('Firing rate (Hz)');
    xlim(sac_bincents([1 end])*dt);
    
    subplot(3,2,3);
    hold on
    plot(sac_bincents*dt,spost_gsac_mod.mods(2).filtK);
    hold on
    plot(sac_bincents*dt,spost_gsac_mod.mods(3).filtK,'r');
    xlabel('Time since sac (s)');
    ylabel('Filter amp');
    legend('Offset','Gain');
    xlim(sac_bincents([1 end])*dt);
    
    subplot(3,2,5); hold on
    plot(sac_bincents*dt,sac_spost_info,'r','linewidth',2);
%     plot(sac_bincents*dt,sta_ssinfo,'b','linewidth',2);
%     line(sac_bincents([1 end])*dt,[gsac_ov_LLinfo gsac_ov_LLinfo],'color','k')
    xlabel('Time since sac (s)');
    ylabel('SS info (bits/spk)');
    legend('Model-based','STA-based');
    xlim(sac_bincents([1 end])*dt);
    
    subplot(3,2,2)
    plot(mod(poss_oris,180),ov_sta,'k','linewidth',2);
    hold on
    for ii = 1:mlrange
        plot(mod(poss_oris,180),ori_sta(mmloc(1)-1+ii,:),'color',cmap(ii,:),'linewidth',1);
    end
    xlim([0 180]);
    xlabel('Orientation (deg)');
    ylabel('STA');
    
    subplot(3,2,4); hold on
    plot(mod(poss_oris,180),mean(ori_avg)/dt,'k','linewidth',2);
    for ii = 1:mlrange
        plot(mod(poss_oris,180),ori_avg(mmloc(1)-1+ii,:)/dt,'color',cmap(ii,:),'linewidth',1);
    end
    xlim([0 180]);
    xlabel('Orientation (deg)');
    ylabel('Firing rate (Hz)');
    
    
    
    fig_width = 8; rel_height = 1.1;
    figufy(f1);
    fname = [fig_dir sprintf('Gsac_orimod_flen_E%d_MU%d.pdf',Expt_num,ss)];
    exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(f1);
        
% %     pause
% %     close all
    end
end

%%
t_ax = (0:flen)*dt;
close all
% for cc = 17
cc = 17;
 f1 = figure();
 subplot(2,2,1)
    plot(sac_bincents*dt,all_mu_arate(cc,:)/dt)
    xlabel('Time (s)');
    ylabel('Rate (Hz)');
    subplot(2,2,3)
    plot(sac_bincents*dt,all_mu_info(cc,:)); hold on
    plot(sac_bincents*dt,all_mu_ss_info(cc,:),'r')
    xlabel('Time (s)');
    ylabel('SS-info bits/spk');
    
    subplot(2,2,[2]);
    ffr = reshape(all_mu_mods(cc).mods(1).filtK,[flen 18]);
    imagesc(poss_oris,t_ax, ffr); set(gca,'ydir','normal'); colorbar
    xlabel('Orientation (deg)');
    ylabel('Lag (s)');
    caxis([-1.5 1.5])
    subplot(2,2,[4]);
    ffr = reshape(all_mu_mods(cc).mods(1).filtK,[flen 18]);
    imagesc(poss_oris,t_ax, ffr); set(gca,'ydir','normal'); colorbar
    xlabel('Orientation (deg)');
    ylabel('Lag (s)');
    caxis([-0.5 0.5])
    cc
    
            fig_width = 8; rel_height = 0.7;
    figufy(f1);
    fname = [fig_dir sprintf('Gsac_orimod_examp_MU%d.pdf',cc)];
    exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%     pause
%     clf
% end

%%
% t_ax = (0:flen)*dt;
% close all
% for cc = 1
cc = 1;
    
f1 =figure();
subplot(2,2,1)
    plot(sac_bincents*dt,all_su_arate(cc,:)/dt)
    xlabel('Time (s)');
    ylabel('Rate (Hz)');
    subplot(2,2,3)
    plot(sac_bincents*dt,all_su_modinfo(cc,:)); hold on
    plot(sac_bincents*dt,all_su_stainfo(cc,:),'r')
    xlabel('Time (s)');
    ylabel('SS-info bits/spk');
    caxis([-2 2])
    subplot(2,2,[2]);
    ffr = reshape(all_su_mods(cc).mods(1).filtK,[flen 18]);
    imagesc(poss_oris,t_ax, ffr); set(gca,'ydir','normal'); colorbar
    xlabel('Orientation (deg)');
    ylabel('Lag (s)');
    subplot(2,2,[4]);
    ffr = reshape(all_su_mods(cc).mods(1).filtK,[flen 18]);
    imagesc(poss_oris,t_ax, ffr); set(gca,'ydir','normal'); colorbar
    xlabel('Orientation (deg)');
    ylabel('Lag (s)');
    caxis([-0.5 0.5]);
    
    
        fig_width = 8; rel_height = 0.7;
    figufy(f1);
    fname = [fig_dir sprintf('Gsac_orimod_examp_SU%d.pdf',cc)];
    exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

    cc
%     pause
%     clf
% end