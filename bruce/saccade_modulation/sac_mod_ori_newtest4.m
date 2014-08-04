clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

fig_dir = '/home/James/Analysis/bruce/saccade_modulation/ori_tuning/';

%run on 86 [0, 90];
Expt_num = 296;
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
flen = 20;
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

%%
poss_oris = unique(all_stim_or);
poss_oris(poss_oris < -90) = [];

n_coarse_bins = 9;
if Expt_num == 277
    bin_edges = linspace(-80,90,n_coarse_bins+1);
else
    bin_edges = linspace(0,180,n_coarse_bins+1);   
end
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);

coarse_stim_or = all_stim_or;
for ii = 1:length(bin_cents)
    cur_set = find(all_stim_or >= bin_edges(ii) & all_stim_or < bin_edges(ii+1));
    coarse_stim_or(cur_set) = bin_cents(ii);
end
coarse_stim_or(all_stim_or >= bin_edges(end)) = bin_cents(end);

% all_stim_or = coarse_stim_or;
% poss_oris = bin_cents;

full_nPix = length(poss_oris);
all_stim_mat = zeros(length(all_t_axis),full_nPix);
for ii = 1:full_nPix
    cur_set = find(all_stim_or == poss_oris(ii));
    all_stim_mat(cur_set,ii) = 1;
end

all_blank_mat = zeros(length(all_t_axis),1);
all_blank_mat(all_stim_or == -1009) = 1;
all_junk_mat = zeros(length(all_t_axis),1);
all_junk_mat(all_stim_or == -1005) = 1;

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
% saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));
saccade_trial_inds = all_trialvec(used_inds(saccade_stop_inds));

Xsac = zeros(NT,n_sac_bins);
Xmsac = zeros(NT,n_sac_bins);
for ii = 1:n_sac_bins
%     cur_sac_target = saccade_start_inds(big_sacs) + sac_bincents(ii);
    cur_sac_target = saccade_stop_inds(big_sacs) + sac_bincents(ii);
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

silent = 0;
if Expt_num == 277
    mod_signs = [1 1 1];
    init_Xtargs = [1 2 3];
else
    mod_signs = [1 1];
    init_Xtargs = [1 2];
end
NL_types = {'lin','lin','lin'};

lambda_custom = 0;
lambda_L2 = 0;
lambda_d2T = lambda_custom;

%create L2 mat
L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,flen*full_nPix);

other_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2);
init_reg_params = NMMcreate_reg_params('lambda_custom',[lambda_custom; 0; 0],'lambda_L2',[lambda_L2; 0; 0]);

cd(anal_dir);

tot_nUnits = length(su_probes) + n_probes;
all_mod_SU = zeros(tot_nUnits,1);
all_mos_SUnum = zeros(tot_nUnits,1);
for ss = 1:n_probes;
    fprintf('Computing base LLs for MU %d of %d\n',ss,n_probes);
    cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
        
        gqm1 = NMMinitialize_model(mod_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        for jj = 2:length(mod_signs)
            gqm1.mods(jj).reg_params = other_reg_params;
        end
        
        gqm1 = NMMfit_filters(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds),[],[],silent,[],L2_mat);
        gqm1 = NMMfit_logexp_spkNL_withoffset(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds));
        
        all_mu_mods(ss) = gqm1;
    end
end

%%
silent = 1;
NL_types = {'lin','lin','lin'};

if Expt_num == 277
%     mod_signs = [1 1 1];
%     init_Xtargs = [1 2 3];
    mod_signs = [1 1 1];
    init_Xtargs = [1 2 3];
else
    mod_signs = [1 1];
    init_Xtargs = [1 2];
end
lambda_custom = 50;
lambda_L2 = 500;
lambda_d2T = lambda_custom;

%create L2 mat
L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,flen*full_nPix);

other_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
init_reg_params = NMMcreate_reg_params('lambda_custom',[lambda_custom],'lambda_L2',[lambda_L2]);
% for ss = 1:length(su_probes)
ss=1;
fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
    cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
    cur_xv_inds = xv_inds(~isnan(all_binned_sua(used_inds(xv_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
        
        gqm1 = NMMinitialize_model(mod_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        for jj = 2:length(mod_signs)
            gqm1.mods(jj).reg_params = other_reg_params;
        end
        gqm1 = NMMfit_filters(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds),[],[],silent,[],L2_mat);
        gqm1 = NMMfit_logexp_spkNL(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds));
        
        xvLL = NMMmodel_eval(gqm1,all_binned_sua(used_inds(cur_xv_inds),ss),get_Xcell_tInds(X,cur_xv_inds));
        
        [~,~,cprate] = NMMmodel_eval(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds));
        all_su_mods(ss) = gqm1;
    end
% end

%%
blag = 2;
flag = 20;
all_SU_tavgs = nan(length(su_probes),length(poss_oris),blag+flag+1);
for ss = 1:length(su_probes)
    fprintf('SU %d of %d\n',ss,length(su_probes));
    cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    cur_Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
    if ~isempty(cur_Robs)
        for ii = 1:length(poss_oris)
            [all_SU_tavgs(ss,ii,:),ori_lags] = get_event_trig_avg_v3(cur_Robs,find(all_stim_or(used_inds(cc_uinds)) == poss_oris(ii)),blag,flag);
        end
    end
    su_avgrate(ss) = mean(cur_Robs);
end
all_SU_Navgs = bsxfun(@minus,all_SU_tavgs,su_avgrate');
all_SU_NRavgs = bsxfun(@rdivide,all_SU_Navgs,std(all_SU_Navgs,[],3));
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
    gqm1 = NMMfit_filters(gqm1,cur_Robs,all_UXmat,[],[],silent,[],L2_mat);

    %% FIT UPSTREAM STIM-MODULATION
    cur_Xsac = Xsac(cc_uinds,:);
%     cur_Xsac = Xmsac(cc_uinds,:);
    any_sac_inds = find(any(cur_Xsac > 0,2));
    tr_sac_inds = any_sac_inds(ismember(any_sac_inds,cur_tr_inds));
    xv_sac_inds = any_sac_inds(ismember(any_sac_inds,cur_xv_inds));
  

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
    mmloc(1) = mmloc(1)-1;
    
    mlrange = diff(mmloc)+1;
    cmap = jet(mlrange);
        
    sta_ssinfo = nanmean(ori_avg.*log2(bsxfun(@rdivide,ori_avg,ori_arate')),2)./ori_arate';

    ov_avg = sum(bsxfun(@times,all_UXmat,cur_Robs))/sum(cur_Robs);
    ov_sta = ov_avg - mean(all_UXmat);
   
    all_su_stainfo(ss,:) = sta_ssinfo;

    %%
% %     sc_Xmat = reshape(all_Xmat(used_inds(cc_uinds),:),length(cc_uinds),flen,[]);
%     sc_Xmat = reshape(all_UXmat,length(cc_uinds),[]);
% %     resh_filt = reshape(gqm1.mods(1).filtK,flen,[]);
%     resh_filt = gqm1.mods(1).filtK;
%     g_tot = nan(length(cc_uinds),length(poss_oris));
%     g_tot = bsxfun(@times,sc_Xmat,resh_filt');
% %     for ii = 1:length(poss_oris)
% %         g_tot(:,ii) = sc_Xmat(:,:,ii)*resh_filt(:,ii);
% %     end
%     
%     Xsac_tot = bsxfun(@times,cur_Xsac,reshape(g_tot,[],1,length(poss_oris)));
%     clear tr_stim
%     tr_stim{1} = [g_tot];
%     tr_stim{2} = cur_Xsac;
%     tr_stim{3} = reshape(Xsac_tot,length(cc_uinds),[]);
%     clear sac_stim_params
%     sac_stim_params(1) = NMMcreate_stim_params(length(poss_oris));
%     sac_stim_params(2) = NMMcreate_stim_params(n_sac_bins);
%     sac_stim_params(3) = NMMcreate_stim_params([n_sac_bins length(poss_oris)]);
%     
%      sac_reg_params = NMMcreate_reg_params('lambda_d2T',50,'boundary_conds',[0 0 0]);
% 
%     mod_signs = [1 1 1];
%     Xtargets = [1 2 3];
%     NL_types = {'lin','lin','lin'};
%     spost_gsac_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%     spost_gsac_mod.mods(1).reg_params = NMMcreate_reg_params();
%     spost_gsac_mod.mods(3).reg_params = NMMcreate_reg_params('lambda_d2T',5,'lambda_L2',0.1);
%     spost_gsac_mod = NMMfit_filters(spost_gsac_mod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds),[],[],silent);
% 
%     %%
% %                 [gainLL,gain_pred_rate] = eval_sacgain_mod( temp_gain_mod, cur_Robs, all_Xmat(used_inds(cc_uinds),:), cur_Xsac);
% 
%     [spost_LL,~,spost_pred_rate,~,~,~,nullLL] = NMMmodel_eval(spost_gsac_mod, cur_Robs,tr_stim);
%     [sac_spost_info,sac_info,sac_arate,sac_spost_LL,sac_nullLL,sac_Nspks] = deal(nan(length(sac_bincents),1));
%     for ii = 1:length(sac_bincents)
%         tempr = find(cur_Xsac(:,ii) == 1);
%         sac_arate(ii) = mean(cur_Robs(tempr));
%         cur_avg_rate = mean(cur_Robs(tempr))*ones(size(tempr));
%         sac_nullLL(ii) = nansum(cur_Robs(tempr).*log2(cur_avg_rate) - cur_avg_rate);
%         sac_Nspks(ii) = sum(cur_Robs(tempr));
%         
%         %                         sac_LL(ii) = nansum(cur_Robs(tempr).*log2(gain_pred_rate(tempr)) - gain_pred_rate(tempr));
%         %                 sac_info(ii) = nanmean(gain_pred_rate(tempr).*log2(gain_pred_rate(tempr)/mean(gain_pred_rate(tempr))))/mean(gain_pred_rate(tempr));
% 
%         sac_spost_LL(ii) = nansum(cur_Robs(tempr).*log2(spost_pred_rate(tempr)) - spost_pred_rate(tempr));
%         sac_spost_info(ii) = nanmean(spost_pred_rate(tempr).*log2(spost_pred_rate(tempr)/mean(spost_pred_rate(tempr))))/mean(spost_pred_rate(tempr));
%     end
%     
%     gsac_spostLL_info = (sac_spost_LL - sac_nullLL)./sac_Nspks;
%     gsac_ov_LLinfo = (spost_LL - nullLL)/log(2);
%     
%     all_su_modinfo(ss,:) = sac_spost_info;
%     all_su_arate(ss,:) = sac_arate;
    
    %%
    flen_tax = (0:(flen-1))*dt + dt/2;
    ffr = reshape([all_su_mods(ss).mods(1).filtK],flen,[]);
    tkern = sqrt(sum(ffr.^2,2));
    ffrn = bsxfun(@rdivide,ffr,tkern);

    f1 = figure;
    
    subplot(2,3,1)
    imagesc(flen_tax,poss_oris,flipud(ffr'));
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    axis tight
    xlabel('Time lag (s)');
    ylabel('Orientation');
    
    subplot(2,3,4)
    imagesc(flen_tax,poss_oris,flipud(ffrn'));
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    axis tight
    xlabel('Time lag (s)');
    ylabel('Orientation');

    subplot(2,3,2); hold on
    plot(sac_bincents*dt,ori_arate/dt,'k');
    for ii = 1:mlrange
        plot(sac_bincents(mmloc(1)-1+ii)*dt,ori_arate(mmloc(1)-1+ii)/dt,'o','color',cmap(ii,:),'linewidth',2);
    end
    xlabel('Time since sac (s)');
    ylabel('Firing rate (Hz)');
    xlim([-0.1 0.15])
    
    
    subplot(2,3,5); hold on
    plot(sac_bincents*dt,sta_ssinfo,'b','linewidth',2);
%     line(sac_bincents([1 end])*dt,[gsac_ov_LLinfo gsac_ov_LLinfo],'color','k')
    xlabel('Time since sac (s)');
    ylabel('SS info (bits/spk)');
    xlim([-0.1 0.15])
        
    subplot(2,3,3); hold on
    imagesc(sac_bincents*dt,poss_oris,ori_avg');
    axis tight
    xlim([-0.1 0.15])
    ylabel('Orientation (deg)');
    xlabel('Time since sac (s)');
    title('Absolute rate');
    
    subplot(2,3,6); hold on
    imagesc(sac_bincents*dt,poss_oris,bsxfun(@rdivide,ori_avg,mean(ori_avg))');
    axis tight
    xlim([-0.1 0.15])
    ylabel('Orientation (deg)');
    xlabel('Time since sac (s)');
    title('Relative rate');
    
    
    fig_width = 11; rel_height = 0.6;
    
    figufy(f1);
    fname = [fig_dir sprintf('GsacE_neworimod4_E%d_SU%d.pdf',Expt_num,ss)];
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
L2_params = create_L2_params([],[1 uflen*full_nPix],[uflen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,uflen*full_nPix);
% L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
% L2_mat = generate_L2_mat(L2_params,flen*full_nPix);

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
    
    gqm1 = NMMinitialize_model(ustim_params,1,NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,cur_Robs,all_UXmat,[],[],silent,[],L2_mat);

    %% FIT UPSTREAM STIM-MODULATION
    cur_Xsac = Xsac(cc_uinds,:);
%     cur_Xsac = Xmsac(cc_uinds,:);
    any_sac_inds = find(any(cur_Xsac > 0,2));
    tr_sac_inds = any_sac_inds(ismember(any_sac_inds,cur_tr_inds));
    xv_sac_inds = any_sac_inds(ismember(any_sac_inds,cur_xv_inds));
  

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
   
    all_mu_stainfo(ss,:) = sta_ssinfo;

    %%
% %     sc_Xmat = reshape(all_Xmat(used_inds(cc_uinds),:),length(cc_uinds),flen,[]);
%     sc_Xmat = reshape(all_UXmat,length(cc_uinds),[]);
% %     resh_filt = reshape(gqm1.mods(1).filtK,flen,[]);
%     resh_filt = gqm1.mods(1).filtK;
%     g_tot = nan(length(cc_uinds),length(poss_oris));
%     g_tot = bsxfun(@times,sc_Xmat,resh_filt');
% %     for ii = 1:length(poss_oris)
% %         g_tot(:,ii) = sc_Xmat(:,:,ii)*resh_filt(:,ii);
% %     end
%     
%     Xsac_tot = bsxfun(@times,cur_Xsac,reshape(g_tot,[],1,length(poss_oris)));
%     clear tr_stim
%     tr_stim{1} = [g_tot];
%     tr_stim{2} = cur_Xsac;
%     tr_stim{3} = reshape(Xsac_tot,length(cc_uinds),[]);
%     clear sac_stim_params
%     sac_stim_params(1) = NMMcreate_stim_params(length(poss_oris));
%     sac_stim_params(2) = NMMcreate_stim_params(n_sac_bins);
%     sac_stim_params(3) = NMMcreate_stim_params([n_sac_bins length(poss_oris)]);
%     
%      sac_reg_params = NMMcreate_reg_params('lambda_d2T',50,'boundary_conds',[0 0 0]);
% 
%     mod_signs = [1 1 1];
%     Xtargets = [1 2 3];
%     NL_types = {'lin','lin','lin'};
%     spost_gsac_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%     spost_gsac_mod.mods(1).reg_params = NMMcreate_reg_params();
%     spost_gsac_mod.mods(3).reg_params = NMMcreate_reg_params('lambda_d2T',5,'lambda_L2',0.1);
%     spost_gsac_mod = NMMfit_filters(spost_gsac_mod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds),[],[],silent);
% 
%     %%
% %                 [gainLL,gain_pred_rate] = eval_sacgain_mod( temp_gain_mod, cur_Robs, all_Xmat(used_inds(cc_uinds),:), cur_Xsac);
% 
%     [spost_LL,~,spost_pred_rate,~,~,~,nullLL] = NMMmodel_eval(spost_gsac_mod, cur_Robs,tr_stim);
%     [sac_spost_info,sac_info,sac_arate,sac_spost_LL,sac_nullLL,sac_Nspks] = deal(nan(length(sac_bincents),1));
%     for ii = 1:length(sac_bincents)
%         tempr = find(cur_Xsac(:,ii) == 1);
%         sac_arate(ii) = mean(cur_Robs(tempr));
%         cur_avg_rate = mean(cur_Robs(tempr))*ones(size(tempr));
%         sac_nullLL(ii) = nansum(cur_Robs(tempr).*log2(cur_avg_rate) - cur_avg_rate);
%         sac_Nspks(ii) = sum(cur_Robs(tempr));
%         
%         %                         sac_LL(ii) = nansum(cur_Robs(tempr).*log2(gain_pred_rate(tempr)) - gain_pred_rate(tempr));
%         %                 sac_info(ii) = nanmean(gain_pred_rate(tempr).*log2(gain_pred_rate(tempr)/mean(gain_pred_rate(tempr))))/mean(gain_pred_rate(tempr));
% 
%         sac_spost_LL(ii) = nansum(cur_Robs(tempr).*log2(spost_pred_rate(tempr)) - spost_pred_rate(tempr));
%         sac_spost_info(ii) = nanmean(spost_pred_rate(tempr).*log2(spost_pred_rate(tempr)/mean(spost_pred_rate(tempr))))/mean(spost_pred_rate(tempr));
%     end
%     
%     gsac_spostLL_info = (sac_spost_LL - sac_nullLL)./sac_Nspks;
%     gsac_ov_LLinfo = (spost_LL - nullLL)/log(2);
%     
%     all_mu_modinfo(ss,:) = sac_spost_info;
%     all_mu_arate(ss,:) = sac_arate;
    
    %%
    flen_tax = (0:(flen-1))*dt + dt/2;
    ffr = reshape([all_mu_mods(ss).mods(1).filtK],flen,[]);
    tkern = sqrt(sum(ffr.^2,2));
    ffrn = bsxfun(@rdivide,ffr,tkern);

    f1 = figure;
    
    subplot(2,3,1)
    imagesc(flen_tax,poss_oris,flipud(ffr'));
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    axis tight
    xlabel('Time lag (s)');
    ylabel('Orientation');
    
    subplot(2,3,4)
    imagesc(flen_tax,poss_oris,flipud(ffrn'));
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    axis tight
    xlabel('Time lag (s)');
    ylabel('Orientation');

    subplot(2,3,2); hold on
    plot(sac_bincents*dt,ori_arate/dt,'k');
    for ii = 1:mlrange
        plot(sac_bincents(mmloc(1)-1+ii)*dt,ori_arate(mmloc(1)-1+ii)/dt,'o','color',cmap(ii,:),'linewidth',2);
    end
    xlabel('Time since sac (s)');
    ylabel('Firing rate (Hz)');
    xlim([-0.1 0.15])
    
    
    subplot(2,3,5); hold on
    plot(sac_bincents*dt,sta_ssinfo,'b','linewidth',2);
%     line(sac_bincents([1 end])*dt,[gsac_ov_LLinfo gsac_ov_LLinfo],'color','k')
    xlim([-0.05 0.2])
    xlabel('Time since sac (s)');
    ylabel('SS info (bits/spk)');
    xlim([-0.1 0.15])
        
    subplot(2,3,3); hold on
    imagesc(sac_bincents*dt,poss_oris,ori_avg');
    axis tight
    xlim([-0.1 0.15])
    ylabel('Orientation (deg)');
    xlabel('Time since sac (s)');
    title('Absolute rate');
    
    subplot(2,3,6); hold on
    imagesc(sac_bincents*dt,poss_oris,bsxfun(@rdivide,ori_avg,mean(ori_avg))');
    axis tight
    xlim([-0.1 0.15])
    ylabel('Orientation (deg)');
    xlabel('Time since sac (s)');
    title('Relative rate');
    
    
    fig_width = 11; rel_height = 0.6;
    
    figufy(f1);
    fname = [fig_dir sprintf('GsacE_neworimod3_E%d_MU%d.pdf',Expt_num,ss)];
    exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(f1);
        
% %     pause
% %     close all
    end
end

