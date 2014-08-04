clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

%run on 86 [0, 90];
Expt_num = 277;
Expt_name = sprintf('M%.3d',Expt_num);
if Expt_num > 280 && Expt_num < 289
    data_dir = ['/media/NTlab_data2/Data/bruce/' Expt_name];
elseif Expt_num >= 289
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end

cd(data_dir);

% load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

has_data = cellfun(@(x) length(x),Expts) > 0;
expt_names(has_data) = cellfun(@(x) x.Header.expname,Expts(has_data),'UniformOutput',false);

cur_block_set = find(strcmp('rls.orXFaRC',expt_names));
cur_block_set(cur_block_set == 26) = [];
% cur_block_set = find(strcmp('rls.orRC',expt_names));
n_blocks = length(cur_block_set);
%%
flen = 15;
full_nPix = 18;
% full_nPix = 9;
stim_fs = 100; %in Hz
dt = 0.01;
Fr = 1;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
uflen = 1;
ustim_params = NIMcreate_stim_params([uflen full_nPix],dt);

min_trial_dur = 0.75;

beg_buffer = 0.15;
end_buffer = 0.05;
trial_dur = 4;

n_probes = 24;

use_right_eye = false;

n_use_blocks = Inf;

sac_backlag = round(0.1/dt);
sac_forlag = round(0.25/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

%% load overall su data
% LOAD REFCLUSTERS
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
fname = [cluster_dir '/final_cluster.mat'];
if exist(fname,'file')
    load(fname);
    SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
    for ii = 1:length(SU_numbers)
        SU_tot_nblocks = sum(SU_ID_mat(:) == SU_numbers(ii));
    end
    fprintf('%d SUs Clustered\n',length(SU_numbers));
    
else
    disp('No Cluster data found.');
end

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
all_stim_mat = zeros(length(all_t_axis),full_nPix);
for ii = 1:full_nPix
    cur_set = find(all_stim_or == poss_oris(ii));
    all_stim_mat(cur_set,ii) = 1;
end

all_Xmat = create_time_embedding(all_stim_mat,stim_params);

[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
u_ks = find(Tinds == 6);
all_UXmat = all_Xmat(:,u_ks);

%% BIN SPIKES FOR MU AND SU
% LOAD REFCLUSTERS
fname = [cluster_dir '/final_cluster.mat'];
if exist(fname,'file')
    load(fname);
    SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
    for ii = 1:length(SU_numbers)
        SU_tot_nblocks = sum(SU_ID_mat(:) == SU_numbers(ii));
    end
    fprintf('%d SUs Clustered\n',length(SU_numbers));
    
else
    disp('No Cluster data found.');
end

%for SU probes
all_binned_sua = nan(length(all_t_axis),length(SU_numbers));
cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
SU_ID_mat = SU_ID_mat(cur_used_blocks,:);
[CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));

su_probes = nan(1,length(SU_numbers));
all_su_spk_times = cell(length(SU_numbers),1);
all_su_spk_inds = cell(length(SU_numbers),1);
for ss = 1:length(SU_numbers)
    used_clust_set = unique(CC(SU_ID_mat==SU_numbers(ss))); %set of clusters used to capture this SU
    SU_block_probes(ss,:) = nan(1,length(cur_block_set));
    cur_su_spk_times = [];
    cur_su_spk_inds = [];
    cur_blocks = [];
    for cc = 1:length(used_clust_set)
        cur_clust = used_clust_set(cc);
        cur_probe = SU_clust_data(cur_clust).probe_num;
        cur_clust_label = SU_clust_data(cur_clust).cluster_label;
        cur_blocks = [cur_blocks; find(SU_ID_mat(:,cur_clust) == SU_numbers(ss))];
        SU_block_probes(ss,cur_blocks) = cur_probe;
        
        all_su_inds = all_clust_ids{cur_probe} == cur_clust_label;
        cur_su_spk_times = all_spk_times{cur_probe}(all_su_inds);
        cur_su_spk_inds = all_spk_inds{cur_probe}(all_su_inds);
        spk_block_inds = round(interp1(all_t_axis,all_blockvec,cur_su_spk_times));
        cur_su_spk_times = cur_su_spk_times(ismember(spk_block_inds,cur_blocks));
        cur_su_spk_inds = cur_su_spk_inds(ismember(spk_block_inds,cur_blocks));
        
        all_su_spk_times{ss} = cat(1,all_su_spk_times{ss},cur_su_spk_times(:));
        all_su_spk_inds{ss} = cat(1,all_su_spk_inds{ss},cur_su_spk_inds(:));
    end
    if ~isempty(all_su_spk_times{ss})
        cur_suahist = histc(all_su_spk_times{ss},all_t_bin_edges);
        cur_suahist(all_bin_edge_pts) = [];
        cur_id_set = ismember(all_blockvec,cur_blocks);
        all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
        su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
    end
end

double_spike_buffer = 3; %number of samples (in either direction) to exclude double spikes from adjacent-probe SUs
all_binned_mua = nan(length(all_t_axis),n_probes);
for cc = 1:n_probes
    %this is the set of blocks where this probe had an SU, and the
    %correspodning SU numbers
    cur_set = find(SU_block_probes == cc);
    if ~isempty(cur_set)
        [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
    else
        cur_SS = [];
    end
    unique_su_nums = unique(cur_SS);
    cur_mua_inds = find(all_clust_ids{cc} >= 1);
    
    %remove spikes from isolated SUs on the same probe from the MU
    for ss = 1:length(unique_su_nums)
        cur_mua_inds(ismember(all_spk_inds{cc}(cur_mua_inds),all_su_spk_inds{unique_su_nums(ss)})) = [];
    end
    
    nearby_probes = [cc-1 cc+1]; nearby_probes(nearby_probes < 1 | nearby_probes > n_probes) = [];
    cur_set = find(ismember(SU_block_probes,nearby_probes));
    if ~isempty(cur_set)
        [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
    else
        cur_SS = [];
    end
    unique_su_nums = unique(cur_SS); %set of SU numbers picked up on adjacent probes
    if ~isempty(unique_su_nums)
        double_spikes = [];
        for ss = 1:length(unique_su_nums)
            cur_blocked_inds = bsxfun(@plus,all_su_spk_inds{unique_su_nums(ss)},-double_spike_buffer:double_spike_buffer);
            double_spikes = [double_spikes; find(ismember(all_spk_inds{cc}(cur_mua_inds),cur_blocked_inds))];
        end
        fprintf('Eliminating %d of %d double spikes in MUA\n',length(double_spikes),length(cur_mua_inds));
        cur_mua_inds(double_spikes) = [];
    end
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end



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
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,zeros(size(cur_block_set)),used_inds,false);

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);


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

%%
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

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
tr_inds = 1:length(used_inds);

%%
backlag = round(0.1/dt);
forwardlag = round(4/dt);
all_norm_binned_mua = bsxfun(@rdivide,all_binned_mua,nanmean(all_binned_mua));
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_start_inds(isnan(all_trial_start_inds)) = 1;
[ov_gsac_trig_avg,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds,backlag,forwardlag);
[ov_gsac_trig_g4,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g4),backlag,forwardlag);
[ov_gsac_trig_g2,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g2),backlag,forwardlag);
[ov_gsac_trig_g3,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g3),backlag,forwardlag);
[ov_gsac_trig_g1,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g1),backlag,forwardlag);


%%
backlag = round(0.1/dt);
forwardlag = round(0.3/dt);
[all_sac_avgs,lags] = get_event_trig_avg(all_norm_binned_mua(used_inds,:),saccade_start_inds(big_sacs),backlag,forwardlag);

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS

mod_signs = [1];
NL_types = {'lin'};
init_d2X = [1];
init_Xtargs = [1];

base_lambda_d2X = 10;
base_lambda_L1 = 1;

%create L2 mat
L2_params = create_L2_params([],[1 uflen*full_nPix],[uflen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,uflen*full_nPix);

init_reg_params = NMMcreate_reg_params('lambda_custom',init_d2X);

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
        
        gqm1 = NMMinitialize_model(ustim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,all_UXmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);

        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,all_UXmat(used_inds(cur_tr_inds),:));
        
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_custom',base_lambda_d2X./var(gint)');
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        gqm1 = NMMfit_filters(gqm1,Robs,all_UXmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);  
        all_mu_mods(ss) = gqm1;
    end
end

%%
mod_signs = [1];
NL_types = {'lin'};
init_d2X = [1];
init_Xtargs = [1];

base_lambda_d2X = 10;
base_lambda_L1 = 1;

%create L2 mat
L2_params = create_L2_params([],[1 uflen*full_nPix],[uflen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,uflen*full_nPix);

init_reg_params = NMMcreate_reg_params('lambda_custom',init_d2X);
for ss = 1:length(su_probes)
    fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
    cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
        
        gqm1 = NMMinitialize_model(ustim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,all_UXmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);

        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,all_UXmat(used_inds(cur_tr_inds),:));
        
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_custom',base_lambda_d2X./var(gint)');
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        gqm1 = NMMfit_filters(gqm1,Robs,all_UXmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);  
        all_su_mods(ss) = gqm1;
    end
end

%%
sac_stim_params(1) = NMMcreate_stim_params([1 1],dt);
sac_stim_params(2) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(3) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(4) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(5) = NMMcreate_stim_params([n_sac_bins 1],dt);

base_lambda_d2T = 1;
base_lambda_L2 = 0.5;

gain_lambda_d2T = 30;
gain_lambda_L2 = 1;

sac_BCs = [0 0 0];
sac_reg_params = NMMcreate_reg_params('lambda_d2T',base_lambda_d2T,'lambda_L2',base_lambda_L2,'boundary_conds',sac_BCs);
close all
for ss = 1:length(su_probes);
    %%
    cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
    
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
    stim_filts = all_su_mods(ss).mods(1).filtK;
    stim_out = all_UXmat(used_inds,:)*stim_filts;
    
    X{1} = stim_out(cur_tr_inds,:);
    Xsac_stim = bsxfun(@times,Xsac,stim_out);
    Xmsac_stim = bsxfun(@times,Xmsac,stim_out);
    
    X{2} = Xsac(cur_tr_inds,:);
    X{3} = Xmsac(cur_tr_inds,:);
    X{4} = Xsac_stim(cur_tr_inds,:);
    X{5} = Xmsac_stim(cur_tr_inds,:);
    
    mod_signs = [1 1];
    Xtargets = [2 3];
    NL_types = {'lin','lin','lin','lin'};
    null_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    null_mod = NMMfit_filters(null_mod,Robs,X,[],[],silent);
    [nullLL, penLL, null_pred_rate] = NMMmodel_eval( null_mod, Robs, X);
    
    mod_signs = [1 1 1];
    Xtargets = [1 2 3];
    NL_types = {'lin','lin','lin','lin'};
    sac_full_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    sac_full_mod.mods(1).reg_params = NMMcreate_reg_params();
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[],silent);
    
    %     filt_weights = abs(sac_full_mod.mods(1).filtK(1));
    %     filt_weights = [filt_weights; filt_weights];
    filt_weights = [1; 1];
    cur_reg_params = NMMcreate_reg_params('lambda_d2T',gain_lambda_d2T*filt_weights.^2,...
        'lambda_L2',gain_lambda_L2*filt_weights.^2,'boundary_conds',[sac_BCs;sac_BCs]);
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,4,zeros(n_sac_bins,1),cur_reg_params(1));
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,5,zeros(n_sac_bins,1),cur_reg_params(2));
    
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[2 3 4 5],silent);
    sac_filters = [sac_full_mod.mods([sac_full_mod.mods(:).Xtarget] > 1).filtK];
    
    subplot(2,1,1)
    plot(sac_bincents*dt,sac_filters(:,1:2));
    subplot(2,1,2)
    plot(sac_bincents*dt,sac_filters(:,3:4))
    pause
    clf
    %     

    end

end

%%
sac_stim_params(1) = NMMcreate_stim_params([1 1],dt);
sac_stim_params(2) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(3) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(4) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(5) = NMMcreate_stim_params([n_sac_bins 1],dt);

base_lambda_d2T = 1;
base_lambda_L2 = 0.5;

gain_lambda_d2T = 30;
gain_lambda_L2 = 1;

sac_BCs = [0 0 0];
sac_reg_params = NMMcreate_reg_params('lambda_d2T',base_lambda_d2T,'lambda_L2',base_lambda_L2,'boundary_conds',sac_BCs);
close all
for ss = 1:n_probes
    %%
    cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
    
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
    stim_filts = all_mu_mods(ss).mods(1).filtK;
    stim_out = all_UXmat(used_inds,:)*stim_filts;
    
    X{1} = stim_out(cur_tr_inds,:);
    Xsac_stim = bsxfun(@times,Xsac,stim_out);
    Xmsac_stim = bsxfun(@times,Xmsac,stim_out);
    
    X{2} = Xsac(cur_tr_inds,:);
    X{3} = Xmsac(cur_tr_inds,:);
    X{4} = Xsac_stim(cur_tr_inds,:);
    X{5} = Xmsac_stim(cur_tr_inds,:);
    
    mod_signs = [1 1];
    Xtargets = [2 3];
    NL_types = {'lin','lin','lin','lin'};
    null_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    null_mod = NMMfit_filters(null_mod,Robs,X,[],[],silent);
    [nullLL, penLL, null_pred_rate] = NMMmodel_eval( null_mod, Robs, X);

    mod_signs = [1 1 1];
    Xtargets = [1 2 3];
    NL_types = {'lin','lin','lin','lin'};
    sac_full_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    sac_full_mod.mods(1).reg_params = NMMcreate_reg_params();
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[],silent);
    
%     filt_weights = abs(sac_full_mod.mods(1).filtK(1));
%     filt_weights = [filt_weights; filt_weights];
filt_weights = [1; 1];
    cur_reg_params = NMMcreate_reg_params('lambda_d2T',gain_lambda_d2T*filt_weights.^2,...
        'lambda_L2',gain_lambda_L2*filt_weights.^2,'boundary_conds',[sac_BCs;sac_BCs]);
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,4,0.01*randn(n_sac_bins,1),cur_reg_params(1));
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,5,0.01*randn(n_sac_bins,1),cur_reg_params(2));
    
    optim_params.optTol = 1e-5;
    optim_params.progTol = 1e-8;
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[2 3 4 5],silent,optim_params);
    sac_filters = [sac_full_mod.mods([sac_full_mod.mods(:).Xtarget] > 1).filtK];
    
   

%     [LL, penLL, pred_rate, G,gint,fgint] = NMMmodel_eval( cur_mod(ss),Robs, all_Xmat(cur_tr_inds,:));
%     g_tot = G - cur_mod(ss).spk_NL_params(1);
%     [post_best_LL] = NMMmodel_eval(sac_full_mod, Robs, X);
%     
%     for ii = 1:length(sac_bincents)
%         fprintf('Testing post-scramble at gsac lag %d of %d\n',ii,length(sac_bincents));
%         to_shuffle = find(Xsac(cur_tr_inds,ii) == 1);
%         shuf_X = all_Xmat(cur_tr_inds,:);
%         shuf_X(to_shuffle,:) = shuf_X(randi(length(cur_tr_inds),length(to_shuffle),1),:);
%         
%         [LL, penLL, pred_rate, G,gint,fgint] = NMMmodel_eval( cur_mod(ss),Robs, shuf_X);
%         g_tot = G - cur_mod(ss).spk_NL_params(1);
%         Xsac_tot = bsxfun(@times,Xsac(cur_tr_inds,:),g_tot);
%         Xmsac_tot = bsxfun(@times,Xmsac(cur_tr_inds,:),g_tot);
%         X{1} = g_tot;
%         X{4} = Xsac_tot;
%         X{5} = Xmsac_tot;
%         
%         [post_shuf_LL(ii)] = NMMmodel_eval( sac_full_mod,Robs, X);
%         
%         n_frames(ii) = length(to_shuffle);
%     end
%     post_gsac_stiminfo(cc,:) = (post_best_LL-post_shuf_LL)./n_frames*length(cur_tr_inds);
% 
    
    ss
    subplot(2,1,1)
    plot(sac_filters(:,1:2));
    subplot(2,1,2)
    plot(sac_filters(:,3:4))
    pause
    clf
    %     

    end

end

%%
time_to_gsac = nan(length(all_t_axis),1);
time_to_msac = nan(length(all_t_axis),1);
max_dist = round(0.4/dt);
rel_tax = (-max_dist:max_dist)';
for ii = 1:length(big_sacs)
    cur_inds = (used_inds(saccade_start_inds(big_sacs(ii))) - max_dist):(used_inds(saccade_start_inds(big_sacs(ii))) + max_dist);
    cur_uset = find(cur_inds > 1 & cur_inds < length(all_t_axis));
    cur_uset(abs(time_to_gsac(cur_inds(cur_uset))) < abs(rel_tax(cur_uset))) = [];
    time_to_gsac(cur_inds(cur_uset)) = rel_tax(cur_uset);
end
for ii = 1:length(micro_sacs)
    cur_inds = (used_inds(saccade_start_inds(micro_sacs(ii))) - max_dist):(used_inds(saccade_start_inds(micro_sacs(ii))) + max_dist);
    cur_uset = find(cur_inds > 1 & cur_inds < length(all_t_axis));
    cur_uset(abs(time_to_gsac(cur_inds(cur_uset))) < abs(rel_tax(cur_uset))) = [];
    time_to_msac(cur_inds(cur_uset)) = rel_tax(cur_uset);
end
time_to_gsac = time_to_gsac*dt;
time_to_msac = time_to_msac*dt;

%%
back_t = 0.15;
for_t = 0.3;
n_Gbins = 25;
Xtick = -(back_t-dt/2):(dt):(for_t+dt/2);
n_sbins = length(Xtick);
addpath('~/James_scripts/TentBasis2D/');
%%
silent = 1;
close all
clear all_sac_filters
for ss = 1:n_probes
    ss
%     cur_Robs = all_binned_mua(used_inds,ss);
    cur_Robs = all_binned_sua(used_inds,ss);
    cur_tr_inds = tr_inds(~isnan(cur_Robs(tr_inds)));
    
    if ~isempty(cur_tr_inds)
%         cur_mod = all_mu_mods(ss);
        cur_mod = all_su_mods(ss);
        Xtargs = [cur_mod.mods(:).Xtarget];
        stim_NL_types = {cur_mod.mods(:).NLtype};
        stim_mod_signs = [cur_mod.mods(Xtargs==1).sign];
        stim_filters = [cur_mod.mods(Xtargs == 1).filtK];
        stim_outs = all_UXmat*stim_filters;
        G = sum(bsxfun(@times,stim_outs,stim_mod_signs),2);
        
        TB_stim = [time_to_gsac(used_inds(cur_tr_inds)) G(used_inds(cur_tr_inds))];
%         TB_stim = [time_to_msac(used_inds(cur_tr_inds)) G(used_inds(cur_tr_inds))];
        Ytick = linspace(my_prctile(TB_stim(:,2),0.1),my_prctile(TB_stim(:,2),100-1),n_Gbins);
        TB = TentBasis2D(Xtick, Ytick);
        
        used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
            TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
        [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
        L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
        TB_fitmod = regGLM_fit(TB_Xmat,cur_Robs(used_data),L2_params,50,[],[],1);
        [LL, penLL, pred_rate, G] = regGLM_eval(TB_fitmod,cur_Robs(used_data),TB_Xmat);
        TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
        bin_areas = TB.GetBinAreas();
        gsac_TB_dist = TB_counts./bin_areas;
        gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
        gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
        
        %INFO CALS
        cur_avg_rate = mean(cur_Robs(used_data));
        marg_gdist = sum(gsac_TB_dist,2);
        marg_sdist = sum(gsac_TB_dist);
        marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
        marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
        ov_info_gsac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        ov_info_gsac_ps = ov_info_gsac/cur_avg_rate;
        
        gsacdep_info = nan(n_sbins,1);
        for tt = 1:n_sbins
            gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
        end
        gsacdep_info_ps = gsacdep_info./marg_gsacrate';
        
        figure
        imagesc(Xtick,Ytick,gsac_TB_rate);
        set(gca,'ydir','normal');
        
        figure
        subplot(2,2,1);
        plot(Xtick,marg_gsacrate)
        line(Xtick([1 end]),[cur_avg_rate cur_avg_rate],'color','k');
        axis tight
        subplot(2,2,2);
        hold on
        plot(Xtick,gsacdep_info);
        line(Xtick([1 end]),[ov_info_gsac ov_info_gsac],'color','k');
        axis tight
        subplot(2,2,3)
        plot(Xtick,gsacdep_info_ps,'r');
        line(Xtick([1 end]),[ov_info_gsac_ps ov_info_gsac_ps],'color','g');
        axis tight
        pause
        close all
    end
end
