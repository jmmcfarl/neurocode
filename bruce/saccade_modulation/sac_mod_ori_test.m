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

min_trial_dur = 0.75;

beg_buffer = 0.15;
end_buffer = 0.05;
trial_dur = 4;

n_probes = 24;

use_right_eye = false;

n_use_blocks = Inf;

sac_backlag = round(0.4/dt);
sac_forlag = round(0.6/dt);
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

big_sacs = big_sacs(~ismember(all_trialvec(used_inds(saccade_start_inds(big_sacs))),g4));


% dir_big_sacs1 = big_sacs(ismember(all_trialvec(used_inds(saccade_start_inds(big_sacs))),g1));
% dir_big_sacs2 = big_sacs(ismember(all_trialvec(used_inds(saccade_start_inds(big_sacs))),g2));
% dir_big_sacs3 = big_sacs(ismember(all_trialvec(used_inds(saccade_start_inds(big_sacs))),g3));
% dir_big_sacs4 = big_sacs(ismember(all_trialvec(used_inds(saccade_start_inds(big_sacs))),g4));
% big_sacs = dir_big_sacs(ismember(dir_big_sacs,outward_sacs));
% big_sacs = dir_big_sacs(ismember(dir_big_sacs,inward_sacs));



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

%     gsac_towards = used_gsacs(ismember(all_trialvec(used_inds(used_gsacs)),towards_trials));
%     gsac_away = used_gsacs(ismember(all_trialvec(used_inds(used_gsacs)),away_trials));

%%
% backlag = round(0.1/dt);
% forwardlag = round(4/dt);
% all_norm_binned_mua = bsxfun(@rdivide,all_binned_mua,nanmean(all_binned_mua));
% all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
% all_trial_start_inds(isnan(all_trial_start_inds)) = 1;
% [ov_gsac_trig_avg,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds,backlag,forwardlag);
% [ov_gsac_trig_g4,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g4),backlag,forwardlag);
% [ov_gsac_trig_g1,tlags] = get_event_trig_avg(all_norm_binned_mua,all_trial_start_inds(g1),backlag,forwardlag);
%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS

mod_signs = [1];
NL_types = {'lin'};
init_d2XT = [1];
init_Xtargs = [1];

base_lambda_d2XT = 10;
base_lambda_L1 = 1;

%create L2 mat
L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,flen*full_nPix);

init_reg_params = NMMcreate_reg_params('lambda_custom',init_d2XT);

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
        
        gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);

        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,all_Xmat(used_inds(cur_tr_inds),:));
        
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_custom',base_lambda_d2XT./var(gint)');
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);  
        cur_mod(ss) = gqm1;
    end
end

%%
mod_signs = [1];
NL_types = {'lin'};
init_d2XT = [1];
init_Xtargs = [1];

base_lambda_d2XT = 10;
base_lambda_L1 = 1;

%create L2 mat
L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,flen*full_nPix);

init_reg_params = NMMcreate_reg_params('lambda_custom',init_d2XT);
for ss = 1:length(su_probes)
    fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
    cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
        
        gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);

        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,all_Xmat(used_inds(cur_tr_inds),:));
        
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_custom',base_lambda_d2XT./var(gint)');
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(used_inds(cur_tr_inds),:),[],[],silent,[],L2_mat);  
        cur_sumod(ss) = gqm1;
    end
end

%%
for ss = 1:n_probes
   cur_filts = reshape(cur_mod(ss).mods(1).filtK,[flen length(poss_oris)]);
   tkern = std(cur_filts,[],2);
   [~,best_lag] = max(tkern);
   ori_tune(ss,:) = cur_filts(best_lag,:);
   ori_tkern(ss,:) = tkern;
end
%%
sac_stim_params(1) = NMMcreate_stim_params([1 1],dt);
sac_stim_params(2) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(3) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(4) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(5) = NMMcreate_stim_params([n_sac_bins 1],dt);

base_lambda_d2T = 100;
base_lambda_L2 = 2;

gain_lambda_d2T = 500;
gain_lambda_L2 = 5;

sac_BCs = [0 0 0];
sac_reg_params = NMMcreate_reg_params('lambda_d2T',base_lambda_d2T,'lambda_L2',base_lambda_L2,'boundary_conds',sac_BCs);
close all
for ss = 1:24;
    ss
    cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
    
    stim_filts = cur_mod(ss).mods(1).filtK;
    stim_out = all_Xmat(used_inds,:)*stim_filts;
    
    X{1} = stim_out(cur_tr_inds,:);
    Xsac_stim = bsxfun(@times,Xsac,stim_out);
    Xmsac_stim = bsxfun(@times,Xmsac,stim_out);
    
    X{2} = Xsac(cur_tr_inds,:);
    X{3} = Xmsac(cur_tr_inds,:);
    X{4} = Xsac_stim(cur_tr_inds,:);
    X{5} = Xmsac_stim(cur_tr_inds,:);
    
%         mod_signs = [1 1];
%     Xtargets = [2 3];
%     NL_types = {'lin','lin','lin','lin'};
     mod_signs = [1];
    Xtargets = [3];
    NL_types = {'lin','lin','lin','lin'};
    null_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    null_mod = NMMfit_filters(null_mod,Robs,X,[],[],silent);
    [nullLL, penLL, null_pred_rate] = NMMmodel_eval( null_mod, Robs, X);

%     mod_signs = [1 1 1];
%     Xtargets = [1 2 3];
    mod_signs = [1 1];
    Xtargets = [1 3];
    NL_types = {'lin','lin','lin','lin'};
    sac_full_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    sac_full_mod.mods(1).reg_params = NMMcreate_reg_params();
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[],silent);
    
    % filt_weights = [std(exc_stim_out) std(sup_stim_out)];
    filt_weights = abs(sac_full_mod.mods(1).filtK(1));
    filt_weights = [filt_weights; filt_weights];
    % filt_weights = [1; 1; 1; 1];
    cur_reg_params = NMMcreate_reg_params('lambda_d2T',gain_lambda_d2T*filt_weights.^2,...
        'lambda_L2',gain_lambda_L2*filt_weights.^2,'boundary_conds',[sac_BCs;sac_BCs]);
%     sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,4,zeros(n_sac_bins,1),cur_reg_params(1));
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},-1,5,zeros(n_sac_bins,1),cur_reg_params(2));
    
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[2 3],silent);
%     sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[2 3 4 5],silent);
    sac_filters = [sac_full_mod.mods([sac_full_mod.mods(:).Xtarget] > 1).filtK];
    
        [LL, penLL, pred_rate] = NMMmodel_eval( sac_full_mod, Robs, X);
    bsac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(big_sacs)));
    msac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(micro_sacs)));
    avg_sac_rate = get_event_trig_avg(pred_rate,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_rate = get_event_trig_avg(pred_rate,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_sac_R = get_event_trig_avg(Robs,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_R = get_event_trig_avg(Robs,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    
    all_sac_filters(ss,:,:) = sac_filters;
    all_sac_prate(ss,:,1) = avg_sac_rate;
    all_sac_prate(ss,:,2) = avg_msac_rate;
    all_sac_mrate(ss,:,1) = avg_sac_R;
    all_sac_mrate(ss,:,2) = avg_msac_R;
    
     ov_avg_rate = mean(Robs);
    null_LLvec = Robs.*log2(null_pred_rate)-null_pred_rate;
    sacmod_LLvec = Robs.*log2(pred_rate) - pred_rate;
    avg_sac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_sac_info_rate = jmm_smooth_1d_cor(avg_sac_info_rate,1)';
     avg_sac_rate = jmm_smooth_1d_cor(avg_sac_rate,1)';
    avg_sac_info = avg_sac_info_rate./avg_sac_rate;
     avg_msac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_info_rate = jmm_smooth_1d_cor(avg_msac_info_rate,1)';
     avg_msac_rate = jmm_smooth_1d_cor(avg_msac_rate,1)';
    avg_msac_info = avg_msac_info_rate./avg_msac_rate;
   
    ov_info = (LL-nullLL)/log(2);
    ov_info_rate = ov_info*ov_avg_rate;
    
    xx = linspace(prctile(stim_out,0.1),prctile(stim_out,99.9),100);
    pp = ksdensity(stim_out,xx);
    pp = pp/sum(pp);
    cur_sac_prate = nan(length(sac_bincents),length(xx));
     cur_msac_prate = nan(length(sac_bincents),length(xx));
   for ii = 1:length(sac_bincents)
%        cur_sac_prate(ii,:) = sac_filters(ii,1) + sac_filters(ii,3)*xx;
%        cur_msac_prate(ii,:) = sac_filters(ii,2) + sac_filters(ii,4)*xx;
       cur_msac_prate(ii,:) = sac_filters(ii,1) + sac_filters(ii,2)*xx;
    end
%     cur_sac_prate = sac_full_mod.spk_NL_params(1) + bsxfun(@plus,cur_sac_prate,xx*sac_full_mod.mods(1).filtK);
%     cur_sac_prate = log(1+exp(cur_sac_prate));
     cur_msac_prate = sac_full_mod.spk_NL_params(1) + bsxfun(@plus,cur_msac_prate,xx*sac_full_mod.mods(1).filtK);
    cur_msac_prate = log(1+exp(cur_msac_prate));
 
%     
%     cur_avg_rate = mean(Robs);
%     marg_gdist = sum(gsac_TB_dist,2);
%     marg_sdist = sum(gsac_TB_dist);
%     marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
%     marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
%     ov_info_gsac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
%     gsacdep_info = nan(n_sac_bins,1);
%     for tt = 1:n_sbins
%         gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
%     end
    

    
    xr = [-0.3 0.5];
    h = figure();
    subplot(3,2,1)
    plot(sac_bincents*dt,sac_filters(:,[1])); xlim(xr);
    ylabel('Sac offset');
    line(xr,[0 0],'color','k','linestyle','--');
    
    subplot(3,2,3)
    plot(sac_bincents*dt,1+sac_filters(:,[2]));xlim(xr);
     ylabel('Sac gain');
   line(xr,[1 1],'color','k','linestyle','--');
   
    subplot(3,2,5); hold on
    plot(sac_bincents*dt,avg_msac_rate);xlim(xr);
    line(xr,[ov_avg_rate ov_avg_rate],'color','k','linestyle','--');
    ylabel('Avg rate');
    
    subplot(3,2,[4]); hold on; 
    imagesc(sac_bincents*dt,xx,log(cur_msac_prate)');
    axis tight
    set(gca,'ydir','normal')
    prc_set = [10:10:90];
    for ii = 1:length(prc_set)
        line(xr,prctile(stim_out,prc_set([ii ii])),'color','w','linewidth',0.5);
    end
    xlim(xr);
    
    subplot(3,2,[2]); hold on
    cur_filts = reshape(stim_filts,[flen length(poss_oris)]);
    imagesc(poss_oris,(1:flen)*dt,cur_filts);
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    set(gca,'ydir','normal'); axis tight
    xlabel('Orientation');
    ylabel('Time lag (s)');

    subplot(3,2,6); hold on
    plot(sac_bincents*dt,avg_msac_info_rate/ov_info_rate);xlim(xr);
    plot(sac_bincents*dt,avg_msac_info/ov_info,'r');xlim(xr);
    line(xr,[1 1],'color','k','linestyle','--');
    ylabel('Relative information');
    
    figufy(h);
    pause
    close
%     fig_width = 8;
%     rel_height = 1.25;
%     outname = sprintf('/home/james/Desktop/lab_meeting_figs/ori_sacmod_MU%d',ss);
%     exportfig(h,outname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    
%     all_rel_info_cont(ss,:) = avg_sac_info_rate/ov_info_rate;
%     all_rel_info_cont_ps(ss,:) = avg_sac_info/ov_info;

%     subplot(3,1,1)
%     plot(sac_bincents*dt,sac_filters(:,[1 2]));
%     subplot(3,1,2)
%     plot(sac_bincents*dt,sac_filters(:,[3 4]));
%     subplot(3,1,3); hold on
%     plot(sac_bincents*dt,avg_sac_rate,sac_bincents*dt,avg_msac_rate,'r');
%     pause
%     clf
end
% avg_rates = nanmean(all_binned_mua(tr_inds,:));
% all_sac_mrate = bsxfun(@rdivide,all_sac_mrate,avg_rates');
%%
close all
for ss = 1:length(su_probes);
    %%
    cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
    
    if ~isempty(cur_tr_inds) & sum(Robs) > 0
    stim_filts = cur_sumod(ss).mods(1).filtK;
    stim_out = all_Xmat(used_inds,:)*stim_filts;
    
    X{1} = stim_out(cur_tr_inds,:);
    Xsac_stim = bsxfun(@times,Xsac,stim_out);
    Xmsac_stim = bsxfun(@times,Xmsac,stim_out);
    
    X{2} = Xsac(cur_tr_inds,:);
    X{3} = Xmsac(cur_tr_inds,:);
    X{4} = Xsac_stim(cur_tr_inds,:);
    X{5} = Xmsac_stim(cur_tr_inds,:);
    
    mod_signs = [1];
    Xtargets = [3];
    NL_types = {'lin','lin','lin','lin'};
    null_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    null_mod = NMMfit_filters(null_mod,Robs,X,[],[],silent);
    [nullLL, penLL, null_pred_rate] = NMMmodel_eval( null_mod, Robs, X);

    mod_signs = [1 1];
    Xtargets = [1 3];
    NL_types = {'lin','lin','lin','lin'};
    sac_full_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    sac_full_mod.mods(1).reg_params = NMMcreate_reg_params();
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[],silent);
    
    % filt_weights = [std(exc_stim_out) std(sup_stim_out)];
    filt_weights = abs(sac_full_mod.mods(1).filtK(1));
    filt_weights = [filt_weights; filt_weights];
    % filt_weights = [1; 1; 1; 1];
    cur_reg_params = NMMcreate_reg_params('lambda_d2T',gain_lambda_d2T*filt_weights.^2,...
        'lambda_L2',gain_lambda_L2*filt_weights.^2,'boundary_conds',[sac_BCs;sac_BCs]);
%     sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,4,zeros(n_sac_bins,1),cur_reg_params(1));
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},-1,5,zeros(n_sac_bins,1),cur_reg_params(2));
    
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[2 3],silent);
    sac_filters = [sac_full_mod.mods([sac_full_mod.mods(:).Xtarget] > 1).filtK];
    
    [LL, penLL, pred_rate] = NMMmodel_eval( sac_full_mod, Robs, X);
    bsac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(big_sacs)));
    msac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(micro_sacs)));
    avg_sac_rate = get_event_trig_avg(pred_rate,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_rate = get_event_trig_avg(pred_rate,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_sac_R = get_event_trig_avg(Robs,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_R = get_event_trig_avg(Robs,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    
    all_sac_filters(ss+n_probes,:,:) = sac_filters;
    all_sac_prate(ss+n_probes,:,1) = avg_sac_rate;
    all_sac_prate(ss+n_probes,:,2) = avg_msac_rate;
    all_sac_mrate(ss+n_probes,:,1) = avg_sac_R;
    all_sac_mrate(ss+n_probes,:,2) = avg_msac_R;

     ov_avg_rate = mean(Robs);
    null_LLvec = Robs.*log2(null_pred_rate)-null_pred_rate;
    sacmod_LLvec = Robs.*log2(pred_rate) - pred_rate;
    avg_sac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
     avg_sac_info_rate = jmm_smooth_1d_cor(avg_sac_info_rate,1)';
     avg_sac_rate = jmm_smooth_1d_cor(avg_sac_rate,1)';
   avg_sac_info = avg_sac_info_rate./avg_sac_rate;
     avg_msac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_info_rate = jmm_smooth_1d_cor(avg_msac_info_rate,1)';
     avg_msac_rate = jmm_smooth_1d_cor(avg_msac_rate,1)';
    avg_msac_info = avg_msac_info_rate./avg_msac_rate;
    
    ov_info = (LL-nullLL)/log(2);
    ov_info_rate = ov_info*ov_avg_rate;
    
    xx = linspace(prctile(stim_out,0.1),prctile(stim_out,99.9),100);
    pp = ksdensity(stim_out,xx);
    pp = pp/sum(pp);
    cur_sac_prate = nan(length(sac_bincents),length(xx));
     cur_msac_prate = nan(length(sac_bincents),length(xx));
   for ii = 1:length(sac_bincents)
%        cur_sac_prate(ii,:) = sac_filters(ii,1) + sac_filters(ii,3)*xx;
%        cur_msac_prate(ii,:) = sac_filters(ii,2) + sac_filters(ii,4)*xx;
       cur_msac_prate(ii,:) = sac_filters(ii,1) + sac_filters(ii,2)*xx;
    end
%     cur_sac_prate = sac_full_mod.spk_NL_params(1) + bsxfun(@plus,cur_sac_prate,xx*sac_full_mod.mods(1).filtK);
%     cur_sac_prate = log(1+exp(cur_sac_prate));
    cur_msac_prate = sac_full_mod.spk_NL_params(1) + bsxfun(@plus,cur_msac_prate,xx*sac_full_mod.mods(1).filtK);
    cur_msac_prate = log(1+exp(cur_msac_prate));
    
    xr = [-0.3 0.5];
%     h = figure();
%     subplot(3,2,1)
%     plot(sac_bincents*dt,sac_filters(:,[1])); xlim(xr); hold on
%     plot(sac_bincents*dt,sac_filters(:,[2]),'r'); 
%     ylabel('Sac offset');
%     line(xr,[0 0],'color','k','linestyle','--');
%     
%     subplot(3,2,3)
%     plot(sac_bincents*dt,1+sac_filters(:,[3]));xlim(xr); hold on
%     plot(sac_bincents*dt,1+sac_filters(:,[4]),'r');
%      ylabel('Sac gain');
%    line(xr,[1 1],'color','k','linestyle','--');
%    
%     subplot(3,2,5); hold on
%     plot(sac_bincents*dt,avg_sac_rate);xlim(xr); hold on
%     plot(sac_bincents*dt,avg_msac_rate,'r');
%     line(xr,[ov_avg_rate ov_avg_rate],'color','k','linestyle','--');
%     ylabel('Avg rate');
%     
%     subplot(3,2,[4]); hold on; 
%     imagesc(sac_bincents*dt,xx,log(cur_sac_prate)');
%     axis tight
%     set(gca,'ydir','normal')
%     prc_set = [10:10:90];
%     for ii = 1:length(prc_set)
%         line(xr,prctile(stim_out,prc_set([ii ii])),'color','w','linewidth',0.5);
%     end
%     xlim(xr);
%     
%     subplot(3,2,[2]); hold on
%     cur_filts = reshape(stim_filts,[flen length(poss_oris)]);
%     imagesc(poss_oris,(1:flen)*dt,cur_filts);
%     ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
%     set(gca,'ydir','normal'); axis tight
%     xlabel('Orientation');
%     ylabel('Time lag (s)');
% 
%     subplot(3,2,6); hold on
%     plot(sac_bincents*dt,avg_sac_info_rate/ov_info_rate);xlim(xr);
%     plot(sac_bincents*dt,avg_sac_info/ov_info,'--');
%     plot(sac_bincents*dt,avg_msac_info_rate/ov_info_rate,'r');
%     plot(sac_bincents*dt,avg_msac_info/ov_info,'r--');
%     line(xr,[1 1],'color','k','linestyle','--');
%     ylabel('Relative information');
%     
%     figufy(h);
    
%     fig_width = 8;
%     rel_height = 1.25;
%     outname = sprintf('/home/james/Desktop/lab_meeting_figs/ori_sacmod_SU%d',ss);
%     exportfig(h,outname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
%     
%     close(h);
%     pause
%     close all



    h = figure();  
    subplot(3,1,1); hold on
%     plot(sac_bincents*dt,avg_sac_rate/dt);xlim(xr); hold on
    plot(sac_bincents*dt,avg_msac_rate/dt,'r');
    legend('Sac','Micro-sac');
    line(xr,[ov_avg_rate ov_avg_rate]/dt,'color','k','linestyle','--');
    ylabel('Avg rate (Hz)');
    xlabel('Time since saccade onset (s)');
    
    subplot(3,1,2); hold on
%     plot(sac_bincents*dt,avg_sac_info_rate/ov_info_rate);xlim(xr);
    plot(sac_bincents*dt,avg_msac_info_rate/ov_info_rate,'r');
    line(xr,[1 1],'color','k','linestyle','--');
    ylabel('Relative information rate');
    
    subplot(3,1,3); hold on
%     plot(sac_bincents*dt,avg_sac_info/ov_info,'b');
    plot(sac_bincents*dt,avg_msac_info/ov_info,'r');
    line(xr,[1 1],'color','k','linestyle','--');
    ylabel('Relative information per spike');

    pause
    cloes all
%     figufy(h);
%     fig_width = 4;
%     rel_height = 2.5;
%     outname = sprintf('/home/james/Desktop/lab_meeting_figs/ori_sacmod_v2_SU%d',ss);
%     exportfig(h,outname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
%     close(h);

    end
end


%%
sac_stim_params(1) = NMMcreate_stim_params([1 1],dt);
sac_stim_params(2) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(3) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(4) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(5) = NMMcreate_stim_params([n_sac_bins 1],dt);

base_lambda_d2T = 100;
base_lambda_L2 = 2;

gain_lambda_d2T = 500;
gain_lambda_L2 = 5;

sac_BCs = [0 0 0];
sac_reg_params = NMMcreate_reg_params('lambda_d2T',base_lambda_d2T,'lambda_L2',base_lambda_L2,'boundary_conds',sac_BCs);
close all
for ss = 1:24;
    ss
    cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
    
    stim_filts = cur_mod(ss).mods(1).filtK;
    stim_out = all_Xmat(used_inds,:)*stim_filts;
    
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
    
    % filt_weights = [std(exc_stim_out) std(sup_stim_out)];
    filt_weights = abs(sac_full_mod.mods(1).filtK(1));
    filt_weights = [filt_weights; filt_weights];
    % filt_weights = [1; 1; 1; 1];
    cur_reg_params = NMMcreate_reg_params('lambda_d2T',gain_lambda_d2T*filt_weights.^2,...
        'lambda_L2',gain_lambda_L2*filt_weights.^2,'boundary_conds',[sac_BCs;sac_BCs]);
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,4,zeros(n_sac_bins,1),cur_reg_params(1));
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},-1,5,zeros(n_sac_bins,1),cur_reg_params(2));
    
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[2 3 4 5],silent);
    sac_filters = [sac_full_mod.mods([sac_full_mod.mods(:).Xtarget] > 1).filtK];
    
        [LL, penLL, pred_rate] = NMMmodel_eval( sac_full_mod, Robs, X);
    bsac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(big_sacs)));
    msac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(micro_sacs)));
    avg_sac_rate = get_event_trig_avg(pred_rate,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_rate = get_event_trig_avg(pred_rate,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_sac_R = get_event_trig_avg(Robs,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_R = get_event_trig_avg(Robs,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    
    all_sac_filters(ss,:,:) = sac_filters;
    all_sac_prate(ss,:,1) = avg_sac_rate;
    all_sac_prate(ss,:,2) = avg_msac_rate;
    all_sac_mrate(ss,:,1) = avg_sac_R;
    all_sac_mrate(ss,:,2) = avg_msac_R;
    
     ov_avg_rate = mean(Robs);
    null_LLvec = Robs.*log2(null_pred_rate)-null_pred_rate;
    sacmod_LLvec = Robs.*log2(pred_rate) - pred_rate;
    avg_sac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_sac_info_rate = jmm_smooth_1d_cor(avg_sac_info_rate,1)';
     avg_sac_rate = jmm_smooth_1d_cor(avg_sac_rate,1)';
    avg_sac_info = avg_sac_info_rate./avg_sac_rate;
     avg_msac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_info_rate = jmm_smooth_1d_cor(avg_msac_info_rate,1)';
     avg_msac_rate = jmm_smooth_1d_cor(avg_msac_rate,1)';
    avg_msac_info = avg_msac_info_rate./avg_msac_rate;
   
    ov_info = (LL-nullLL)/log(2);
    ov_info_rate = ov_info*ov_avg_rate;
    
    xx = linspace(prctile(stim_out,0.1),prctile(stim_out,99.9),100);
    pp = ksdensity(stim_out,xx);
    pp = pp/sum(pp);
    cur_sac_prate = nan(length(sac_bincents),length(xx));
     cur_msac_prate = nan(length(sac_bincents),length(xx));
   for ii = 1:length(sac_bincents)
       cur_sac_prate(ii,:) = sac_filters(ii,1) + sac_filters(ii,3)*xx;
       cur_msac_prate(ii,:) = sac_filters(ii,2) + sac_filters(ii,4)*xx;
    end
    cur_sac_prate = sac_full_mod.spk_NL_params(1) + bsxfun(@plus,cur_sac_prate,xx*sac_full_mod.mods(1).filtK);
    cur_sac_prate = log(1+exp(cur_sac_prate));
     cur_msac_prate = sac_full_mod.spk_NL_params(1) + bsxfun(@plus,cur_msac_prate,xx*sac_full_mod.mods(1).filtK);
    cur_msac_prate = log(1+exp(cur_msac_prate));
   
    xr = [-0.3 0.5];
    h = figure();
    subplot(3,2,1)
    plot(sac_bincents*dt,sac_filters(:,[2])); xlim(xr);
    ylabel('Sac offset');
    line(xr,[0 0],'color','k','linestyle','--');
    
    subplot(3,2,3)
    plot(sac_bincents*dt,1+sac_filters(:,[4]));xlim(xr);
     ylabel('Sac gain');
   line(xr,[1 1],'color','k','linestyle','--');
   
    subplot(3,2,5); hold on
    plot(sac_bincents*dt,avg_msac_rate);xlim(xr);
    line(xr,[ov_avg_rate ov_avg_rate],'color','k','linestyle','--');
    ylabel('Avg rate');
    
    subplot(3,2,[4]); hold on; 
    imagesc(sac_bincents*dt,xx,log(cur_msac_prate)');
    axis tight
    set(gca,'ydir','normal')
    prc_set = [10:10:90];
    for ii = 1:length(prc_set)
        line(xr,prctile(stim_out,prc_set([ii ii])),'color','w','linewidth',0.5);
    end
    xlim(xr);
    
    subplot(3,2,[2]); hold on
    cur_filts = reshape(stim_filts,[flen length(poss_oris)]);
    imagesc(poss_oris,(1:flen)*dt,cur_filts);
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    set(gca,'ydir','normal'); axis tight
    xlabel('Orientation');
    ylabel('Time lag (s)');

    subplot(3,2,6); hold on
    plot(sac_bincents*dt,avg_msac_info_rate/ov_info_rate);xlim(xr);
    plot(sac_bincents*dt,avg_msac_info/ov_info,'r');xlim(xr);
    line(xr,[1 1],'color','k','linestyle','--');
    ylabel('Relative information');
    
    figufy(h);
    
    fig_width = 8;
    rel_height = 1.25;
    outname = sprintf('/home/james/Desktop/lab_meeting_figs/ori_msacmod_MU%d',ss);
    exportfig(h,outname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    
end

%%
close all
for ss = 1:length(su_probes);
    %%
    cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
    
    if ~isempty(cur_tr_inds) & sum(Robs) > 0
    stim_filts = cur_sumod(ss).mods(1).filtK;
    stim_out = all_Xmat(used_inds,:)*stim_filts;
    
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
    
    % filt_weights = [std(exc_stim_out) std(sup_stim_out)];
    filt_weights = abs(sac_full_mod.mods(1).filtK(1));
    filt_weights = [filt_weights; filt_weights];
    % filt_weights = [1; 1; 1; 1];
    cur_reg_params = NMMcreate_reg_params('lambda_d2T',gain_lambda_d2T*filt_weights.^2,...
        'lambda_L2',gain_lambda_L2*filt_weights.^2,'boundary_conds',[sac_BCs;sac_BCs]);
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,4,zeros(n_sac_bins,1),cur_reg_params(1));
    sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},-1,5,zeros(n_sac_bins,1),cur_reg_params(2));
    
    sac_full_mod = NMMfit_filters(sac_full_mod,Robs,X,[],[2 3 4 5],silent);
    sac_filters = [sac_full_mod.mods([sac_full_mod.mods(:).Xtarget] > 1).filtK];
    
    [LL, penLL, pred_rate] = NMMmodel_eval( sac_full_mod, Robs, X);
    bsac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(big_sacs)));
    msac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(micro_sacs)));
    avg_sac_rate = get_event_trig_avg(pred_rate,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_rate = get_event_trig_avg(pred_rate,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_sac_R = get_event_trig_avg(Robs,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_R = get_event_trig_avg(Robs,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    
    all_sac_filters(ss+n_probes,:,:) = sac_filters;
    all_sac_prate(ss+n_probes,:,1) = avg_sac_rate;
    all_sac_prate(ss+n_probes,:,2) = avg_msac_rate;
    all_sac_mrate(ss+n_probes,:,1) = avg_sac_R;
    all_sac_mrate(ss+n_probes,:,2) = avg_msac_R;

     ov_avg_rate = mean(Robs);
    null_LLvec = Robs.*log2(null_pred_rate)-null_pred_rate;
    sacmod_LLvec = Robs.*log2(pred_rate) - pred_rate;
    avg_sac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
     avg_sac_info_rate = jmm_smooth_1d_cor(avg_sac_info_rate,1)';
     avg_sac_rate = jmm_smooth_1d_cor(avg_sac_rate,1)';
   avg_sac_info = avg_sac_info_rate./avg_sac_rate;
     avg_msac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
    avg_msac_info_rate = jmm_smooth_1d_cor(avg_msac_info_rate,1)';
     avg_msac_rate = jmm_smooth_1d_cor(avg_msac_rate,1)';
    avg_msac_info = avg_msac_info_rate./avg_msac_rate;
    
    ov_info = (LL-nullLL)/log(2);
    ov_info_rate = ov_info*ov_avg_rate;
    
    xx = linspace(prctile(stim_out,0.1),prctile(stim_out,99.9),100);
    pp = ksdensity(stim_out,xx);
    pp = pp/sum(pp);
    cur_sac_prate = nan(length(sac_bincents),length(xx));
     cur_msac_prate = nan(length(sac_bincents),length(xx));
   for ii = 1:length(sac_bincents)
       cur_sac_prate(ii,:) = sac_filters(ii,1) + sac_filters(ii,3)*xx;
       cur_msac_prate(ii,:) = sac_filters(ii,2) + sac_filters(ii,4)*xx;
    end
    cur_sac_prate = sac_full_mod.spk_NL_params(1) + bsxfun(@plus,cur_sac_prate,xx*sac_full_mod.mods(1).filtK);
    cur_sac_prate = log(1+exp(cur_sac_prate));
    cur_msac_prate = sac_full_mod.spk_NL_params(1) + bsxfun(@plus,cur_msac_prate,xx*sac_full_mod.mods(1).filtK);
    cur_msac_prate = log(1+exp(cur_msac_prate));
    
    xr = [-0.3 0.5];
    h = figure();
    subplot(3,2,1)
    plot(sac_bincents*dt,sac_filters(:,[2])); xlim(xr);
    ylabel('Sac offset');
    line(xr,[0 0],'color','k','linestyle','--');
    
    subplot(3,2,3)
    plot(sac_bincents*dt,1+sac_filters(:,[4]));xlim(xr);
     ylabel('Sac gain');
   line(xr,[1 1],'color','k','linestyle','--');
   
    subplot(3,2,5); hold on
    plot(sac_bincents*dt,avg_sac_rate);xlim(xr);
    line(xr,[ov_avg_rate ov_avg_rate],'color','k','linestyle','--');
    ylabel('Avg rate');
    
    subplot(3,2,[4]); hold on; 
    imagesc(sac_bincents*dt,xx,log(cur_sac_prate)');
    axis tight
    set(gca,'ydir','normal')
    prc_set = [10:10:90];
    for ii = 1:length(prc_set)
        line(xr,prctile(stim_out,prc_set([ii ii])),'color','w','linewidth',0.5);
    end
    xlim(xr);
    
    subplot(3,2,[2]); hold on
    cur_filts = reshape(stim_filts,[flen length(poss_oris)]);
    imagesc(poss_oris,(1:flen)*dt,cur_filts);
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    set(gca,'ydir','normal'); axis tight
    xlabel('Orientation');
    ylabel('Time lag (s)');

    subplot(3,2,6); hold on
    plot(sac_bincents*dt,avg_sac_info_rate/ov_info_rate);xlim(xr);
    plot(sac_bincents*dt,avg_sac_info/ov_info,'r');xlim(xr);
    line(xr,[1 1],'color','k','linestyle','--');
    ylabel('Relative information');
    
    figufy(h);
    
%     fig_width = 8;
%     rel_height = 1.25;
%     outname = sprintf('/home/james/Desktop/lab_meeting_figs/ori_msacmod_SU%d',ss);
%     exportfig(h,outname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    
%     close(h);
%     pause
%     close all
    end
end



%% 
[OO,LL] = meshgrid(poss_oris,1:flen);
use_lag = find(LL == 4);

X{1} = bsxfun(@times,Xsac,reshape(all_Xmat(used_inds,use_lag),[length(used_inds) 1 full_nPix]));
X{1} = reshape(X{1},length(used_inds),n_sac_bins*full_nPix);
X{2} = Xsac;

max_sac = max(Xsac,[],2);

cur_stim_params(1) = NMMcreate_stim_params([n_sac_bins full_nPix],dt);
cur_stim_params(2) = NMMcreate_stim_params([n_sac_bins 1],dt);

mod_signs = [1 1];
NL_types = {'lin','lin'};
init_d2XT = [1; 0];
init_Xtargs = [1 2];

base_lambda_d2XT = 5;
base_lambda_L1 = 0;

%create L2 mat
L2_params = create_L2_params([],[1 n_sac_bins*full_nPix],[n_sac_bins full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,n_sac_bins*full_nPix);

init_reg_params = NMMcreate_reg_params('lambda_custom',init_d2XT);

backlag = round(0.1/dt);
forwardlag = round(0.4/dt);

cd(anal_dir);

tot_nUnits = length(su_probes) + n_probes;
all_mod_SU = zeros(tot_nUnits,1);
all_mos_SUnum = zeros(tot_nUnits,1);
silent = 1;
for ss = 1:n_probes;
    fprintf('Computing base LLs for MU %d of %d\n',ss,n_probes);
    cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
    cur_tr_inds(max_sac(cur_tr_inds) < 1) = [];
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
        
        gqm1 = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds),[],[],silent,[],L2_mat);

        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds));
        
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_custom',base_lambda_d2XT./var(gint)');
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        gqm1 = NMMfit_filters(gqm1,Robs,get_Xcell_tInds(X,cur_tr_inds),[],[],silent,[],L2_mat);  
        sac_mod(ss) = gqm1;
    end
    
    cur_sac_set = find(ismember(cur_tr_inds,saccade_start_inds));
    sta(ss,:) = get_event_trig_avg(Robs,cur_sac_set,backlag,forwardlag);
    
end

%%
for ss = 1:length(su_probes);
    fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
    all_mod_SU(ss+n_probes) = su_probes(ss);
    all_mod_SUnum(ss+n_probes) = SU_numbers(ss);
    cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
    cur_xv_inds = xv_inds(~isnan(all_binned_sua(used_inds(xv_inds),ss)));
    tr_NT = length(cur_tr_inds);
    xv_NT = length(cur_xv_inds);
    Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
    Robsxv = all_binned_sua(used_inds(cur_xv_inds),ss);
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
        
        tr_X{1} = all_Xmat(used_inds(cur_tr_inds),use_kInds);
        tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
        xv_X{1} = all_Xmat(used_inds(cur_xv_inds),use_kInds);
        xv_X{2} = Xblock(used_inds(cur_xv_inds),:);
        if use_sac_kerns
            tr_X{3} = Xsac(cur_tr_inds,:);
            tr_X{4} = Xmsac(cur_tr_inds,:);
            xv_X{3} = Xsac(cur_xv_inds,:);
            xv_X{4} = Xmsac(cur_xv_inds,:);
            null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
        else
            null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
        end
        null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
        
        gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p);
        
        tr_X{1} = all_Xmat_us(used_inds(cur_tr_inds),use_kInds_up);
        xv_X{1} = all_Xmat_us(used_inds(cur_xv_inds),use_kInds_up);
        %spatial up-sampling of filter estimates
        base_filts = reshape([gqm1.mods(find(init_Xtargs == 1)).filtK],[flen use_nPix n_squared_filts+1]);
        if spatial_usfac == 2
            base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
            for ii = 1:use_nPix
                base_filts_up(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                base_filts_up(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
            end
        elseif spatial_usfac == 1
            base_filts_up = base_filts;
        else
            error('unsupported')
        end
        base_filts_up = reshape(base_filts_up,use_nPix_us*flen,	n_squared_filts+1);
        
        init_filts{end} = gqm1.mods(find(init_Xtargs==2)).filtK;
        for ii = 1:n_squared_filts+1
            init_filts{ii} = base_filts_up(:,ii);
        end
        gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
        gqm2.spk_NL_params(1) = gqm1.spk_NL_params(1);
        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
        gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
        gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        if use_sac_kerns
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
        end
        gqm2 = NMMfit_filters(gqm2,Robs,tr_X,[],[],silent);
        
        xvLL = NMMmodel_eval(gqm2,Robsxv,xv_X);
        null_xvLL(ss+n_probes) = NMMmodel_eval(null_mod,Robsxv,xv_X(2:end));
        all_mod_xvLLimp(ss+n_probes) = (xvLL-null_xvLL(ss+n_probes))/log(2);
        
        %now refit model using all (usable) data
        cur_tr_inds = sort([cur_tr_inds; cur_xv_inds]);
        Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
        tr_X{1} = all_Xmat_us(used_inds(cur_tr_inds),use_kInds_up);
        tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
        if use_sac_kerns
            tr_X{3} = Xsac(cur_tr_inds,:);
            tr_X{4} = Xmsac(cur_tr_inds,:);
        end
        null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
        all_nullmod(ss+n_probes) = null_mod;
        
        gqm2 = NMMfit_filters(gqm2,Robs,tr_X,[],[],silent);
        all_mod_fits(ss+n_probes) = gqm2;
        all_mod_fits_withspkNL(ss+n_probes) = NMMfit_logexp_spkNL(gqm2,Robs,tr_X);
        
        [LL,~,pred_rate] = NMMmodel_eval(all_mod_fits_withspkNL(ss+n_probes),Robs,tr_X);
        [null_LL(ss+n_probes),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
        [all_mod_R2(ss+n_probes),all_mod_dev(ss+n_probes),all_null_dev(ss+n_probes)] = pseudo_r2(Robs,pred_rate,null_prate);
        all_mod_LLimp(ss+n_probes) = (LL-null_LL(ss+n_probes))/log(2);
    else
        all_mod_LLimp(ss+n_probes) = nan;
        all_mod_xvLLimp(ss+n_probes) = nan;
        all_mod_R2(ss+n_probes) = nan;
        null_LL(ss+n_probes) = nan;
        null_xvLL(ss+n_probes) = nan;
    end
end

