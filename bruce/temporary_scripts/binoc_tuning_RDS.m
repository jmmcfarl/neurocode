clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_name = 'G086';
Expt_num = str2num(Expt_name(2:end));
if Expt_num > 280 && Expt_num < 289
    data_dir = ['/media/NTlab_data2/Data/bruce/' Expt_name];
elseif Expt_num >= 289
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    system(['mkdir ' data_dir]);
end
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

n_probes = 96;
min_trial_dur = 0.3;
dt = 0.01;
%%
temp = cellfun(@(x) x.Header.expname,Expts,'uniformoutput',false);

% load overall su data
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

%%
cur_block_set = find(strcmp(temp,'rds.dxXce'));

n_blocks = length(cur_block_set);


%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_me = [];
all_trial_ce = [];
all_trial_dx = [];
all_trial_st = [];
all_trial_o2 = [];
all_trial_cb = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
cur_spkind_offset = 0;

trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
for ee = 1:n_blocks;
    fprintf('Block %d of %d\n',ee,n_blocks);
    cur_block = cur_block_set(ee);
    
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end

    fname = [data_dir sprintf('/Expt%dClusterTimes.mat',cur_block)];
    load(fname);
    for cc = 1:n_probes
       all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times'); 
    end
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
%     if length(un_ids) < length(trial_ids)
%         rpt_trials = true;
%         fprintf('Warning, repeat trial inds detected!\n');
%         use_trials = [];
%     else
        use_trials = find(trial_durs >= min_trial_dur);
%     end
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    trial_ce = [Expts{cur_block}.Trials(:).ce];
    all_trial_ce = cat(1,all_trial_ce,trial_ce(use_trials)');
    
    trial_me = [Expts{cur_block}.Trials(:).me];
    all_trial_me = cat(1,all_trial_me,trial_me(use_trials)');

    trial_dx= [Expts{cur_block}.Trials(:).dx];
    all_trial_dx = cat(1,all_trial_dx,trial_dx(use_trials)');

        trial_cb= [Expts{cur_block}.Trials(:).cb];
    all_trial_cb = cat(1,all_trial_cb,trial_cb(use_trials)');

    trial_o2= [Expts{cur_block}.Trials(:).o2];
    all_trial_o2 = cat(1,all_trial_o2,trial_o2(use_trials)');

    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_t_edges = [trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt))];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        all_t_axis = [all_t_axis; cur_t_axis'];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges'];
        all_tsince_start = [all_tsince_start; cur_tsince_start'];
        all_blockvec = [all_blockvec; ones(length(cur_t_axis),1)*ee];
        all_trialvec = [all_trialvec; ones(length(cur_t_axis),1)*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    end
    trial_cnt = trial_cnt + n_trials;
end

%%
all_trial_ce(all_trial_ce == 1 & all_trial_o2 ~= 3 & all_trial_dx == 0) = nan;

%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = n_probes;
% if strcmp(rec_type,'LP')
%     clust_params.exclude_adjacent = true;
% else
    clust_params.exclude_adjacent = false;
% end
[all_binned_mua,all_binned_sua,Clust_data] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;
SU_isodist = Clust_data.SU_isodists;
SU_Lratio = Clust_data.SU_Lratios;

%% DEFINE DATA USED FOR ANALYSIS
beg_buffer = 0.1;
end_buffer = 0.05;
trial_dur = 0.4;
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);


%%
poss_ces = unique(all_trial_ce);
poss_dxs = unique(all_trial_dx);
ov_avg_rates = mean(all_binned_mua(used_inds,:))/dt;
ov_avg_SUrates = mean(all_binned_sua(used_inds,:))/dt;
%%
clear avg_ce_rates avg_ce_SUrates
for ii = 1:length(poss_ces)
    cur_trials = find(all_trial_ce == poss_ces(ii));
    cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_ce_rates(ii,:) = mean(all_binned_mua(cur_inds,:))/dt;
    avg_ce_SUrates(ii,:) = mean(all_binned_sua(cur_inds,:))/dt;
end
avg_ce_rel_rates = bsxfun(@rdivide,avg_ce_rates,ov_avg_rates);
avg_ce_rel_SUrates = bsxfun(@rdivide,avg_ce_SUrates,ov_avg_SUrates);

%%
clear avg_dx_rates* avg_dx_SUrates*
for ii = 1:length(poss_dxs)
    cur_trials = find(all_trial_dx == poss_dxs(ii) & all_trial_ce == 1 & all_trial_me == 0);
%     cur_trials = find(all_trial_dx == poss_dxs(ii) & all_trial_ce == 1 & all_trial_me == 0);
%     cur_trials = find(all_trial_dx == poss_dxs(ii) & all_trial_ce == 1);
    cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_dx_rates(ii,:) = mean(all_binned_mua(cur_inds,:))/dt;
    avg_dx_SUrates(ii,:) = mean(all_binned_sua(cur_inds,:))/dt;
    std_dx_rates(ii,:) = std(all_binned_mua(cur_inds,:))/dt/sqrt(length(cur_trials));
    std_dx_SUrates(ii,:) = std(all_binned_sua(cur_inds,:))/dt/sqrt(length(cur_trials));
    n_examps(ii) = length(cur_trials);
    
    cur_trials = find(all_trial_dx == poss_dxs(ii) & all_trial_ce == -1 & all_trial_me == 0);
%     cur_trials = find(all_trial_dx == poss_dxs(ii) & all_trial_ce == -1);
    cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_Ndx_rates(ii,:) = mean(all_binned_mua(cur_inds,:))/dt;
    avg_Ndx_SUrates(ii,:) = mean(all_binned_sua(cur_inds,:))/dt;
end
avg_dx_rel_rates = bsxfun(@rdivide,avg_dx_rates,ov_avg_rates);
avg_dx_rel_SUrates = bsxfun(@rdivide,avg_dx_SUrates,ov_avg_SUrates);
avg_Ndx_rel_rates = bsxfun(@rdivide,avg_Ndx_rates,ov_avg_rates);
avg_Ndx_rel_SUrates = bsxfun(@rdivide,avg_Ndx_SUrates,ov_avg_SUrates);

%%
[~,best_ori] = max(avg_ori_rates);

for cc = 1:n_probes
    for ii = 1:length(poss_sfs)
    cur_trials = find(all_trial_or == poss_ors(best_ori(cc)) & all_trial_sf == poss_sfs(ii));
     cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_sf_rates_bestori(ii,cc) = mean(all_binned_mua(cur_inds,cc))/dt;
   end
end
avg_sf_rel_rates_bestori = bsxfun(@rdivide,avg_sf_rates_bestori,ov_avg_rates);