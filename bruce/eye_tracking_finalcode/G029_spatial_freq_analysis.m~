clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 29;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir');
    system(['mkdir ' data_dir]);
end
cd(data_dir);

load(sprintf('%sExpts.mat',Expt_name)); %load in Expts struct

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

n_probes = 96;
min_trial_dur = 0.3;
dt = 0.01;

%%
temp = cellfun(@(x) x.Header.expname,Expts,'uniformoutput',false);
cur_block_set = find(strcmp(temp,'grating.orXsf'));
n_blocks = length(cur_block_set);

%% LOAD REFCLUSTERS
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

all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_or = [];
all_trial_tf = [];
all_trial_sf = [];
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
    fprintf('Block %d of %d\n',ee,n_blocks);
    cur_block = cur_block_set(ee);
    
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end

%     fname = [data_dir sprintf('/Expt%dClusterTimes.mat',cur_block)];
%     load(fname);
%     for cc = 1:n_probes
%        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times'); 
%     end
    
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
    
    trial_or = [Expts{cur_block}.Trials(:).or];
%     trial_or = trial_or(id_inds);
    all_trial_or = cat(1,all_trial_or,trial_or(use_trials)');
    
    trial_sf= [Expts{cur_block}.Trials(:).sf];
%     trial_sf = trial_sf(id_inds);
    all_trial_sf = cat(1,all_trial_sf,trial_sf(use_trials)');
    
%     trial_tf= [Expts{cur_block}.Trials(:).tf];
%     trial_tf = trial_tf(id_inds);
%     all_trial_tf = cat(1,all_trial_tf,trial_tf(use_trials)');
        
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


%% BIN SPIKES FOR MU AND SU
%for SU probes
fprintf('Using %d SUs\n',length(SU_numbers));
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
    if ~isempty(used_clust_set)
        for cc = 1:length(used_clust_set)
            cur_clust = used_clust_set(cc);
            cur_probe = SU_clust_data(cur_clust).probe_num;
            cur_clust_label = SU_clust_data(cur_clust).cluster_label;
            cur_blocks = [cur_blocks find(SU_ID_mat(:,cur_clust) == SU_numbers(ss))];
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
        
        cur_suahist = histc(all_su_spk_times{ss},all_t_bin_edges);
        cur_suahist(all_bin_edge_pts) = [];
        cur_id_set = ismember(all_blockvec,cur_blocks);
        all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
        su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
    end
end

%for only-MU probes
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
    
    %remove spikes from isolated SUs from the MU
    for ss = 1:length(unique_su_nums)
       cur_mua_inds(ismember(all_spk_inds{cc}(cur_mua_inds),all_su_spk_inds{unique_su_nums(ss)})) = []; 
    end
    
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end

%% DEFINE DATA USED FOR ANALYSIS
beg_buffer = 0.1;
end_buffer = 0.05;
trial_dur = 0.4;
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);


%%
poss_sfs = unique(all_trial_sf);
poss_ors = unique(all_trial_or);
ov_avg_rates = nanmean(all_binned_mua(used_inds,:))/dt;
ov_avg_rates_SU = nanmean(all_binned_sua(used_inds,:))/dt;
bad_su = find(all(isnan(all_binned_sua(used_inds,:))));
%% FOR MUA
for ii = 1:length(poss_ors)
    cur_trials = find(all_trial_or == poss_ors(ii));
    cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_ori_rates(ii,:) = mean(all_binned_mua(cur_inds,:))/dt;
end
avg_ori_rel_rates = bsxfun(@rdivide,avg_ori_rates,ov_avg_rates);

%
for ii = 1:length(poss_sfs)
    cur_trials = find(all_trial_sf == poss_sfs(ii));
    cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_sf_rates(ii,:) = mean(all_binned_mua(cur_inds,:))/dt;
end
avg_sf_rel_rates = bsxfun(@rdivide,avg_sf_rates,ov_avg_rates);

%
[avg_rate_best_ori,best_ori] = max(avg_ori_rates);

for cc = 1:n_probes
    for ii = 1:length(poss_sfs)
    cur_trials = find(all_trial_or == poss_ors(best_ori(cc)) & all_trial_sf == poss_sfs(ii));
     cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_sf_rates_bestori(ii,cc) = mean(all_binned_mua(cur_inds,cc))/dt;
   end
end
avg_sf_rel_rates_bestori = bsxfun(@rdivide,avg_sf_rates_bestori,avg_rate_best_ori);
[~,best_sf_loc] = max(avg_sf_rel_rates_bestori);
best_sf = poss_sfs(best_sf_loc);
%% FOR SUA
for ii = 1:length(poss_ors)
    cur_trials = find(all_trial_or == poss_ors(ii));
    cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_ori_rates_SU(ii,:) = nanmean(all_binned_sua(cur_inds,:))/dt;
end
avg_ori_rel_rates_SU = bsxfun(@rdivide,avg_ori_rates_SU,ov_avg_rates_SU);

%
for ii = 1:length(poss_sfs)
    cur_trials = find(all_trial_sf == poss_sfs(ii));
    cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_sf_rates_SU(ii,:) = nanmean(all_binned_sua(cur_inds,:))/dt;
end
avg_sf_rel_rates_SU = bsxfun(@rdivide,avg_sf_rates_SU,ov_avg_rates_SU);

%
[avg_rate_best_ori_SU,best_ori_SU] = max(avg_ori_rates_SU);

for cc = 1:length(SU_numbers)
    for ii = 1:length(poss_sfs)
    cur_trials = find(all_trial_or == poss_ors(best_ori_SU(cc)) & all_trial_sf == poss_sfs(ii));
     cur_inds = used_inds(ismember(all_trialvec(used_inds),cur_trials));
    avg_sf_rates_bestori_SU(ii,cc) = nanmean(all_binned_sua(cur_inds,cc))/dt;
   end
end
avg_sf_rel_rates_bestori_SU = bsxfun(@rdivide,avg_sf_rates_bestori_SU,avg_rate_best_ori_SU);
[~,best_sf_loc_SU] = max(avg_sf_rel_rates_bestori_SU);
best_sf_SU = poss_sfs(best_sf_loc_SU);
best_sf_SU(bad_su) = nan;