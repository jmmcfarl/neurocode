clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

%run on 86 [0, 90];
Expt_num = 95;
Expt_name = sprintf('G%.3d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

% load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

anal_dir = ['~/Analysis/bruce/' Expt_name '/stim_mods/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = 'ori_tuning';

has_data = cellfun(@(x) length(x),Expts) > 0;
expt_names(has_data) = cellfun(@(x) x.Header.expname,Expts(has_data),'UniformOutput',false);

cur_block_set = find(strcmp('rls.orXme',expt_names));
n_blocks = length(cur_block_set);
%%

min_trial_dur = 0.4;
n_probes = 96;
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_blockvec = [];
all_trialvec = [];
all_trial_me = [];
all_trial_or = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];

all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
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

    trial_or = [Expts{cur_block}.Trials(:).or];
    trial_me = [Expts{cur_block}.Trials(:).me];
    
    all_trial_or = cat(1,all_trial_or,trial_or(use_trials)');
    all_trial_me = cat(1,all_trial_me,trial_me(use_trials)');
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    
        n_trials = length(use_trials);
    trial_cnt = trial_cnt + n_trials;
end

%% load overall su data
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

%% BIN SPIKES FOR MU AND SU
%for SU probes
fprintf('Using %d SUs\n',length(SU_numbers));
cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
SU_ID_mat = SU_ID_mat(cur_used_blocks,:);
[CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));

su_probes = nan(1,length(SU_numbers));
all_su_spk_times = cell(length(SU_numbers),1);
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
            
            all_su_spk_times{ss} = cat(1,all_su_spk_times{ss},cur_su_spk_times(:));
        end
        su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
    end
end

%%
beg_buff = 0.1;
end_buff = 0;
su_spk_cnts = nan(length(SU_numbers),trial_cnt);
mu_spk_cnts = nan(n_probes,trial_cnt);
for tt = 1:trial_cnt
    fprintf('Trial %d/%d\n',tt,trial_cnt);
   win_start = all_trial_start_times(tt) + beg_buff;
   win_end = all_trial_end_times(tt) - end_buff;
    for ss = 1:length(SU_numbers)
        su_spk_cnts(ss,tt) = sum(all_su_spk_times{ss} >= win_start & all_su_spk_times{ss} <= win_end);
    end
    for ss = 1:n_probes
        mu_spk_cnts(ss,tt) = sum(all_spk_times{ss} >= win_start & all_spk_times{ss} <= win_end);
    end
end

all_trial_durs = all_trial_end_times - all_trial_start_times;
su_trial_rate = bsxfun(@rdivide,su_spk_cnts,all_trial_durs');
mu_trial_rate = bsxfun(@rdivide,mu_spk_cnts,all_trial_durs');

su_avg_rate = mean(su_trial_rate,2);
mu_avg_rate = mean(mu_trial_rate,2);

su_rel_rate = bsxfun(@rdivide,su_trial_rate,su_avg_rate);
mu_rel_rate = bsxfun(@rdivide,mu_trial_rate,mu_avg_rate);

%%
poss_me = [-1 1];
poss_oris = unique(all_trial_or);
for mm = 1:length(poss_me)
    for oo = 1:length(poss_oris)
        cur_set = find(all_trial_me == poss_me(mm) & all_trial_or == poss_oris(oo));
        su_tuning(mm,oo,:) = mean(su_rel_rate(:,cur_set),2);
        mu_tuning(mm,oo,:) = mean(mu_rel_rate(:,cur_set),2);
    end
end
ov_su_tuning = squeeze(mean(su_tuning));
ov_mu_tuning = squeeze(mean(mu_tuning));

[~,pref_su_ori] = nanmax(ov_su_tuning);
[~,pref_mu_ori] = nanmax(ov_mu_tuning);
pref_su_ori = poss_oris(pref_su_ori);
pref_mu_ori = poss_oris(pref_mu_ori);

%%
for ss = 1:length(SU_numbers)
   SU_data(ss).avg_rate = su_avg_rate(ss);
   SU_data(ss).avg_ori_tuning = ov_su_tuning(:,ss);
   SU_data(ss).tot_ori_tuning = squeeze(su_tuning(:,:,ss));
   SU_data(ss).tot_spks = sum(su_spk_cnts(ss,:));
end
for ss = 1:n_probes
   MU_data(ss).avg_rate = mu_avg_rate(ss);
   MU_data(ss).avg_ori_tuning = ov_mu_tuning(:,ss);
    MU_data(ss).tot_ori_tuning = squeeze(mu_tuning(:,:,ss));
    MU_data(ss).tot_spks = sum(mu_spk_cnts(ss,:));
end

%%
% load ~/Data/bruce/general_array_data/array_pos_data.mat
% all_mu_oris = nan(10,10);
% for ii = 1:96
%     all_mu_oris(Y_pos(ii),X_pos(ii)) = pref_mu_ori(ii);
% end

%%
cd(anal_dir);
save(anal_name,'SU_data','MU_data');
