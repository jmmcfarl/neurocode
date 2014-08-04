clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

%run on 86 [0, 90];
Expt_num = 296;
Expt_name = sprintf('M%.3d',Expt_num);
% data_dir = ['~/Data/bruce/' Expt_name];
data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
cd(data_dir);

% load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

anal_dir = ['~/Analysis/bruce/' Expt_name '/stim_mods/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = 'ori_tuning';

has_data = cellfun(@(x) length(x),Expts) > 0;
expt_names(has_data) = cellfun(@(x) x.Header.expname,Expts(has_data),'UniformOutput',false);

% cur_block_set = find(strcmp('rls.orXme',expt_names));
cur_block_set = find(strcmp('rls.orRC',expt_names));
if Expt_num == 296
    cur_block_set(1) = [];
end
n_blocks = length(cur_block_set);

%%
dt = 0.01;
min_trial_dur = 0.4;
n_probes = 24;
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_blockvec = [];
all_trialvec = [];
all_trial_me = [];
all_trial_or = [];
all_stim_or = [];
all_stim_times = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_bin_edge_pts = [];
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
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        n_frames = length(cur_stim_times);
        if n_frames > 0
            cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt];
        end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_stim_or = Expts{cur_block}.Trials(use_trials(tt)).or;
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        
        if ~isempty(all_stim_times)
            if any(cur_stim_times+cur_toffset < all_stim_times(end))
                fprintf('Warn trial %d\n',tt);
            end
        end
        all_stim_times = [all_stim_times; cur_stim_times + cur_toffset];
        all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
        all_stim_or = [all_stim_or; cur_stim_or];
        all_tsince_start = [all_tsince_start; cur_tsince_start];
        all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    end
    trial_cnt = trial_cnt + n_trials;

    
        n_trials = length(use_trials);
    trial_cnt = trial_cnt + n_trials;
end

%% BIN SPIKES FOR MU AND SU
rec_type = 'UA';
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

%%
flen = 12;
poss_oris = unique(all_stim_or);
poss_oris(poss_oris < -90) = [];

full_nPix = length(poss_oris);
all_stim_mat = zeros(length(all_t_axis),full_nPix);
for ii = 1:full_nPix
    cur_set = find(all_stim_or == poss_oris(ii));
    all_stim_mat(cur_set,ii) = 1;
end

all_blank_mat = zeros(length(all_t_axis),1);
all_blank_mat(all_stim_or == -1009) = 1;

stim_params = NMMcreate_stim_params([flen full_nPix],dt);
uflen = 1;
ustim_params = NMMcreate_stim_params([uflen full_nPix],dt);
bstim_params = NMMcreate_stim_params([flen 1],dt);

all_Xmat = create_time_embedding(all_stim_mat,stim_params);
all_Xmat_blank = create_time_embedding(all_blank_mat,bstim_params);
%%
trial_dur = 2;
beg_buffer = 0.15;
end_buffer = 0;
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);

%%
X{1} = all_Xmat(used_inds,:);
X{2} = all_Xmat_blank(used_inds,:);

mod_stim_params(1) = stim_params;
mod_stim_params(2) = bstim_params;

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS

mod_signs = [1 1];
NL_types = {'lin','lin','lin'};
init_Xtargs = [1 2];

lambda_custom = 10;
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
    cur_tr_inds = find(~isnan(all_binned_mua(used_inds,ss)));
    tr_NT = length(cur_tr_inds);
    Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
    if ~isempty(cur_tr_inds) && sum(Robs) > 0
        
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

lambda_custom = 10;
lambda_L1 = 0;

%create L2 mat
L2_params = create_L2_params([],[1 flen*full_nPix],[flen full_nPix],2,3,[Inf -1]);
L2_mat = generate_L2_mat(L2_params,flen*full_nPix);


init_reg_params = NMMcreate_reg_params('lambda_custom',[lambda_custom; 0; 0]);
for ss = 1:length(su_probes)
    fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
    cur_tr_inds = find(~isnan(all_binned_sua(used_inds,ss)));
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


