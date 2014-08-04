clear all
close all

min_trial_dur = 2;
trial_dur = 4;
use_nPix = 20;
dt = 0.01;
flen = 12;
xv_frac = 0;
n_probes = 24;

data_dir = '~/James_scripts/bruce/online_processing/';
cd(data_dir);

%%
monName = 'lem';
exp_name = 'M266';
block_nums = [3];
Expts = {};
for bb = 1:length(block_nums)
    dat_name = [pwd sprintf('/%s%s.%d.mat',monName,exp_name,block_nums(bb))];
    [a,cur_Expt] = APlaySpkFile(dat_name,'nospikes','noerrs');
    Expts = {Expts{:} cur_Expt{:}};
end
load('./stim_data.mat');
load('./expt_data.mat');

%%
fprintf('Computing prep data\n');
trial_cnt = 0;
cur_toffset = 0;

all_stim_times = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
trial_toffset = zeros(length(block_nums),1);
for ee = 1:length(Expts);
    fprintf('Block %d of %d\n',ee,length(Expts));
    cur_block_num = block_nums(ee);
    
    Fr = Expts{ee}.Stimvals.Fr;
    
    % LOAD SPIKING DATA
    spike_data_name = [data_dir sprintf('/Block%d_Clusters.mat',cur_block_num)];
    if exist(spike_data_name,'file')
        load(spike_data_name);
        
        for cc = 1:n_probes
            all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
            all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
            all_clust_iso_dists(cc,ee,:) = Clusters{cc}.iso_dists;
            all_clust_Lratios(cc,ee,:) = Clusters{cc}.Lratios;
        end
        
    else
        fprintf('No spiking data found for block %d\n',cur_block_num);
    end
    clear Clusters
    
    trial_start_times = [Expts{ee}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{ee}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{ee}.Trials(:).dur]/1e4;
    trial_ids = [Expts{ee}.Trials(:).id];
    
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
    
    trial_Se = [Expts{ee}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    
    fname = sprintf('Expt%d_stim.mat',cur_block_num);
    load(fname);
    buffer_pix = floor((expt_npix(ee) - use_nPix)/2);
    cur_use_pix = (1:use_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{ee}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if length(cur_stim_times) == 1
            cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
            cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            
            all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_stim_times = [all_stim_times; cur_stim_times' + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
            
        end
    end
    trial_cnt = trial_cnt + n_trials;
    trial_toffset(ee) = all_t_bin_edges(end);
    cur_toffset = trial_toffset(ee);
end

%% BIN SPIKES FOR MU AND SU
all_binned_spikes = nan(length(all_t_axis),n_probes);
su_binned_spikes = [];
for cc = 1:n_probes
    cur_mua = find(all_clust_ids{cc} >= 1);
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_spikes(:,cc) = cur_spkhist;
end

%% CREATE STIMULUS MATRIX
stim_params = NIMcreate_stim_params([flen use_nPix],dt);
all_Xmat = create_time_embedding(all_stim_mat,stim_params);

% DEFINE DATA USED FOR ANALYSIS
beg_buffer = 0.2;
end_buffer = 0.1;
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);

% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

n_xv_trials = round(xv_frac*nuse_trials);
xv_trials = randperm(nuse_trials);
xv_trials(n_xv_trials+1:end) = [];
xv_trials = use_trials(xv_trials);
tr_trials = setdiff(use_trials,xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

tr_X{1} = abs(all_Xmat(used_inds(tr_inds),:));
xv_X{1} = abs(all_Xmat(used_inds(xv_inds),:));
tr_X{2} = all_Xmat(used_inds(tr_inds),:);
xv_X{2} = all_Xmat(used_inds(xv_inds),:);

%%
init_stim_params(1) = NMMcreate_stim_params([flen use_nPix],dt);
init_stim_params(2) = NMMcreate_stim_params([flen use_nPix],dt);
silent = 1;
base_lambda_d2XT = 5;
base_lambda_L1 = 0;
init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-9;

n_stim_filts = 2;
mod_signs = ones(1,n_stim_filts);
% NL_types = [{'lin'} repmat({'quad'},1,n_stim_filts-1)];
NL_types = [{'lin' 'lin'}];
init_d2XT = [ones(n_stim_filts,1)];
init_L2 = [zeros(n_stim_filts,1);];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT);
init_Xtargs = [1; 2];

for ss = 1:n_probes;
    fprintf('Fitting model for MU %d of %d\n',ss,n_probes);
    Robs = all_binned_spikes(used_inds(tr_inds),ss);
    Robsxv = all_binned_spikes(used_inds(xv_inds),ss);
            
    gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p);

    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,tr_X);
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
    gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
    
    all_mod_fits(ss) = gqm1;
    
    [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMmodel_eval(all_mod_fits(ss),Robs,tr_X);
    all_mod_LLimp(ss) = (LL-nullLL)/log(2);
    
    if xv_frac > 0
        [xvLL, penLL, pred_rate, G, gint, fgint, xvnullLL] = NMMmodel_eval(all_mod_fits(ss),Robsxv,xv_X);
        all_mod_xvLLimp(ss) = (xvLL-xvnullLL)/log(2);
    end
    
end

%%
for ss = 1:n_probes
    all_filts(ss,:,:) = reshape(all_mod_fits(ss).mods(1).filtK,flen,use_nPix);
end
temp_profiles = squeeze(std(all_filts,[],3));
space_profiles = squeeze(std(all_filts,[],2));
[~,best_lags] = max(temp_profiles,[],2);
[~,best_pos] = max(space_profiles,[],2);
for ss = 1:n_probes
    best_temp_profiles(ss,:) = squeeze(all_filts(ss,:,best_pos(ss)));
    best_space_profiles(ss,:) = squeeze(all_filts(ss,best_lags(ss),:));
end
figure
imagesc(0:dt*1e3:(flen-1)*dt*1e3,1:n_probes,best_temp_profiles);
xlabel('Time (ms)','fontsize',12);
ylabel('Probe','fontsize',12);

figure
imagesc(1:use_nPix,1:n_probes,best_space_profiles);
xlabel('Bar position','fontsize',12);
ylabel('Probe','fontsize',12);