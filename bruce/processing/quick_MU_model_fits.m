clear all
close all

monName = 'jbe';
Expt_name = 'M010';
data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
% data_dir = ['~/Data/bruce/' Expt_name];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
rec_type = 'LP';

use_SUs = false;

rec_number = 1;
% use_block_range =1:27;

if rec_number > 1
    cluster_dir = [cluster_dir sprintf('/rec%d',rec_number)];
end

%%

min_trial_dur = 2;
trial_dur = 4;
use_nPix = 20;
dt = 0.01;
flen = 12;
Fr = 1;
n_probes = 24;
spatial_usfac = 1;

%%
cd(data_dir);
ExptFileName = strcat(monName,Expt_name,'Expts.mat');
fprintf('Loading existing file %s\n',ExptFileName);
load(ExptFileName);
cd stims/
load('./stim_data.mat');
load('./expt_data.mat');

%%
% include_expts = {'rls.froNoise'};
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa','rls.FaXwi','rls.FaXwiXimi','rls.AllSacB','rls.froNoise'};
usable_expts = find(cellfun(@(x) ~isempty(x),Expts));
expt_names = cellfun(@(x) x.Header.expname,Expts(usable_expts),'uniformoutput',0);

cur_block_set = usable_expts(cellfun(@(x) any(strcmp(x,include_expts)),expt_names));

full_nPix = unique(expt_npix(cur_block_set));
if length(full_nPix) > 1
    warning('multiple npix detected');
    full_nPix = min(full_nPix);
end

eds = cellfun(@(x) x.Stimvals.ed,Expts);
if exist('use_block_range','var')
    cur_block_set(~ismember(cur_block_set,use_block_range)) = [];
end

n_blocks = length(cur_block_set);

%% Spatial resolution
all_dws = cellfun(@(x) x.Stimvals.dw,Expts(cur_block_set));
base_sp_dx = mode(all_dws);
if length(unique(all_dws)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac; %model dx in deg
ov_RF_pos = Expts{cur_block_set(1)}.Stimvals.rf(1:2);

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_wi = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_trial_rptframes = [];

trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
for ee = 1:n_blocks;
    fprintf('Block %d of %d;  UNMATCHED EXPT TYPE\n',ee,n_blocks);
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
%     trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
    if buffer_pix == -1
        for ii = 1:length(left_stim_mats)
            left_stim_mats{ii} = [zeros(size(left_stim_mats{ii},1),1) left_stim_mats{ii} zeros(size(left_stim_mats{ii},1),1)];
        end
        buffer_pix = 0;
    end
    cur_use_pix = (1:full_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    cur_nrpt_frames = zeros(n_trials,1);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if isfield(Expts{cur_block}.Trials(use_trials(tt)),'rptframes')
            cur_nrpt_frames(tt) = length(Expts{cur_block}.Trials(use_trials(tt)).rptframes);
        end
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
                cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
            end
        end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times+cur_toffset < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
    all_trial_rptframes = [all_trial_rptframes; cur_nrpt_frames];
    
    %need to keep track of block time offsets for LP recordings
    if strcmp(rec_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

%% BIN SPIKES FOR MU AND SU
if ~use_SUs
    all_binned_mua = get_quickbinned_mu(all_spk_times,all_clust_ids,all_t_axis,all_t_bin_edges,all_bin_edge_pts,clust_params);
else
clust_params.n_probes = n_probes;
%     if strcmp(rec_type,'LP')
%         clust_params.exclude_adjacent = true;
%     else
        clust_params.exclude_adjacent = false;
%     end
    [all_binned_mua,all_binned_sua,Clust_data,all_su_spk_times,~,all_mu_spk_times] = ...
        get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
        all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
    SU_probes = Clust_data.SU_probes;
    SU_numbers = Clust_data.SU_numbers;
end
%%
full_nPix_us = spatial_usfac*full_nPix;
if spatial_usfac > 1
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:spatial_usfac
            all_stimmat_up(:,spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
        end
    end
elseif spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
end

%select subset of pixels used for model fitting
buffer_pix = floor((full_nPix - use_nPix)/2);
[Xinds_up,~] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
cur_use_pix = (1/spatial_usfac:1/spatial_usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);

%% CREATE STIMULUS MATRIX
stim_params = NIMcreate_stim_params([flen full_nPix_us],dt);
all_Xmat = create_time_embedding(all_stimmat_up,stim_params);
all_Xmat = all_Xmat(:,use_kInds_up);

%%
% DEFINE DATA USED FOR ANALYSIS
beg_buffer = 0.2;
end_buffer = 0.05;
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);

% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

%%
tr_X{2} = abs(all_Xmat(used_inds,:));
tr_X{1} = all_Xmat(used_inds,:);

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
init_Xtargs = [2; 1];

for ss = 1:n_probes;
    fprintf('Fitting model for MU %d of %d\n',ss,n_probes);
    Robs = all_binned_mua(used_inds,ss);
    
    gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p);
    
    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,tr_X);
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
    gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
    
    all_mod_fits(ss) = gqm1;
    
    [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMmodel_eval(all_mod_fits(ss),Robs,tr_X);
    all_mod_LLimp(ss) = (LL-nullLL)/log(2);
    
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

%%
if use_SUs
    
    n_SUs = size(all_binned_sua,2);
    
    silent = 1;
    base_lambda_d2XT = 5;
    base_lambda_L1 = 0;
    init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-9;
    
    n_stim_filts = 3;
    mod_signs = ones(1,n_stim_filts);
    NL_types = [{'lin' 'quad' 'quad'}];
    init_d2XT = [ones(n_stim_filts,1)];
    init_L2 = [zeros(n_stim_filts,1);];
    init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT);
    
    for ss = 1:n_SUs;
        fprintf('Fitting model for SU %d of %d\n',ss,n_SUs);
        Robs = all_binned_sua(used_inds,ss);
        uinds = find(~isnan(Robs));
        if ~isempty(uinds)
            
            gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params);
            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],uinds,silent,init_optim_p);
            
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,tr_X);
            gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
            gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],uinds,silent);
            
            all_SU_fits(ss) = gqm1;
            
            [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMmodel_eval(all_SU_fits(ss),Robs,tr_X);
            all_SU_LLimp(ss) = (LL-nullLL)/log(2);
            
        end
    end
end