clear all
close all

data_dir = '/media/NTlab_data3/Data/bruce/M320/';
cd(data_dir);
load lemM320Expts.mat
cluster_dir = '~/Analysis/bruce/M320/clustering';
%%
expt_names = cellfun(@(x) x.Header.expname,Expts,'uniformoutput',0);
cur_block_set = find(strcmp(expt_names,'rls.smoothPursuit'));
n_blocks = length(cur_block_set);
params.n_probes = 24;
params.min_trial_dur = 4;
params.dt = 0.01;
params.Fr = 1;
rec_block_range = nan;
rec_type = 'LP';

nf = 400; %fixing number of frames at 400
%%

trial_cnt = 0;

%this is a rediculous way of grabbing all this info...
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_stim_times = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_exvals = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(params.n_probes,1);
all_clust_ids = cell(params.n_probes,1);
all_spk_inds = cell(params.n_probes,1);
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
for ee = 1:n_blocks;
    fprintf('Block %d of %d;  UNMATCHED EXPT TYPE\n',ee,n_blocks);
    cur_block = cur_block_set(ee);
    
    %get spiking data from this block
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:params.n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end
    
    %get trial start and end times
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
        use_trials = find(trial_durs >= params.min_trial_dur);
    end
    
    %add in time offset 
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    %get seed values for each used trial
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
    
    %this stores all trial-varying experimental values
    if isfield(Expts{cur_block}.Trials,'exvals')
        exvals = reshape([Expts{cur_block}.Trials(:).exvals],length(Expts{cur_block}.Trials(1).exvals),[]);
        trial_exvals = exvals(:,id_inds)';
    else
        trial_exvals = nan(length(id_inds),3);
    end
    all_trial_exvals = cat(1,all_trial_exvals,trial_exvals(use_trials,:));
    
    %load in stimulus data
%     fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
%     load(fname);
%     buffer_pix = floor((expt_npix(cur_block) - params.full_nPix)/2);
%     if buffer_pix == -1 %in case there is slightly less pixels than desired by full_nPix, just buffer by a couple zeros
%         for ii = 1:length(left_stim_mats)
%             left_stim_mats{ii} = [zeros(size(left_stim_mats{ii},1),1) left_stim_mats{ii} zeros(size(left_stim_mats{ii},1),1)];
%         end
%         buffer_pix = 0;
%     end
%     cur_use_pix = (1:params.full_nPix) + buffer_pix; %use central pixels
    
    %cycle over trials
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4; %stimulus frame onset times
%         n_frames = size(left_stim_mats{use_trials(tt)},1);
        n_frames = nf;
%         if n_frames > 0
            if length(cur_stim_times) == 1 %if only the first stim time is stored
                cur_stim_times = (cur_stim_times:params.dt*params.Fr:(cur_stim_times + (n_frames-1)*params.dt*params.Fr))'; %create time axis
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
                cur_t_edges = [cur_stim_times; cur_stim_times(end) + params.dt*params.Fr];
            end
%         end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        %         cur_t_edges_up = cur_t_edges(1):up_dt:cur_t_edges(end);
        %         cur_t_axis_up = 0.5*cur_t_edges_up(1:end-1) + 0.5*cur_t_edges_up(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt)); %time since trial onset
        
%         if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > params.min_trial_dur/params.dt %if using this trial
%             use_frames = min(length(cur_stim_times),n_frames); %number of used frames
%             cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix)); %stimulus matrix
            
%             if ~isempty(all_stim_times)
%                 if any(cur_stim_times+cur_toffset < all_stim_times(end))
%                     fprintf('Warn trial %d\n',tt);
%                 end
%             end
            
            %cat values
            all_stim_times = [all_stim_times; cur_stim_times + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
%             all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
%         end
    end
    trial_cnt = trial_cnt + n_trials;
    
    %need to keep track of block time offsets for LP recordings
    if strcmp(rec_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = params.n_probes;
if strcmp(rec_type,'LP')
    clust_params.exclude_adjacent = true;
else
    clust_params.exclude_adjacent = false;
end
[all_binned_mua,all_binned_sua,Clust_data,all_su_spk_times,~,all_mu_spk_times] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params,rec_block_range);
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

%%
sm_win = 5;
all_sm_mua = nan(size(all_binned_mua));
all_sm_sua = nan(size(all_binned_sua));
for ii = 1:params.n_probes
    ii
all_sm_mua(:,ii) = jmm_smooth_1d_cor(all_binned_mua(:,ii),sm_win);
end
for ii = 1:length(SU_probes)
    ii
all_sm_sua(:,ii) = jmm_smooth_1d_cor(all_binned_sua(:,ii),sm_win);
end
%% get trial-based spike responses
unique_trials = unique(all_trialvec);

n_trials = length(unique_trials);

tbt_SUA = nan(nf,n_trials,length(SU_numbers));
tbt_MUA = nan(nf,n_trials,params.n_probes);
for nn = 1:n_trials
   cur_inds = find(all_trialvec == unique_trials(nn));
   if length(cur_inds) ~= nf
       error('parse issue')
   end
   tbt_SUA(:,nn,:) = all_sm_sua(cur_inds,:);
   tbt_MUA(:,nn,:) = all_sm_mua(cur_inds,:);
end

tbt_SUA = bsxfun(@rdivide,tbt_SUA,reshape(nanmean(all_binned_sua),[1 1 length(SU_numbers)]));
tbt_MUA = bsxfun(@rdivide,tbt_MUA,reshape(nanmean(all_binned_mua),[1 1 params.n_probes]));

%%
static_trials = find(all_trial_exvals(:,1) == 0);
drift_trials = find(all_trial_exvals(:,1) > 0);

static_SU_psths = squeeze(nanmean(tbt_SUA(:,static_trials,:),2));
drift_SU_psths = squeeze(nanmean(tbt_SUA(:,drift_trials,:),2));
static_MU_psths = squeeze(nanmean(tbt_MUA(:,static_trials,:),2));
drift_MU_psths = squeeze(nanmean(tbt_MUA(:,drift_trials,:),2));

%%
drift_dirs = unique(all_trial_exvals(:,2));
for ii = 1:length(drift_dirs)
   cur_trials = find(all_trial_exvals(:,2) == drift_dirs(ii) & all_trial_exvals(:,1) > 0);
spec_drift_SU_psths(:,:,ii) = squeeze(nanmean(tbt_SUA(:,cur_trials,:),2));
spec_drift_MU_psths(:,:,ii) = squeeze(nanmean(tbt_MUA(:,cur_trials,:),2));
   cur_trials = find(all_trial_exvals(:,2) == drift_dirs(ii) & all_trial_exvals(:,1) == 0);
spec_static_SU_psths(:,:,ii) = squeeze(nanmean(tbt_SUA(:,cur_trials,:),2));
spec_static_MU_psths(:,:,ii) = squeeze(nanmean(tbt_MUA(:,cur_trials,:),2));
end
%%
use_SUs = setdiff(1:length(SU_probes),7);

load ~/Analysis/bruce/M320/ET_final_imp/PARDRIFT_eyetrack_Rinit_Cprior_ori100.mat dit*
best_mods = dit_mods{end};
%%
close all
tax = (1:nf)*params.dt;

f1 = figure();
f2 = figure();
for ii = 1:length(SU_probes)
    fprintf('Neuron %d\n',ii);
    
    figure(f1);
    clf
    plot(tax,squeeze(spec_drift_SU_psths(:,ii,1)),'b','linewidth',1);
    hold on
    plot(tax,squeeze(spec_drift_SU_psths(:,ii,2)),'m','linewidth',1);
    plot(tax,squeeze(spec_drift_SU_psths(:,ii,3)),'g','linewidth',1);
    plot(tax,squeeze(spec_drift_SU_psths(:,ii,4)),'r','linewidth',1);

        plot(tax,squeeze(spec_static_SU_psths(:,ii,1)),'b--','linewidth',1);
    plot(tax,squeeze(spec_static_SU_psths(:,ii,2)),'m--','linewidth',1);
    plot(tax,squeeze(spec_static_SU_psths(:,ii,3)),'g--','linewidth',1);
    plot(tax,squeeze(spec_static_SU_psths(:,ii,4)),'r--','linewidth',1);

    
    plot(tax,static_SU_psths(:,ii),'k','linewidth',2);

%     if ~isempty(ModData(ii+params.n_probes).bestGQM)
%     plot_NMM_filters_1d(ModData(ii+params.n_probes).bestGQM,[],[],[],f2);
    if ~isempty(best_mods(ii+params.n_probes))
    plot_NMM_filters_1d(best_mods(ii+params.n_probes),[],[],[],f2);
    else
        figure(f2)
        clf
    end
    pause
end

%%
close all
tax = (1:nf)*params.dt;

f1 = figure();
for ii = 1:params.n_probes
    fprintf('Neuron %d\n',ii);
    
    figure(f1);
    clf
    plot(tax,squeeze(spec_drift_MU_psths(:,ii,1)),'b','linewidth',1);
    hold on
    plot(tax,squeeze(spec_drift_MU_psths(:,ii,2)),'m','linewidth',1);
    plot(tax,squeeze(spec_drift_MU_psths(:,ii,3)),'g','linewidth',1);
    plot(tax,squeeze(spec_drift_MU_psths(:,ii,4)),'r','linewidth',1);
    
    plot(tax,static_MU_psths(:,ii),'k','linewidth',2);

%     if ~isempty(ModData(ii+params.n_probes).bestGQM)
%     plot_NMM_filters_1d(ModData(ii+params.n_probes).bestGQM,[],[],[],f2);
%     else
%         figure(f2)
%         clf
%     end
    pause
end


%%
or_tuning = find(strcmp(expt_names,'rls.orRC'));
cur_block_set = or_tuning;
n_blocks = length(cur_block_set);

nf = 200;
params.min_trial_dur = 2;

%this is a rediculous way of grabbing all this info...
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_stim_times = [];
all_stim_or = [];
all_blockvec = [];
all_trialvec = [];
all_trial_or = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(params.n_probes,1);
all_clust_ids = cell(params.n_probes,1);
all_spk_inds = cell(params.n_probes,1);
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
for ee = 1:n_blocks;
    fprintf('Block %d of %d;  UNMATCHED EXPT TYPE\n',ee,n_blocks);
    cur_block = cur_block_set(ee);
    
    %get spiking data from this block
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:params.n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end
    
    %get trial start and end times
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
        use_trials = find(trial_durs >= params.min_trial_dur);
    end
    
    %add in time offset 
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    %get seed values for each used trial
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
    
    %load in stimulus data
%     fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
%     load(fname);
%     buffer_pix = floor((expt_npix(cur_block) - params.full_nPix)/2);
%     if buffer_pix == -1 %in case there is slightly less pixels than desired by full_nPix, just buffer by a couple zeros
%         for ii = 1:length(left_stim_mats)
%             left_stim_mats{ii} = [zeros(size(left_stim_mats{ii},1),1) left_stim_mats{ii} zeros(size(left_stim_mats{ii},1),1)];
%         end
%         buffer_pix = 0;
%     end
%     cur_use_pix = (1:params.full_nPix) + buffer_pix; %use central pixels
    
    %cycle over trials
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4; %stimulus frame onset times
        cur_stim_or = Expts{cur_block}.Trials(use_trials(tt)).or;
%         n_frames = size(left_stim_mats{use_trials(tt)},1);
        n_frames = nf;
%         if n_frames > 0
%             if length(cur_stim_times) == 1 %if only the first stim time is stored
                cur_stim_times = (cur_stim_times:params.dt*params.Fr:(cur_stim_times + (n_frames-1)*params.dt*params.Fr))'; %create time axis
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
                cur_t_edges = [cur_stim_times; cur_stim_times(end) + params.dt*params.Fr];
%             end
%         end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        %         cur_t_edges_up = cur_t_edges(1):up_dt:cur_t_edges(end);
        %         cur_t_axis_up = 0.5*cur_t_edges_up(1:end-1) + 0.5*cur_t_edges_up(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt)); %time since trial onset
        
%         if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > params.min_trial_dur/params.dt %if using this trial
%             use_frames = min(length(cur_stim_times),n_frames); %number of used frames
%             cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix)); %stimulus matrix
            
%             if ~isempty(all_stim_times)
%                 if any(cur_stim_times+cur_toffset < all_stim_times(end))
%                     fprintf('Warn trial %d\n',tt);
%                 end
%             end
            
            %cat values
            all_stim_times = [all_stim_times; cur_stim_times + cur_toffset];
            all_stim_or = [all_stim_or; cur_stim_or(:)];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
%             all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
%         end
    end
    trial_cnt = trial_cnt + n_trials;
    
    %need to keep track of block time offsets for LP recordings
    if strcmp(rec_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end


%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = params.n_probes;
if strcmp(rec_type,'LP')
    clust_params.exclude_adjacent = true;
else
    clust_params.exclude_adjacent = false;
end
[all_binned_mua,all_binned_sua,Clust_data,all_su_spk_times,~,all_mu_spk_times] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params,rec_block_range);
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

all_binned_mua = bsxfun(@rdivide,all_binned_mua,nanmean(all_binned_mua));
all_binned_sua = bsxfun(@rdivide,all_binned_sua,nanmean(all_binned_sua));


%%

un_oris = unique(all_stim_or);
for_win = 10;
back_win = 10;

or_ta_SU = nan(for_win+back_win + 1,length(un_oris),length(SU_probes));
or_ta_MU = nan(for_win+back_win + 1,length(un_oris),params.n_probes);
for ii = 1:length(un_oris)
   cur_set = find(all_stim_or == un_oris(ii));
   or_ta_SU(:,ii,:) = get_event_trig_avg_v3(all_binned_sua,cur_set,back_win,for_win);
   or_ta_MU(:,ii,:) = get_event_trig_avg_v3(all_binned_mua,cur_set,back_win,for_win);
end










