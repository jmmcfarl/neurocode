%
clear all

global Expt_name monk_name bar_ori rec_type

Expt_name = 'M014';
monk_name = 'jbe';
bar_ori = 40; %bar orientation to use (only for UA or single-ori-LP recs)
rec_number = 1;

Expt_num = str2num(Expt_name(2:end));

data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
load(Edata_file);

%is this a laminar probe or utah array rec?
firste = find(cellfun(@(x) ~isempty(x),Expts),1);
if strcmp(Expts{firste}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{firste}.Header.DataType,'Spike2')
    rec_type = 'LP';
end


if strcmp(Expt_name,'M011')
    rec_block_range =1:21; %M011
elseif strcmp(Expt_name,'M012')
    rec_block_range = 1:27;
else
    rec_block_range = nan;
end

%%
if strcmp(rec_type,'LP')
    params.n_probes = 24;
elseif strcmp(rec_type,'UA')
    params.n_probes = 96;
end
if strcmp(monk_name,'jbe')
    params.good_coils = [1 0];
        params.use_coils = [0 0];
elseif strcmp(monk_name,'lem');
    params.good_coils = [1 1];
    params.use_coils = [1 1];
end

load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
if rec_number > 1
    rec_block_range = setdiff(1:length(Expts),rec_block_range);
    cluster_dir = [cluster_dir sprintf('/rec%d',rec_number)];
end

%%
min_trial_dur = 0.75;
stim_fs = 100; %in Hz
dt = 0.005;
Fr = 1;

%exclude data at beginning and end of each trial
beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

%% SELECT BLOCKS FOR ANALYSIS
expt_names = cellfun(@(x) x.Header.expname,Expts,'uniformoutput',0);
cur_block_set = find(strcmp(expt_names,'grating.tf'));
n_blocks = length(cur_block_set);

grating_sf = Expts{cur_block_set(1)}.Stimvals.sf;
grating_or = Expts{cur_block_set(1)}.Stimvals.or;
nf = Expts{cur_block_set(1)}.Stimvals.nf;
%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_tf = [];
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
    fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,n_blocks);
    cur_block = cur_block_set(ee);
    
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:params.n_probes
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
    
    trial_tf = [Expts{cur_block}.Trials(:).tf];
    trial_tf = trial_tf(id_inds);
    all_trial_tf = cat(1,all_trial_tf,trial_tf(:));
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        if length(cur_t_axis) ~= nf/stim_fs/dt
            error('Incomplete trial, havent accounted for this yet');
        end
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
        all_tsince_start = [all_tsince_start; cur_tsince_start];
        all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
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
trial_mua = reshape(all_binned_mua,nf/stim_fs/dt,[],params.n_probes);
trial_sua = reshape(all_binned_sua,nf/stim_fs/dt,[],length(SU_numbers));

%% LOAD LFP DATA
cd(data_dir)
if strcmp(rec_type,'LP')
    Fs = 1000;
    dsf = 5;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    new_niqf = Fsd/2;
    [bb,aa] = butter(2,0.8*new_niqf/niqf);
    use_lfps = 1:1:params.n_probes;
elseif strcmp(rec_type,'UA')
    Fs = 400.0032;
    dsf = 2;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    new_niqf = Fsd/2;
    [bb,aa] = butter(4,0.8*new_niqf/niqf);
    use_lfps = SU_probes;
    %     use_lfps = [28 SU_probes];
end

full_lfps = [];
full_lfp_taxis = [];
cur_toffset = 0;
ublock_set = 1:length(cur_block_set);
for ee = ublock_set
    
    fprintf('Loading LFPs, Expt %d of %d\n',ee,length(cur_block_set));
    
    if strcmp(rec_type,'LP')
        fname = sprintf('%sM%.3dA.%d.lfp.mat',monk_name,Expt_num,cur_block_set(ee));
        load(fname);
        
        tlens = arrayfun(@(X) length(X.ftime),LFP.Trials);
        bad_trials = find(tlens == 0);
        LFP.Trials(bad_trials) = [];
        n_trials(ee) = length(LFP.Trials);
        lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
        lfp_trial_ends = [LFP.Trials(:).End]/1e4;
        expt_lfp_t_axis = [];
        expt_lfps = [];
        for tt = 1:n_trials(ee)
            %         tt
            cur_npts = size(LFP.Trials(tt).LFP,1);
            cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
            cur_t_axis = (lfp_trial_starts(tt):1/Fs:cur_t_end(tt)) + cur_toffset;
            
            if ~isempty(expt_lfp_t_axis)
                cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
            else
                cur_sp = 1;
            end
            cur_t_axis = cur_t_axis(cur_sp:end);
            if length(cur_t_axis) > 50 & size(LFP.Trials(tt).LFP,2) == params.n_probes
                cur_LFP = double([LFP.Trials(tt).LFP]);
                cur_LFP = cur_LFP(cur_sp:end,use_lfps);
                if dsf > 1
                    cur_LFP = filtfilt(bb,aa,cur_LFP);
                    cur_LFP = downsample(cur_LFP,dsf);
                    cur_t_axis = downsample(cur_t_axis,dsf);
                end
                if size(cur_LFP,2) == length(use_lfps)
                    expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
                    expt_lfps = [expt_lfps; cur_LFP];
                end
            end
        end
    else
        lfp_fname = sprintf('Expt%d_LFP.mat',cur_block_set(ee));
        load(lfp_fname);
        
        cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
        if dsf > 1
            cur_lfps = filtfilt(bb,aa,cur_lfps);
            expt_lfps = downsample(cur_lfps,dsf);
            expt_lfp_t_axis = downsample(lfp_t_ax',dsf);
        end
    end
    
    cur_uset = find(all_blockvec == ee);
    if ~isempty(cur_uset)
        uinds = find(expt_lfp_t_axis >= all_t_axis(cur_uset(1)) & expt_lfp_t_axis <= all_t_axis(cur_uset(end)));
        full_lfps = cat(1,full_lfps,expt_lfps(uinds,:));
        full_lfp_taxis = cat(1,full_lfp_taxis,expt_lfp_t_axis(uinds));
    end
    cur_toffset = trial_toffset(ee);
end

%%
interp_lfps = interp1(full_lfp_taxis,full_lfps,all_t_axis);
trial_lfps = reshape(interp_lfps,nf/stim_fs/dt,[],params.n_probes);
trial_lfps(isnan(trial_lfps)) = 0;
%%
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))
tf = 4;
trial_set = find(all_trial_tf == tf);

params.Fs = 1/dt;
params.tapers = [2 3];
params.trialave = 1;

clear S
for ii = 1:params.n_probes
   [S(ii,:),f] = mtspectrumc(squeeze(trial_lfps(:,trial_set,ii)),params); 
end

%%
tf = 8;
trial_set = find(all_trial_tf == tf);

[bb,aa] = butter(2,[tf/1.5 tf*1.5]/(new_niqf));

filtered_lfps = nan(nf/stim_fs/dt,length(trial_set),params.n_probes);
lfp_phase = nan(nf/stim_fs/dt,length(trial_set),params.n_probes);
for ii = 1:length(trial_set)
   temp = filtfilt(bb,aa,trial_lfps(:,trial_set(ii),:)); 
   filtered_lfps(:,ii,:) = temp;
   lfp_phase(:,ii,:) = angle(hilbert(squeeze(temp)));
end

%%
% use_lfp_ch = 12;
% cur_trial_sua = trial_sua(:,trial_set,:);
% sorted_sua = cur_trial_sua;
% cur_trial_mua = trial_mua(:,trial_set,:);
% sorted_mua = cur_trial_mua;
% for ii = 1:nf
%     [a,b] = sort(lfp_phase(ii,:,use_lfp_ch));
%     sorted_sua(ii,:,:) = cur_trial_sua(ii,b,:);
%     sorted_mua(ii,:,:) = cur_trial_mua(ii,b,:);
% end


%%
avg_mua = squeeze(mean(trial_mua,3));
sm_avg_mua = avg_mua;
for ii = 1:length(all_trial_tf)
    sm_avg_mua(:,ii) = jmm_smooth_1d_cor(avg_mua(:,ii),2);
end

cur_avg_mua = sm_avg_mua(:,trial_set);

frame_period = round(1/tf/dt);
chunks_per_trial = floor(nf/stim_fs/dt/frame_period)-1;
pts_per_trial = chunks_per_trial*frame_period;
avg_mua_fp = reshape(cur_avg_mua(frame_period + (1:pts_per_trial),:),frame_period,[]);

use_lfp_ch = 15;
avg_phase_signal = reshape(lfp_phase(frame_period + (1:pts_per_trial),:,use_lfp_ch),frame_period,[]);
cent_pt = round(frame_period/2);
central_phase = avg_phase_signal(cent_pt,:);
[a,phase_sort] = sort(central_phase);


n_phase_bins = 20;
bin_edges = round(linspace(1,size(avg_mua_fp,2),n_phase_bins+1));
binned_avg_phase = nan(n_phase_bins,frame_period);
binned_avg_mua = nan(n_phase_bins,frame_period);
for ii = 1:n_phase_bins
    cur_range = bin_edges(ii):bin_edges(ii+1);
    binned_avg_phase(ii,:) = mean(avg_phase_signal(:,phase_sort(cur_range)),2);
    binned_avg_mua(ii,:) = mean(avg_mua_fp(:,phase_sort(cur_range)),2);
end

hangle = angle(hilbert(binned_avg_phase'))';
[~,binned_central_phase_frames] = max(hangle,[],2);
% binned_central_phase_frames = round(binned_central_phase/(2*pi)*frame_period);
binned_avg_phase_aligned = binned_avg_phase;
binned_avg_mua_aligned = binned_avg_mua;
for ii = 1:n_phase_bins
    binned_avg_phase_aligned(ii,:) = circshift(binned_avg_phase(ii,:)',-binned_central_phase_frames(ii))';
    binned_avg_mua_aligned(ii,:) = circshift(binned_avg_mua(ii,:)',-binned_central_phase_frames(ii))';
end

%%
avg_sua = trial_sua(:,:,1);
sm_avg_sua = avg_sua;
for ii = 1:length(all_trial_tf)
    sm_avg_sua(:,ii) = jmm_smooth_1d_cor(avg_sua(:,ii),2);
end

% cur_avg_sua = sm_avg_sua(:,trial_set);
cur_avg_sua = avg_sua(:,trial_set);
frame_period = round(1/tf/dt);
chunks_per_trial = floor(nf/stim_fs/dt/frame_period)-1;
pts_per_trial = chunks_per_trial*frame_period;
avg_sua_fp = reshape(cur_avg_sua(frame_period + (1:pts_per_trial),:),frame_period,[]);

use_lfp_ch = 15;
avg_phase_signal = reshape(lfp_phase(frame_period + (1:pts_per_trial),:,use_lfp_ch),frame_period,[]);
cent_pt = round(frame_period/2);
central_phase = avg_phase_signal(cent_pt,:);
[a,phase_sort] = sort(central_phase);

hangle = angle(hilbert(avg_phase_signal));
[~,central_phase_frames] = max(hangle);
% [~,central_phase_frames] = max(avg_phase_signal);
% binned_central_phase_frames = round(binned_central_phase/(2*pi)*frame_period);
avg_phase_aligned = avg_phase_signal;
avg_sua_aligned = avg_sua_fp;
for ii = 1:size(hangle,2)
    avg_phase_aligned(:,ii) = circshift(avg_phase_signal(:,ii),-central_phase_frames(ii));
    avg_sua_aligned(:,ii) = circshift(avg_sua_fp(:,ii),-central_phase_frames(ii));
end



% n_phase_bins = 20;
% bin_edges = round(linspace(1,size(avg_mua_fp,2),n_phase_bins+1));
% binned_avg_phase = nan(n_phase_bins,frame_period);
% binned_avg_sua = nan(n_phase_bins,frame_period);
% for ii = 1:n_phase_bins
%     cur_range = bin_edges(ii):bin_edges(ii+1);
%     binned_avg_phase(ii,:) = mean(avg_phase_signal(:,phase_sort(cur_range)),2);
%     binned_avg_sua(ii,:) = mean(avg_sua_fp(:,phase_sort(cur_range)),2);
% end
% 
% hangle = angle(hilbert(binned_avg_phase'))';
% [~,binned_central_phase_frames] = max(hangle,[],2);
% % binned_central_phase_frames = round(binned_central_phase/(2*pi)*frame_period);
% binned_avg_phase_aligned = binned_avg_phase;
% binned_avg_sua_aligned = binned_avg_sua;
% for ii = 1:n_phase_bins
%     binned_avg_phase_aligned(ii,:) = circshift(binned_avg_phase(ii,:)',-binned_central_phase_frames(ii))';
%     binned_avg_sua_aligned(ii,:) = circshift(binned_avg_sua(ii,:)',-binned_central_phase_frames(ii))';
% end
% 
