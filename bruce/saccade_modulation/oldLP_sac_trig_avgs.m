clear all
close all

Expt_name = 'M239';
Expt_num = str2num(Expt_name(2:end));
switch Expt_num
    case 232
        bar_ori = 50;
    case 235 
        bar_ori = 30;
    case 239
        bar_ori = 130;
end

Expt_name = sprintf('M%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load ./random_bar_eyedata_ftime.mat bar_expts

load(['lemM' num2str(Expt_num) 'Expts.mat']);
load ./bar_params.mat

save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end
n_probes = 24;
expt_type = 'LP';

ignore_blocks = [];
%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.005;

min_trial_dur = 1;
beg_buffer = 0.2;
end_buffer = 0.05;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.005/dt);

if strcmp(Expt_name,'G081')
    trial_dur = 2;
else
    trial_dur = 4;
end

%% LOAD EXPTS STRUCT
cd(data_dir)
if Expt_name(1) == 'G'
load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
else
load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
end


%% LOAD OVERALL SU DATA
% LOAD REFCLUSTERS
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
fname = [cluster_dir '/final_cluster.mat'];
if exist(fname,'file')
    load(fname);
    SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
    for ii = 1:length(SU_numbers)
        SU_tot_nblocks(ii) = sum(SU_ID_mat(:) == SU_numbers(ii));
    end
    fprintf('%d SUs Clustered\n',length(SU_numbers));
    
else
    disp('No Cluster data found.');
end

%% PARSE EXPTS STRUCT
cur_block_set = bar_expts;
cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

n_blocks = length(cur_block_set);

gsac_amp = Expts{cur_block_set(1)}.Stimvals.Fs;
gsac_thresh = gsac_amp/3;

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_ttill_end = [];
all_blockvec = [];
all_trialvec = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_trial_blocknums = [];
all_bin_edge_pts = [];

all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:length(cur_block_set);
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
    if length(un_ids) < length(trial_ids)
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);

    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        cur_ttill_end = trial_end_times(use_trials(tt)) - cur_t_axis;
        
        all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
        all_tsince_start = [all_tsince_start; cur_tsince_start];
        all_ttill_end = [all_ttill_end; cur_ttill_end];
        all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    end
    trial_cnt = trial_cnt + n_trials;
    if strcmp(expt_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

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

% double_spike_buffer = 3; %number of samples (in either direction) to exclude double spikes from adjacent-probe SUs
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
%     cur_mua_inds = find(all_clust_ids{cc} >= 1);
    cur_mua_inds = find(all_clust_ids{cc} >= 0);
    
%     %remove spikes from isolated SUs on the same probe from the MU
%     for ss = 1:length(unique_su_nums)
%         cur_mua_inds(ismember(all_spk_inds{cc}(cur_mua_inds),all_su_spk_inds{unique_su_nums(ss)})) = [];
%     end
% 
%     nearby_probes = [cc-1 cc+1]; nearby_probes(nearby_probes < 1 | nearby_probes > n_probes) = [];
%     cur_set = find(ismember(SU_block_probes,nearby_probes));
%     if ~isempty(cur_set)
%         [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
%     else
%         cur_SS = [];
%     end
%     unique_su_nums = unique(cur_SS); %set of SU numbers picked up on adjacent probes
%     if ~isempty(unique_su_nums)
%         double_spikes = [];
%         for ss = 1:length(unique_su_nums)
%             cur_blocked_inds = bsxfun(@plus,all_su_spk_inds{unique_su_nums(ss)},-double_spike_buffer:double_spike_buffer);
%             double_spikes = [double_spikes; find(ismember(all_spk_inds{cc}(cur_mua_inds),cur_blocked_inds))];
%         end
%         fprintf('Eliminating %d of %d double spikes in MUA\n',length(double_spikes),length(cur_mua_inds));
%         cur_mua_inds(double_spikes) = [];
%     end
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end

%% SMOOTH AND NORMALIZE SPIKING DATA
all_mua_rate = nan(size(all_binned_mua));
mua_block_mean_rates = nan(length(cur_block_set),n_probes);
mua_block_n_spikes = nan(length(cur_block_set),n_probes);
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
    for cc = 1:n_probes
        all_mua_rate(cur_block_inds,cc) = jmm_smooth_1d_cor(all_binned_mua(cur_block_inds,cc),mua_sm_sig);
    end
    mua_block_mean_rates(ee,:) = mean(all_mua_rate(cur_block_inds,:));
    mua_block_n_spikes(ee,:) = sum(all_binned_mua(cur_block_inds,:));
    end
end

all_sua_rate = nan(size(all_binned_sua));
sua_block_mean_rates = nan(length(cur_block_set),length(SU_numbers));
sua_block_n_spikes = nan(length(cur_block_set),length(SU_numbers));
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
    for ss = 1:length(SU_numbers)
        all_sua_rate(cur_block_inds,ss) = jmm_smooth_1d_cor(all_binned_sua(cur_block_inds,ss),sua_sm_sig);
    end
    sua_block_mean_rates(ee,:) = mean(all_sua_rate(cur_block_inds,:));
    sua_block_n_spikes(ee,:) = sum(all_binned_sua(cur_block_inds,:));
    end
end

%normalized firing rates (smoothed)
all_sua_rate_norm = nan(size(all_sua_rate));
all_mua_rate_norm = nan(size(all_mua_rate));
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
    all_sua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_sua_rate(cur_block_inds,:),...
        sua_block_mean_rates(ee,:));
    all_mua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_mua_rate(cur_block_inds,:),...
        mua_block_mean_rates(ee,:));
    end
end

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & all_ttill_end >= end_buffer);

%% PROCESS EYE TRACKING DATA
if strcmp(expt_type,'LP')
    expt_bar_ori = ones(size(cur_block_set))*bar_ori;
end
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori,used_inds);

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);

par_thresh = 4;
orth_thresh = 1.2;
[out_of_range] = detect_bad_fixation(corrected_eye_vals_interp,all_trialvec,used_inds,par_thresh,orth_thresh);
fract_out = length(out_of_range)/length(used_inds);
fprintf('Eliminating %.4f of data out of window\n',fract_out);
used_inds(ismember(used_inds,out_of_range)) = [];
NT = length(used_inds);

%% PROCESS SACCADE STATS

%interpolate saccade start times and get rid of saccades that aren't within
%the t-axis binning
sac_start_times = [saccades(:).start_time];
sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
sac_start_inds(isnan(sac_start_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(sac_start_inds)');
bad_sacs = find(isnan(sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
saccades(bad_sacs) = [];
sac_start_inds(bad_sacs) = [];

sac_stop_times = [saccades(:).stop_time];
sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
sac_stop_inds(isnan(sac_stop_inds)) = 1;

sac_peak_times = [saccades(:).peak_time];
sac_peak_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_peak_times));
sac_peak_inds(isnan(sac_peak_inds)) = 1;

max_sac_dur = 0.1;
sac_durs = [saccades(:).duration];
is_blink = sac_durs > max_sac_dur;

isis = [saccades(:).isi];
is_sacburst = false(length(saccades),1);
is_sacburst(isis(1:end-1) < 0.15 | isis(2:end) < 0.15) = true;

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps < 1;

sac_deltaX = [saccades(:).post_Lx] - [saccades(:).pre_Lx];
sac_deltaY = [saccades(:).post_Ly] - [saccades(:).pre_Ly];
sac_postX = [saccades(:).post_Lx];
sac_postY = [saccades(:).post_Ly];

%calculate amplitude of saccades along the guided saccade path
delta_sacpar = abs(sac_deltaX);
is_gsac = delta_sacpar' >= gsac_thresh;


%% PICK OUT SACCADES FOR ANALYSIS

%define which saccades to use
used_msacs = find(is_micro & ~is_blink & ismember(sac_start_inds,used_inds));
used_gsacs = find(is_gsac' & ~is_blink & ismember(sac_start_inds,used_inds));

%microsacs excluding bursts
non_burst_msacs = used_msacs(~is_sacburst(used_msacs));
burst_msacs = used_msacs(is_sacburst(used_msacs));


%% COMPUTE TRIG AVGS FOR MUA
nboot = [];
%set trial numbers to Inf so they don't get included in trig averaging
used_trialvec = ones(size(all_trialvec))*Inf;
used_trialvec(used_inds) = all_trialvec(used_inds);
clear mua_data
for cc = 1:n_probes
    
    fprintf('Computing trig avgs for MUA %d of %d\n',cc,n_probes);
    
    %general averages
    [mua_data(cc).msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).gsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0);
    
    [mua_data(cc).msac_end_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_stop_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).gsac_end_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_stop_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0);
        
    mua_data(cc).avg_rate = mean(all_binned_mua(used_inds,cc));
    mua_data(cc).tot_nspikes = sum(all_binned_mua(used_inds,cc));
    
end

%% COMPUTE TRIG AVGS FOR SUA
nboot = 100;
clear sua_data
for ss = 1:length(SU_numbers)
    
    fprintf('Computing trig avgs for SU %d of %d\n',ss,length(SU_numbers));
    
    cur_use_inds = used_inds(~isnan(all_binned_sua(used_inds,ss)));
    used_trialvec = ones(size(all_trialvec))*Inf;
    used_trialvec(cur_use_inds) = all_trialvec(cur_use_inds);
    
    cur_used_blocks = find(~isnan(SU_block_probes(ss,:)));
    cur_inds = find(ismember(all_blockvec,cur_used_blocks));
    if length(cur_used_blocks) >= 3
        sua_data(ss).used = true;
                
        %general averages
        [sua_data(ss).msac_avg,lags,sua_data(ss).msac_sem,sua_data(ss).nused_msac] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).gsac_avg,lags,sua_data(ss).gsac_sem,sua_data(ss).nused_gsac] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        
        [sua_data(ss).msac_end_avg,lags,sua_data(ss).msac_end_sem,sua_data(ss).nused_msac] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_stop_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).gsac_end_avg,lags,sua_data(ss).gsac_end_sem,sua_data(ss).nused_gsac] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_stop_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
    
        sua_data(ss).avg_rate = mean(all_binned_sua(cur_use_inds,ss));
        sua_data(ss).tot_nspikes = sum(all_binned_sua(cur_use_inds,ss));
        %start and stop trig avgs
        %micro sac burst and non-burst
    else
        sua_data(ss).used = false;
    end
end

%%
cd(save_dir)
sname = 'sac_trig_avg_data';
save(sname,'sua_data','mua_data','lags','dt');
%%
