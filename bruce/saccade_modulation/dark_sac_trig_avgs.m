clear all
close all

Expt_name = 'G099';

Expt_num = str2num(Expt_name(2:end));
data_dir = ['/media/NTlab_data2/Data/bruce/' Expt_name];
cluster_dir = data_dir;

save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end

if Expt_name(1) == 'G'
    n_probes = 96;
    expt_type = 'UA';
elseif Expt_name(1) == 'M'
    n_probes = 24;
    expt_type = 'LP';
    switch Expt_num
        case 266
            bar_ori = 80;
        case 270
            bar_ori = 60;
        case 275
            bar_ori = 135;
        case 277
            bar_ori = 70;
    end
    
else
    error('Unrecognized experiment name');
end

if Expt_num == 275 || Expt_num == 277
    rpt_seed = 1001; %M275 M277
else
    rpt_seed = 1e4; %m270 and 266
end

if strcmp(Expt_name,'M270')
    ignore_blocks = [5];
elseif strcmp(Expt_name,'G087')
    ignore_blocks = [15]; %only 6 trials and causes problems
elseif strcmp(Expt_name,'G093')
    ignore_blocks = [28]; %only 6 trials and causes problems
else
    ignore_blocks = [];
end

use_coils = [1 0];

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.005;

min_trial_dur = 1;
beg_buffer = 0.25;
end_buffer = 0.05;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.01/dt);

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

%% PARSE EXPTS STRUCT
is_sim_msac_expt = false;
include_expts = 'none.fc';

has_data = cellfun(@(x) length(x),Expts) > 0;
expt_names(has_data) = cellfun(@(x) x.Header.expname,Expts(has_data),'UniformOutput',false);
expt_backstim(has_data) = cellfun(@(x) x.Stimvals.Bs,Expts(has_data),'UniformOutput',false);
expt_imback = strcmp(expt_backstim,'image');
expt_dd(has_data) = cellfun(@(x) x.Stimvals.dd,Expts(has_data),'UniformOutput',true);
expt_bar_ori(has_data) = cellfun(@(x) x.Stimvals.or,Expts(has_data),'UniformOutput',true);
expt_Fr(has_data) = cellfun(@(x) x.Stimvals.Fr,Expts(has_data),'UniformOutput',true);
expt_sac_dir(has_data) = mod(cellfun(@(x) x.Stimvals.Fa,Expts(has_data),'UniformOutput',true),180);
expt_sac_amp(has_data) = cellfun(@(x) x.Stimvals.Fs,Expts(has_data),'UniformOutput',true);
expt_ce(has_data) = cellfun(@(x) x.Stimvals.ce,Expts(has_data),'UniformOutput',true);
expt_ijump(has_data) = cellfun(@(x) x.Stimvals.ijump,Expts(has_data),'UniformOutput',true);

included_type(has_data) = ismember(expt_names(has_data),include_expts);
cur_block_set = find(included_type);

cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

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
all_trial_back = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_trial_blocknums = [];
all_bin_edge_pts = [];

all_spk_times = cell(n_probes,1);
% all_spk_inds = cell(n_probes,1);
% all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:length(cur_block_set);
    cur_block = cur_block_set(ee);
    fprintf('Expt %s Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_name,ee,length(cur_block_set));
    
    %     fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
%     load(fname,'Clusters');
    fname = [cluster_dir sprintf('/Expt%dClusterTimes.mat',cur_block)];
    load(fname);
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times' + cur_toffset);
%         all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
%         all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
%         all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
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
    
    cur_t_edges = (trial_start_times(1):dt:trial_end_times(end))';
    cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
    all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
    all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
    all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
    all_trialvec = [all_trialvec; ones(size(cur_t_axis))*ee];
        
%         n_trials = length(use_trials);
%     for tt = 1:n_trials
%         cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
%         cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
%         
%         cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
%         cur_ttill_end = trial_end_times(use_trials(tt)) - cur_t_axis;
%         
%         all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
%         all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
%         all_tsince_start = [all_tsince_start; cur_tsince_start];
%         all_ttill_end = [all_ttill_end; cur_ttill_end];
%         all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
%         all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
%         all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
%     end
%     trial_cnt = trial_cnt + n_trials;
    if strcmp(expt_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

trial_start_inds = [1; find(diff(all_trialvec) > 0)+1];

%% BIN SPIKES FOR MU AND SU
% % LOAD REFCLUSTERS
% fname = [cluster_dir '/final_cluster.mat'];
% if exist(fname,'file')
%     load(fname);
%     SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
%     for ii = 1:length(SU_numbers)
%         SU_tot_nblocks = sum(SU_ID_mat(:) == SU_numbers(ii));
%     end
%     fprintf('%d SUs Clustered\n',length(SU_numbers));
%     
% else
%     disp('No Cluster data found.');
% end

% %for SU probes
% all_binned_sua = nan(length(all_t_axis),length(SU_numbers));
% cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
% SU_ID_mat = SU_ID_mat(cur_used_blocks,:);
% [CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));
% 
% su_probes = nan(1,length(SU_numbers));
% all_su_spk_times = cell(length(SU_numbers),1);
% all_su_spk_inds = cell(length(SU_numbers),1);
% for ss = 1:length(SU_numbers)
%     used_clust_set = unique(CC(SU_ID_mat==SU_numbers(ss))); %set of clusters used to capture this SU
%     SU_block_probes(ss,:) = nan(1,length(cur_block_set));
%     cur_su_spk_times = [];
%     cur_su_spk_inds = [];
%     cur_blocks = [];
%     for cc = 1:length(used_clust_set)
%         cur_clust = used_clust_set(cc);
%         cur_probe = SU_clust_data(cur_clust).probe_num;
%         cur_clust_label = SU_clust_data(cur_clust).cluster_label;
%         cur_blocks = [cur_blocks; find(SU_ID_mat(:,cur_clust) == SU_numbers(ss))];
%         SU_block_probes(ss,cur_blocks) = cur_probe;
%         
%         all_su_inds = all_clust_ids{cur_probe} == cur_clust_label;
%         cur_su_spk_times = all_spk_times{cur_probe}(all_su_inds);
%         cur_su_spk_inds = all_spk_inds{cur_probe}(all_su_inds);
%         spk_block_inds = round(interp1(all_t_axis,all_blockvec,cur_su_spk_times));
%         cur_su_spk_times = cur_su_spk_times(ismember(spk_block_inds,cur_blocks));   
%         cur_su_spk_inds = cur_su_spk_inds(ismember(spk_block_inds,cur_blocks));
%         
%         all_su_spk_times{ss} = cat(1,all_su_spk_times{ss},cur_su_spk_times(:));
%         all_su_spk_inds{ss} = cat(1,all_su_spk_inds{ss},cur_su_spk_inds(:));
%     end
%     if ~isempty(all_su_spk_times{ss})
%     cur_suahist = histc(all_su_spk_times{ss},all_t_bin_edges);
%     cur_suahist(all_bin_edge_pts) = [];
%     cur_id_set = ismember(all_blockvec,cur_blocks);
%     all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
%     su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
%     end
% end

% double_spike_buffer = 3; %number of samples (in either direction) to exclude double spikes from adjacent-probe SUs
all_binned_mua = nan(length(all_t_axis),n_probes);
for cc = 1:n_probes
    %this is the set of blocks where this probe had an SU, and the
    %correspodning SU numbers
%     cur_set = find(SU_block_probes == cc);
%     if ~isempty(cur_set)
%         [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
%     else
%         cur_SS = [];
%     end
%     unique_su_nums = unique(cur_SS);
%     cur_mua_inds = find(all_clust_ids{cc} >= 0);
    cur_mua_inds = 1:length(all_spk_times{cc});
    
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

% all_sua_rate = nan(size(all_binned_sua));
% sua_block_mean_rates = nan(length(cur_block_set),length(SU_numbers));
% sua_block_n_spikes = nan(length(cur_block_set),length(SU_numbers));
% for ee = 1:length(cur_block_set)
%     cur_block_inds = find(all_blockvec==ee);
%     if ~isempty(cur_block_inds)
%         for ss = 1:length(SU_numbers)
%             all_sua_rate(cur_block_inds,ss) = jmm_smooth_1d_cor(all_binned_sua(cur_block_inds,ss),sua_sm_sig);
%         end
%         sua_block_mean_rates(ee,:) = mean(all_sua_rate(cur_block_inds,:));
%         sua_block_n_spikes(ee,:) = sum(all_binned_sua(cur_block_inds,:));
%     end
% end

%normalized firing rates (smoothed)
% all_sua_rate_norm = nan(size(all_sua_rate));
% all_mua_rate_norm = nan(size(all_mua_rate));
% for ee = 1:length(cur_block_set)
%     cur_block_inds = find(all_blockvec==ee);
%     if ~isempty(cur_block_inds)
% %         all_sua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_sua_rate(cur_block_inds,:),...
% %             sua_block_mean_rates(ee,:));
%         all_mua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_mua_rate(cur_block_inds,:),...
%             mua_block_mean_rates(ee,:));
%     end
% end
all_mua_rate_norm = bsxfun(@rdivide,all_mua_rate,mean(all_mua_rate));

%%
% niqf = 1/2/dt;
% lcf = 0.1;
% [b,a] = butter(2,[lcf/niqf],'high');
% 
% for ii = 1:n_probes
%     all_mua_rate_norm(:,ii) = filtfilt(b,a,all_mua_rate_norm(:,ii));
% end

%% DEFINE DATA USED FOR ANALYSIS
% used_inds = find(all_tsince_start >= beg_buffer & all_ttill_end >= end_buffer);
used_inds = (1:length(all_t_axis))';

%% PROCESS EYE TRACKING DATA
if strcmp(expt_type,'LP')
    expt_bar_ori = ones(size(expt_names))*bar_ori;
end
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,use_coils);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data_v2(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,zeros(length(cur_block_set)),used_inds);

[saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params);
orig_saccades = saccades;

% par_thresh = 10;
% orth_thresh = 10;
% [out_of_range] = detect_bad_fixation_v2(corrected_eye_vals_interp,all_trialvec,used_inds,par_thresh,orth_thresh,et_params.use_coils);
% % [out_of_range] = detect_bad_fixation(corrected_eye_vals_interp,all_trialvec,used_inds,par_thresh,orth_thresh,et_params.use_coils);
% fract_out = length(out_of_range)/length(used_inds);
% fprintf('Eliminating %.4f of data out of window\n',fract_out);
% used_inds(ismember(used_inds,out_of_range)) = [];
% NT = length(used_inds);

out_of_range = (abs(corrected_eye_vals_interp(used_inds,1)) > 8 | abs(corrected_eye_vals_interp(used_inds,2) > 8));
fract_out = sum(out_of_range)/length(used_inds);
fprintf('Eliminating %.4f of data out of window\n',fract_out);

%%
dsf = 2; %lfps originally sampled at 400Hz
use_lfps = 1:2:96;
lcf = 0;

all_V = [];
all_V_taxis = [];
for bb = cur_block_set
    %load lfps
    lfp_fname = sprintf('Expt%d_LFP.mat',bb);
    load(lfp_fname);
    Fs = lfp_params.Fsd;
    cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
    if lcf > 0
        [filt_b,filt_a] = butter(2,lcf/(Fs/2),'high');
        for ii = 1:length(use_lfps)
            cur_lfps(:,ii) = filtfilt(filt_b,filt_a,cur_lfps(:,ii));
        end
    end
    if dsf > 1
        [filt_b,filt_a] = butter(4,0.8/dsf,'low');
        for ii = 1:length(use_lfps)
            cur_lfps(:,ii) = filtfilt(filt_b,filt_a,cur_lfps(:,ii));
        end
        cur_lfps = downsample(cur_lfps,dsf);
        cur_lfp_t = downsample(lfp_t_ax,dsf);
    else
        cur_lfp_t = lfp_t_ax;
    end
    all_V = cat(1,all_V,cur_lfps);
    all_V_taxis = [all_V_taxis; cur_lfp_t'];
end

interp_V = interp1(all_V_taxis,all_V,all_t_axis);

%% PROCESS SACCADE STATS

%interpolate saccade start times and get rid of saccades that aren't within
%the t-axis binning
sac_start_times = [saccades(:).start_time];
sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
sac_start_inds(isnan(sac_start_inds)) = 1;
sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).stop_time]));
sac_error = abs(sac_start_times - all_t_axis(sac_start_inds)');
usac_set = find(ismember(sac_start_inds,used_inds) & ismember(sac_stop_inds,used_inds) & sac_error <= dt & sac_start_times >= all_t_axis(1) & [saccades(:).stop_time] <= all_t_axis(end));
saccades = saccades(usac_set);

sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).start_time]));
sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).stop_time]));
sac_peak_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).peak_time]));

% used_sac_set = find(ismember(sac_start_inds,used_inds) & ismember(sac_stop_inds,used_inds));

sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));

micro_thresh = 1;
max_sac_dur = 0.1;
sac_burst_isi = 0.15;
sacburst_set = find([saccades(:).isi] < sac_burst_isi | [saccades(:).next_isi] < sac_burst_isi);
micro_set = find([saccades(:).amplitude] < micro_thresh & [saccades(:).duration] < max_sac_dur);
micro_set(ismember(micro_set,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'
% sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);
sac_deltaX = [saccades(:).amplitude];
gsac_thresh = 1;
gsac_set = find(abs(sac_deltaX) > gsac_thresh & [saccades(:).duration] < max_sac_dur);

msac_dirs = [saccades(micro_set).direction];
large_msacs = micro_set([saccades(micro_set).amplitude] > 0.5);
small_msacs = micro_set([saccades(micro_set).amplitude] < 0.5);

% n_msac_bins = 8;
% msac_bin_edges = linspace(-pi,pi,n_msac_bins+1);
% msac_bin_cents = 0.5*msac_bin_edges(1:end-1) + 0.5*msac_bin_edges(2:end);

% dark_blocks = [4 5 7 8];
dark_blocks = [3 4 6 7];
light_blocks = setdiff(cur_block_set,dark_blocks);

light_sacs = find(ismember(all_blockvec(sac_start_inds),light_blocks));
dark_sacs = find(ismember(all_blockvec(sac_start_inds),dark_blocks));

light_msacs = intersect(micro_set,light_sacs);
light_gsacs = intersect(gsac_set,light_sacs);
dark_msacs = intersect(micro_set,dark_sacs);
dark_gsacs = intersect(gsac_set,dark_sacs);

n_amp_bins = 4;
usac_set = find([saccades(:).duration] < max_sac_dur);
sac_bin_edges = prctile([saccades(usac_set).amplitude],linspace(0,100,n_amp_bins+1));
sac_bin_cents = 0.5*sac_bin_edges(1:end-1) + 0.5*sac_bin_edges(2:end);

%% COMPUTE TRIG AVGS FOR MUA
nboot = [];

% %general averages
% [mua_data.msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(micro_set),backlag,forwardlag,nboot);
% [mua_data.gsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(gsac_set),backlag,forwardlag,nboot);

[light_msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(light_msacs),backlag,forwardlag,nboot);
[light_gsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(light_gsacs),backlag,forwardlag,nboot);
[dark_msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(dark_msacs),backlag,forwardlag,nboot);
[dark_gsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(dark_gsacs),backlag,forwardlag,nboot);

clear *sac_amp_avg 
for ii = 1:n_amp_bins
   cur_set = usac_set([saccades(usac_set).amplitude] >= sac_bin_edges(ii) & [saccades(usac_set).amplitude] < sac_bin_edges(ii+1)); 
    sac_amp_avg(ii,:,:) = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(cur_set),backlag,forwardlag,nboot);

   cur_sacs = cur_set(ismember(cur_set,light_sacs));
    light_sac_amp_avg(ii,:,:) = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(cur_sacs),backlag,forwardlag,nboot);
   cur_sacs = cur_set(ismember(cur_set,dark_sacs));
    dark_sac_amp_avg(ii,:,:) = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(cur_sacs),backlag,forwardlag,nboot);

end

[light_msac_LFP,lags] = get_event_trig_avg(interp_V,sac_start_inds(light_msacs),backlag,forwardlag,nboot);
[light_gsac_LFP,lags] = get_event_trig_avg(interp_V,sac_start_inds(light_gsacs),backlag,forwardlag,nboot);
[dark_msac_LFP,lags] = get_event_trig_avg(interp_V,sac_start_inds(dark_msacs),backlag,forwardlag,nboot);
[dark_gsac_LFP,lags] = get_event_trig_avg(interp_V,sac_start_inds(dark_gsacs),backlag,forwardlag,nboot);

%%
cmap = jet(n_amp_bins);
% figure; hold on
% for ii = 1:n_amp_bins
%     shadedErrorBar(lags*dt,mean(sac_amp_avg(ii,:,uset),3),std(sac_amp_avg(ii,:,uset),[],3)/sqrt(length(uset)),{'color',cmap(ii,:)})
%     
% end
figure; 
subplot(2,1,1)
hold on
for ii = 1:n_amp_bins
    shadedErrorBar(lags*dt,mean(light_sac_amp_avg(ii,:,uset),3),std(light_sac_amp_avg(ii,:,uset),[],3)/sqrt(length(uset)),{'color',cmap(ii,:)})
end
subplot(2,1,2)
hold on
for ii = 1:n_amp_bins
    shadedErrorBar(lags*dt,mean(dark_sac_amp_avg(ii,:,uset),3),std(dark_sac_amp_avg(ii,:,uset),[],3)/sqrt(length(uset)),{'color',cmap(ii,:)})
end

%%
%%
Fs = 1000;
dsf = 4;
Fsd = Fs/dsf;
niqf = Fs/2;
[bb,aa] = butter(2,[1 80]/niqf);

full_lfps = [];
full_lfp_taxis = [];
cur_toffset = 0;
for ee = 1:length(cur_block_set);
    % for ee = 1:3
    fprintf('Loading LFPs, Expt %d of %d\n',ee,length(cur_block_set));
    fname = sprintf('lemM%dA.%d.lfp.mat',Expt_num,cur_block_set(ee));
    load(fname);
    
    Fs = 1/LFP.Header.CRsamplerate;
    
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
        if length(cur_t_axis) > 50
        cur_LFP = [LFP.Trials(tt).LFP];
        cur_LFP = cur_LFP(cur_sp:end,:);
        cur_LFP = filtfilt(bb,aa,cur_LFP);
        
        cur_LFP = downsample(cur_LFP,dsf);
        cur_t_axis = downsample(cur_t_axis,dsf);
        
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
        expt_lfps = [expt_lfps; cur_LFP];
        %         plot(expt_lfp_t_axis)
        %         pause
        %         clf
        end
    end
    
    cur_uset = find(all_blockvec == ee);
    uinds = find(expt_lfp_t_axis >= all_t_axis(cur_uset(1)) & expt_lfp_t_axis <= all_t_axis(cur_uset(end)));
    full_lfps = cat(1,full_lfps,expt_lfps(uinds,:));
    full_lfp_taxis = cat(1,full_lfp_taxis,expt_lfp_t_axis(uinds));
    
    cur_toffset = trial_toffset(ee);
end
