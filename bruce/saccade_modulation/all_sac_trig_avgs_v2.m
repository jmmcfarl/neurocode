clear all
close all

Expt_name = 'G091';
% Expt_name = 'M232';

Expt_num = str2num(Expt_name(2:end));
data_dir = ['~/Data/bruce/' Expt_name];

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

TBT_expts = {'M275','M277'}; %list of expts where conditions are interleaved
is_TBT_expt = false;

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
if ~strcmp(Expt_name,'G081')
    load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
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
is_sim_msac_expt = false;
if strcmp(Expt_name,'G093')
    include_expts = {'rls.FaXwi','rls.FaXwiXimi'};
elseif any(strcmp(Expt_name,{'G095','M266','M270'}))
    include_expts = {'rls.Fa','rls.FaXimi','rls.FaXFaXFs'};
    is_sim_msac_expt = true;
elseif strcmp(Expt_name,'G081')
    include_expts = {'grating.OpXseRC','grating.OpRC'};
elseif any(strcmp(Expt_name,TBT_expts))
    include_expts = {'rls.AllSac','rls.imiXFa'};
    is_TBT_expt = true;
elseif strcmp(Expt_name,'G091')
    include_expts = {'rls.Fa','rls.orRC'};
else
    include_expts = {'rls.Fa', 'rls.FaXimi'};
end

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

if strcmp(expt_type,'UA')
    if strcmp(Expt_name,'G081')
        cur_block_set = find(included_type & expt_Fr == 1 & (expt_bar_ori == 0 | expt_bar_ori == 90 | expt_bar_ori == 45 | expt_bar_ori == 135));
    else
        cur_block_set = find(included_type & expt_Fr == 1 & (expt_bar_ori == 0 | expt_bar_ori == 90));
    end
else
    cur_block_set = find(included_type & expt_ce == 1 & expt_Fr == 1);
end
if strcmp(Expt_name,'G081')
    expt_has_ds = (expt_ijump==0)';
end
expt_has_ds(isnan(expt_has_ds)) = 0;

cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

sim_sac_expts = find(expt_has_ds(cur_block_set) ~= 1);
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');
or_expts = find(strcmp('rls.orRC',expt_names(cur_block_set)));
expt_bar_ori(cur_block_set(or_expts)) = 0;

poss_orth_expts = find(mod(expt_bar_ori(cur_block_set) - expt_sac_dir(cur_block_set),180) == 90);
if ~isempty(poss_orth_expts)
    fprintf('Warning, possible orthoganol saccade expts detected\n');
end

gsac_amp = unique(expt_sac_amp(cur_block_set([grayback_gs_expts; imback_gs_expts])));
if length(gsac_amp) > 1
    error('Multiple guided sac amps detected');
end
gsac_thresh = gsac_amp/3;

if strcmp(Expt_name,'G081')
    sim_sac_times = [0.7 1.4];
else
    sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
end

if is_TBT_expt
    sim_sac_expts = []; imback_gs_expts = []; grayback_gs_expts = [];
end

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
all_trial_wi = [];
all_trial_back = [];
all_trial_Ff = [];
all_trial_se = [];
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
    if ismember(ee,grayback_gs_expts)
        fprintf('Expt %s Block %d of %d; grayback GS, ori:%d\n',Expt_name,ee,length(cur_block_set),expt_bar_ori(cur_block));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Expt %s Block %d of %d; imback GS, ori:%d\n',Expt_name,ee,length(cur_block_set),expt_bar_ori(cur_block));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Expt %s Block %d of %d; SimSac, ori:%d\n',Expt_name,ee,length(cur_block_set),expt_bar_ori(cur_block));
    else
        fprintf('Expt %s Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_name,ee,length(cur_block_set));
    end
    
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

    if strcmp(Expt_name,'G093')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    end
    if isfield(Expts{cur_block}.Trials(1),'Fs')
        trial_Fs = [Expts{cur_block}.Trials(:).Fs];
    else
        trial_Fs = nan(1,length(trial_durs));
    end

    if is_TBT_expt
       trial_back = strcmp('image',{Expts{cur_block}.Trials(:).Bs});
       trial_back = trial_back(id_inds);
       
       trial_Ff = [Expts{cur_block}.Trials(:).Ff];
       trial_Ff = trial_Ff(id_inds);
       
       all_trial_back = cat(1,all_trial_back,trial_back(use_trials)');
       all_trial_Ff = cat(1,all_trial_Ff,trial_Ff(use_trials)');
    end
    
    if isfield(Expts{cur_block}.Trials(1),'se')
        trial_se = [Expts{cur_block}.Trials(:).se];
        trial_se = trial_se(id_inds);
    else
        trial_se = nan(size(id_inds));
    end
    all_trial_se = cat(1,all_trial_se,trial_se(use_trials)');
    
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

trial_start_inds = [1; find(diff(all_trialvec) > 0)+1];
simsac_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),sim_sac_expts));
grayback_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),grayback_gs_expts));
imback_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),imback_gs_expts));

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

%for G093 use only data where stripe width is 2 deg
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi == un_wi_vals(2));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end

%get rid of repeat trials
rpt_trials = find(all_trial_se == rpt_seed);
used_inds(ismember(all_trialvec(used_inds),rpt_trials)) = [];
%% PROCESS EYE TRACKING DATA
if strcmp(expt_type,'LP')
    expt_bar_ori = ones(size(expt_names))*bar_ori;
end
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,use_coils);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data_v2(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

[saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params);
orig_saccades = saccades;

par_thresh = 4;
orth_thresh = 1.2;
[out_of_range] = detect_bad_fixation_v2(corrected_eye_vals_interp,all_trialvec,used_inds,par_thresh,orth_thresh,et_params.use_coils);
% [out_of_range] = detect_bad_fixation(corrected_eye_vals_interp,all_trialvec,used_inds,par_thresh,orth_thresh,et_params.use_coils);
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
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);
gsac_set = find(abs(sac_deltaX) > gsac_thresh & [saccades(:).duration] < max_sac_dur);

outsacs = gsac_set(abs(sac_prepos(1,gsac_set)) < abs(sac_postpos(1,gsac_set)));
insacs = gsac_set(abs(sac_prepos(1,gsac_set)) > abs(sac_postpos(1,gsac_set)));
% negsacs = outsacs(sac_postpos(1,outsacs) < 0);
% possacs = outsacs(sac_postpos(1,outsacs) > 0);
negsacs = gsac_set(sac_postpos(1,gsac_set) < sac_prepos(1,gsac_set));
possacs = gsac_set(sac_postpos(1,gsac_set) > sac_prepos(1,gsac_set));
out_pos_sacs = intersect(outsacs,possacs);
out_neg_sacs = intersect(outsacs,negsacs);
in_pos_sacs = intersect(insacs,possacs);
in_neg_sacs = intersect(insacs,negsacs);

%compile indices of simulated saccades
all_sim_sacs = [];
if is_TBT_expt
    sim_sac_trials = find(all_trial_Ff == 0);
    sim_trial_inds = find(ismember(all_trialvec,sim_sac_trials));
    sim_sacs = cell(length(sim_sac_times),1);
    for ii = 1:length(sim_sac_times)
        sim_sacs{ii} = sim_trial_inds(all_tsince_start(sim_trial_inds(1:end-1)) < sim_sac_times(ii) & ...
            all_tsince_start(sim_trial_inds(2:end)) >= sim_sac_times(ii));
        all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
    end
else
    sim_expt_inds = find(ismember(all_blockvec,sim_sac_expts));
    sim_sacs = cell(length(sim_sac_times),1);
    for ii = 1:length(sim_sac_times)
        sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
            all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
        all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
    end
end
all_sim_sacs = sort(all_sim_sacs(ismember(all_sim_sacs,used_inds)));

if is_TBT_expt
    grayback_gs_trials = find(all_trial_back == 0);
    imback_gs_trials = find(all_trial_back  == 1);
    
    gback_sacs = find(ismember(all_trialvec(sac_start_inds),grayback_gs_trials));
    iback_sacs = find(ismember(all_trialvec(sac_start_inds),imback_gs_trials));
else
    gback_sacs = find(ismember(all_blockvec(sac_start_inds),grayback_gs_expts));
    iback_sacs = find(ismember(all_blockvec(sac_start_inds),imback_gs_expts));
    ss_sacs = find(ismember(all_blockvec(sac_start_inds),sim_sac_expts));
end
or_sacs = find(ismember(all_blockvec(sac_start_inds),or_expts));

sac_oris = expt_bar_ori(cur_block_set(all_blockvec(sac_start_inds)));
sim_sac_oris = expt_bar_ori(cur_block_set(all_blockvec(all_sim_sacs)));
sac_dds = expt_dd(cur_block_set(all_blockvec(sac_start_inds)));
sac_ce = expt_ce(cur_block_set(all_blockvec(sac_start_inds)));

gray_msac_set = intersect(gback_sacs,micro_set);
or_msac_set = intersect(micro_set,or_sacs);
msac_set = setdiff(micro_set,or_sacs);

msac_dirs = [saccades(micro_set).direction];
large_msacs = micro_set([saccades(micro_set).amplitude] > 0.5);
small_msacs = micro_set([saccades(micro_set).amplitude] < 0.5);

n_msac_bins = 8;
msac_bin_edges = linspace(-pi,pi,n_msac_bins+1);
msac_bin_cents = 0.5*msac_bin_edges(1:end-1) + 0.5*msac_bin_edges(2:end);

%% SACCADE TIMING ANALYSIS
sac_sm = round(0.025/dt);
binned_msacs = hist(sac_start_inds(micro_set),1:length(all_t_axis));
binned_gsacs = hist(sac_start_inds(gsac_set),1:length(all_t_axis));
binned_msac_sm = jmm_smooth_1d_cor(binned_msacs,sac_sm);
binned_gsac_sm = jmm_smooth_1d_cor(binned_gsacs,sac_sm);

%trial averages (USE ALL OF TRIAL HERE)
[gen_data.msac_simsac_trial_avg,trial_lags] = get_event_trig_avg(binned_msac_sm,simsac_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.msac_grayback_trial_avg,trial_lags] = get_event_trig_avg(binned_msac_sm,grayback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.msac_imback_trial_avg,trial_lags] = get_event_trig_avg(binned_msac_sm,imback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.gsac_simsac_trial_avg,trial_lags] = get_event_trig_avg(binned_gsac_sm,simsac_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.gsac_grayback_trial_avg,trial_lags] = get_event_trig_avg(binned_gsac_sm,grayback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
[gen_data.gsac_imback_trial_avg,trial_lags] = get_event_trig_avg(binned_gsac_sm,imback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);

maxlag = round(0.5/dt);
[gen_data.msac_acorr,acorr_lags] = xcov(binned_msac_sm,maxlag,'coeff');
[gen_data.gsac_acorr,acorr_lags] = xcov(binned_gsac_sm,maxlag,'coeff');
[gen_data.msac_gsac_xcorr,acorr_lags] = xcov(binned_gsac_sm,binned_msac_sm,maxlag,'coeff');

%% COMPUTE TRIG AVGS FOR MUA
nboot = [];
%set trial numbers to Inf so they don't get included in trig averaging
used_trialvec = ones(size(all_trialvec))*Inf;
used_trialvec(used_inds) = all_trialvec(used_inds);
clear mua_data

fprintf('Computing trig avgs for MUA\n',cc,n_probes);
%trial averages (USE ALL OF TRIAL HERE)
[mua_data.simsac_trial_avg,trial_lags] = get_event_trig_avg(all_mua_rate_norm,simsac_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0);
[mua_data.grayback_trial_avg,trial_lags] = get_event_trig_avg(all_mua_rate_norm,grayback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0);
[mua_data.imback_trial_avg,trial_lags] = get_event_trig_avg(all_mua_rate_norm,imback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0);

%general averages
[mua_data.msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.simsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,all_sim_sacs,backlag,forwardlag,nboot,used_trialvec,0);

[mua_data.msac_or_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(or_msac_set),backlag,forwardlag,nboot,used_trialvec,0);

[mua_data.large_msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(large_msacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.small_msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(small_msacs),backlag,forwardlag,nboot,used_trialvec,0);
for ii = 1:n_msac_bins
    cur_msac_set = micro_set(msac_dirs >= msac_bin_edges(ii) & msac_dirs < msac_bin_edges(ii+1));
    [mua_data.msac_dir_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(cur_msac_set),backlag,forwardlag,nboot,used_trialvec,0);
end

[mua_data.msac_end_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_stop_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_end_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_stop_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);

%background dependent
[mua_data.msac_gray_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(intersect(micro_set,gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_gray_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(intersect(gsac_set,gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.msac_im_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(intersect(micro_set,iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_im_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(intersect(gsac_set,iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);

%sac-location dependent
[mua_data.gsac_out_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(outsacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_in_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(insacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_neg_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(negsacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_pos_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(possacs),backlag,forwardlag,nboot,used_trialvec,0);

[mua_data.gsac_outpos_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(out_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_inpos_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(in_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_outneg_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(out_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_inneg_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(in_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);

%ce/dd dependent
[mua_data.msac_sparse_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(gray_msac_set(sac_dds(gray_msac_set) == 12)),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.msac_dense_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(gray_msac_set(sac_dds(gray_msac_set) == 67)),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.msac_monoc_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(gray_msac_set(sac_ce(gray_msac_set) == 1)),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.msac_binoc_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(gray_msac_set(sac_ce(gray_msac_set) == 0)),backlag,forwardlag,nboot,used_trialvec,0);

if Expt_name(1) == 'G'
    %saccade direction dependent
    poss_sac_oris = [0 90 45 135];
    mua_data.poss_sac_oris = poss_sac_oris;
    for ii = 1:length(poss_sac_oris)
        
        [mua_data.msac_dir_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(micro_set(sac_oris(micro_set) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data.gsac_dir_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(gsac_set(sac_oris(gsac_set) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data.simsac_dir_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate_norm,all_sim_sacs(sim_sac_oris == poss_sac_oris(ii)),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data.gsac_outdir_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(outsacs(sac_oris(outsacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data.gsac_indir_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(insacs(sac_oris(insacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        
        [mua_data.gsac_dir_outpos_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate,sac_start_inds(out_pos_sacs(sac_oris(out_pos_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data.gsac_dir_inpos_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate,sac_start_inds(in_pos_sacs(sac_oris(in_pos_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data.gsac_dir_outneg_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate,sac_start_inds(out_neg_sacs(sac_oris(out_neg_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data.gsac_dir_inneg_avg(ii,:,:),lags] = get_event_trig_avg(all_mua_rate,sac_start_inds(in_neg_sacs(sac_oris(in_neg_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        
        cur_stim_inds = used_inds(ismember(all_blockvec(used_inds),find(expt_bar_ori(cur_block_set) ==poss_sac_oris(ii))));
        mua_data.dir_avg_rates = mean(all_binned_mua(cur_stim_inds,:));
        mua_data.dir_std_rates = std(all_binned_mua(cur_stim_inds,:));
    end
end

mua_data.avg_rates = mean(all_binned_mua(used_inds,:));
mua_data.tot_nspikes = sum(all_binned_mua(used_inds,:));
    
%%
% close all
% for ss = 1:n_probes
%     subplot(2,1,1)
%     plot(lags*dt,mua_data.gsac_dir_outpos_avg(1,:,ss)/dt)
%     hold on
%     plot(lags*dt,mua_data.gsac_dir_outneg_avg(1,:,ss)/dt,'r')
%     plot(lags*dt,mua_data.gsac_dir_inpos_avg(1,:,ss)/dt,'k')
%     plot(lags*dt,mua_data.gsac_dir_inneg_avg(1,:,ss)/dt,'g')
%     
%     subplot(2,1,2)
%     plot(lags*dt,mua_data.gsac_dir_outpos_avg(2,:,ss)/dt)
%     hold on
%     plot(lags*dt,mua_data.gsac_dir_outneg_avg(2,:,ss)/dt,'r')
%     plot(lags*dt,mua_data.gsac_dir_inpos_avg(2,:,ss)/dt,'k')
%     plot(lags*dt,mua_data.gsac_dir_inneg_avg(2,:,ss)/dt,'g')
% 
%     pause
%     clf
% end
% 


%% COMPUTE TRIG AVGS FOR SUA
nboot = 100;
clear sua_data

%general averages
[sua_data.msac_avg,lags,sua_data.msac_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_avg,lags,sua_data.gsac_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.simsac_avg,lags,sua_data.simsac_std] = get_event_trig_avg_v2(all_sua_rate_norm,all_sim_sacs,backlag,forwardlag,nboot,used_trialvec,0);

[sua_data.msac_end_avg,lags,sua_data.msac_end_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_stop_inds(micro_set),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_end_avg,lags,sua_data.gsac_end_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_stop_inds(gsac_set),backlag,forwardlag,nboot,used_trialvec,0);

%background dependent
[sua_data.msac_gray_avg,lags,sua_data.msac_gray_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(intersect(micro_set,gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_gray_avg,lags,sua_data.gsac_gray_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(intersect(gsac_set,gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.msac_im_avg,lags,sua_data.msac_im_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(intersect(micro_set,iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_im_avg,lags,sua_data.gsac_im_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(intersect(gsac_set,iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);

%sac-location dependent
[sua_data.gsac_out_avg,lags,sua_data.gsac_out_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(outsacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_in_avg,lags,sua_data.gsac_in_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(insacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_neg_avg,lags,sua_data.gsac_neg_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(negsacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_pos_avg,lags,sua_data.gsac_pos_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(possacs),backlag,forwardlag,nboot,used_trialvec,0);

[sua_data.gsac_outpos_avg,lags,sua_data.gsac_outpos_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(out_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_inpos_avg,lags,sua_data.gsac_inpos_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(in_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_outneg_avg,lags,sua_data.gsac_outneg_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(out_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.gsac_inneg_avg,lags,sua_data.gsac_inneg_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(in_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);

%ce/dd dependent
[sua_data.msac_sparse_avg,~,sua_data.msac_sparse_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(gray_msac_set(sac_dds(gray_msac_set) == 12)),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.msac_dense_avg,lags,sua_data.msac_dense_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(gray_msac_set(sac_dds(gray_msac_set) == 67)),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.msac_monoc_avg,lags,sua_data.msac_monoc_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(gray_msac_set(sac_ce(gray_msac_set) == 1)),backlag,forwardlag,nboot,used_trialvec,0);
[sua_data.msac_binoc_avg,lags,sua_data.msac_binoc_std] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(gray_msac_set(sac_ce(gray_msac_set) == 0)),backlag,forwardlag,nboot,used_trialvec,0);

if Expt_name(1) == 'G'
    %saccade direction dependent
    poss_sac_oris = [0 90 45 135];
    sua_data.poss_sac_oris = poss_sac_oris;
    for ii = 1:length(poss_sac_oris)
        
        [sua_data.msac_dir_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(micro_set(sac_oris(micro_set) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [sua_data.gsac_dir_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(gsac_set(sac_oris(gsac_set) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [sua_data.simsac_dir_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate_norm,all_sim_sacs(sim_sac_oris == poss_sac_oris(ii)),backlag,forwardlag,nboot,used_trialvec,0);
        [sua_data.gsac_outdir_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(outsacs(sac_oris(outsacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [sua_data.gsac_indir_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(insacs(sac_oris(insacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
%         
%         [sua_data.gsac_dir_outpos_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(out_pos_sacs(sac_oris(out_pos_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
%         [sua_data.gsac_dir_inpos_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(in_pos_sacs(sac_oris(in_pos_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
%         [sua_data.gsac_dir_outneg_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(out_neg_sacs(sac_oris(out_neg_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
%         [sua_data.gsac_dir_inneg_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate_norm,sac_start_inds(in_neg_sacs(sac_oris(in_neg_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [sua_data.gsac_dir_outpos_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate,sac_start_inds(out_pos_sacs(sac_oris(out_pos_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [sua_data.gsac_dir_inpos_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate,sac_start_inds(in_pos_sacs(sac_oris(in_pos_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [sua_data.gsac_dir_outneg_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate,sac_start_inds(out_neg_sacs(sac_oris(out_neg_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        [sua_data.gsac_dir_inneg_avg(ii,:,:),lags] = get_event_trig_avg_v2(all_sua_rate,sac_start_inds(in_neg_sacs(sac_oris(in_neg_sacs) == poss_sac_oris(ii))),backlag,forwardlag,nboot,used_trialvec,0);
        
        cur_stim_inds = used_inds(ismember(all_blockvec(used_inds),find(expt_bar_ori(cur_block_set) ==poss_sac_oris(ii))));
        neg_inds = cur_stim_inds(corrected_eye_vals_interp(cur_stim_inds) < -gsac_thresh);
        pos_inds = cur_stim_inds(corrected_eye_vals_interp(cur_stim_inds) > gsac_thresh);
        cent_inds = cur_stim_inds(corrected_eye_vals_interp(cur_stim_inds) > gsac_thresh);
        sua_data.dir_avg_rates(:,ii) = nanmean(all_binned_sua(cur_stim_inds,:));
        sua_data.dir_std_rates(:,ii) = nanstd(all_binned_sua(cur_stim_inds,:));
    end
end

sua_data.avg_rates = nanmean(all_binned_sua(used_inds,:));
sua_data.tot_nspikes = nansum(all_binned_sua(used_inds,:));
sua_data.n_used_blocks = sum(~isnan(SU_block_probes),2);

%%
cd(save_dir)
sname = 'sac_trig_avg_data2';
save(sname,'sua_data','mua_data','lags','dt');

%%
close all
for ss = 1:length(sua_data.avg_rates)
    subplot(2,1,1)
    plot(lags*dt,sua_data.gsac_dir_outpos_avg(1,:,ss)/dt)
    hold on
    plot(lags*dt,sua_data.gsac_dir_outneg_avg(1,:,ss)/dt,'r')
    plot(lags*dt,sua_data.gsac_dir_inpos_avg(1,:,ss)/dt,'k')
    plot(lags*dt,sua_data.gsac_dir_inneg_avg(1,:,ss)/dt,'g')
    
    subplot(2,1,2)
    plot(lags*dt,sua_data.gsac_dir_outpos_avg(2,:,ss)/dt)
    hold on
    plot(lags*dt,sua_data.gsac_dir_outneg_avg(2,:,ss)/dt,'r')
    plot(lags*dt,sua_data.gsac_dir_inpos_avg(2,:,ss)/dt,'k')
    plot(lags*dt,sua_data.gsac_dir_inneg_avg(2,:,ss)/dt,'g')

    pause
    clf
end

%%
figure;
for ss = 1:length(su_probes)
   shadedErrorBar(lags*dt,sua_data.gsac_outpos_avg(:,ss),sua_data.gsac_outpos_std(:,ss)); hold on
    hold on
   shadedErrorBar(lags*dt,sua_data.gsac_inpos_avg(:,ss),sua_data.gsac_inpos_std(:,ss),{'color','b'}); hold on
   shadedErrorBar(lags*dt,sua_data.gsac_outneg_avg(:,ss),sua_data.gsac_outneg_std(:,ss),{'color','r'}); hold on
   shadedErrorBar(lags*dt,sua_data.gsac_inneg_avg(:,ss),sua_data.gsac_inneg_std(:,ss),{'color','g'}); hold on
    pause
    clf
    
end

%%
[mua_data.gsac_outpos_avg,lags,mua_data.gsac_outpos_std] = get_event_trig_avg_v2(all_mua_rate_norm,sac_start_inds(out_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_inpos_avg,lags,mua_data.gsac_inpos_std] = get_event_trig_avg_v2(all_mua_rate_norm,sac_start_inds(in_pos_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_outneg_avg,lags,mua_data.gsac_outneg_std] = get_event_trig_avg_v2(all_mua_rate_norm,sac_start_inds(out_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);
[mua_data.gsac_inneg_avg,lags,mua_data.gsac_inneg_std] = get_event_trig_avg_v2(all_mua_rate_norm,sac_start_inds(in_neg_sacs),backlag,forwardlag,nboot,used_trialvec,0);

%%
figure;
for ss = 1:n_probes
   shadedErrorBar(lags*dt,mua_data.gsac_outpos_avg(:,ss),mua_data.gsac_outpos_std(:,ss)); hold on
    hold on
   shadedErrorBar(lags*dt,mua_data.gsac_inpos_avg(:,ss),mua_data.gsac_inpos_std(:,ss),{'color','b'}); hold on
   shadedErrorBar(lags*dt,mua_data.gsac_outneg_avg(:,ss),mua_data.gsac_outneg_std(:,ss),{'color','r'}); hold on
   shadedErrorBar(lags*dt,mua_data.gsac_inneg_avg(:,ss),mua_data.gsac_inneg_std(:,ss),{'color','g'}); hold on
    pause
    clf
    
end

close