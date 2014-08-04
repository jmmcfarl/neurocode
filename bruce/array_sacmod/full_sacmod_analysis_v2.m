clear all
close all

% Expt_name = 'G085';
Expt_name = 'M275';
full_data_dir = ['/media/NTlab_data1/Data/bruce/' Expt_name];
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
else
    error('Unrecognized experiment name');
end

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.005;

min_trial_dur = 1;
beg_buffer = 0.25;
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
if strcmp(Expt_name,'G093')
    include_expts = {'rls.FaXwi','rls.FaXwiXimi'};
elseif any(strcmp(Expt_name,{'G095','M266','M270'}))
    include_expts = {'rls.Fa','rls.FaXimi','rls.FaXFaXFs'};
elseif strcmp(Expt_name,'G081')
    include_expts = {'grating.OpXseRC','grating.OpRC'};
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
    cur_block_set = find(included_type & expt_ce == 1 & expt_Fr == 1 & (expt_bar_ori == 0 | expt_bar_ori == 90 | expt_bar_ori == 45 | expt_bar_ori == 135));
    else
    cur_block_set = find(included_type & expt_ce == 1 & expt_Fr == 1 & (expt_bar_ori == 0 | expt_bar_ori == 90));
    end
else
    cur_block_set = find(included_type & expt_ce == 1 & expt_Fr == 1);
end
if strcmp(Expt_name,'M270')
    cur_block_set(cur_block_set == 5) = [];
end
if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
end
if strcmp(Expt_name,'G093')
    cur_block_set(cur_block_set ==  28) = []; %only 6 trials and causes problems
end

if strcmp(Expt_name,'G081')
    expt_has_ds = (expt_ijump==0)';
end
expt_has_ds(isnan(expt_has_ds)) = 0;

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

hori_expts = find(expt_bar_ori(cur_block_set) == 0);
ver_expts = find(expt_bar_ori(cur_block_set) == 90);
 
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
all_trial_start_times = [];
all_trial_end_times = [];
all_trial_blocknums = [];
all_bin_edge_pts = [];
trial_toffset = zeros(length(cur_block_set),1);
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
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        cur_ttill_end = trial_end_times(use_trials(tt)) - cur_t_axis;
        
        all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
        all_tsince_start = [all_tsince_start; cur_ttill_end];
        all_ttill_end = [all_ttill_end; cur_tsince_start];
        all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    end
    trial_cnt = trial_cnt + n_trials;
    if strcmp(expt_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
    end
end

trial_start_inds = [1; find(diff(all_trialvec) > 0)+1];
simsac_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),sim_sac_expts));
grayback_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),grayback_gs_expts));
imback_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),imback_gs_expts));

%% LOAD IN SPIKING DATA
fprintf('Loading spike data\n');
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
for ee = 1:length(cur_block_set);
    cur_block = cur_block_set(ee);
    if ee > 1
        cur_toffset = trial_toffset(ee-1);
    else
        cur_toffset = 0;
    end
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end
end

%% BIN SPIKING DATA
fprintf('Binning spike data\n');

su_used_blocks = false(length(cur_block_set),length(SU_numbers));
%for MU
all_binned_mua = nan(length(all_t_axis),n_probes);
for cc = 1:n_probes
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc},all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end

%for SU probes
fprintf('Using %d SUs\n',length(SU_numbers));
all_binned_sua = nan(length(all_t_axis),length(SU_numbers));
cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
SU_ID_mat = SU_ID_mat(cur_used_blocks,:); %matrix telling which cluster in each block is assigned a given SU number
[CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));
for ss = 1:length(SU_numbers)
    used_clust_set = unique(CC(SU_ID_mat==ss)); %set of clusters used to capture this SU
    SU_block_probes(ss,:) = nan(1,length(cur_block_set));
    for cc = 1:length(used_clust_set)
        cur_clust = used_clust_set(cc);
        cur_probe = SU_clust_data(cur_clust).probe_num;
        cur_clust_label = SU_clust_data(cur_clust).cluster_label;
        cur_blocks = find(SU_ID_mat(:,cur_clust) == ss);
        SU_block_probes(ss,cur_blocks) = cur_probe;
        
        all_su_inds = all_clust_ids{cur_probe} == cur_clust_label;
        all_su_spk_times = all_spk_times{cur_probe}(all_su_inds);
        spk_block_inds = round(interp1(all_t_axis,all_blockvec,all_su_spk_times));
        all_su_spk_times = all_su_spk_times(ismember(spk_block_inds,cur_blocks));
        
        cur_suahist = histc(all_su_spk_times,all_t_bin_edges);
        cur_suahist(all_bin_edge_pts) = [];
        cur_id_set = ismember(all_blockvec,cur_blocks);
        all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
    end
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

%% PROCESS EYE TRACKING DATA
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

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

%compile indices of simulated saccades
all_sim_sacs = [];
sim_expt_inds = find(ismember(all_blockvec,sim_sac_expts));
sim_sacs = cell(length(sim_sac_times),1);
for ii = 1:length(sim_sac_times)
    sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
        all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
    all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
end

%% PICK OUT SACCADES FOR ANALYSIS

%define which saccades to use
used_msacs = find(is_micro & ~is_blink & ismember(sac_start_inds,used_inds));
used_gsacs = find(is_gsac' & ~is_blink & ismember(sac_start_inds,used_inds));
used_simsacs = find(ismember(all_sim_sacs,used_inds));

used_gback_sacs = find(ismember(all_blockvec(sac_start_inds),grayback_gs_expts));
used_iback_sacs = find(ismember(all_blockvec(sac_start_inds),imback_gs_expts));

sac_oris = expt_bar_ori(cur_block_set(all_blockvec(sac_start_inds)));

sim_sac_oris = expt_bar_ori(cur_block_set(all_blockvec(all_sim_sacs)));

%microsacs excluding bursts
non_burst_msacs = used_msacs(~is_sacburst(used_msacs));
burst_msacs = used_msacs(is_sacburst(used_msacs));

hor_sacs = find(sac_oris == 0);
ver_sacs = find(sac_oris == 90);
hor_sim_sacs = find(sim_sac_oris == 0);
ver_sim_sacs = find(sim_sac_oris == 90);

gsac_to_left = hor_sacs(ismember(hor_sacs,used_gsacs(sac_postX(used_gsacs) < -gsac_amp/2)));
gsac_to_right = hor_sacs(ismember(hor_sacs,used_gsacs(sac_postX(used_gsacs) > gsac_amp/2)));
gsac_to_top = ver_sacs(ismember(ver_sacs,used_gsacs(sac_postX(used_gsacs) < -gsac_amp/2)));
gsac_to_bottom = ver_sacs(ismember(ver_sacs,used_gsacs(sac_postX(used_gsacs) > gsac_amp/2)));
msac_left = hor_sacs(ismember(hor_sacs,used_msacs(sac_postX(used_msacs) < -gsac_amp/2)));
msac_right = hor_sacs(ismember(hor_sacs,used_msacs(sac_postX(used_msacs) > gsac_amp/2)));
msac_top = ver_sacs(ismember(ver_sacs,used_msacs(sac_postX(used_msacs) < -gsac_amp/2)));
msac_bottom = ver_sacs(ismember(ver_sacs,used_msacs(sac_postX(used_msacs) > gsac_amp/2)));

%% SACCADE TIMING ANALYSIS
sac_sm = round(0.025/dt);
binned_msacs = hist(sac_start_inds(used_msacs),1:length(all_t_axis));
binned_gsacs = hist(sac_start_inds(used_gsacs),1:length(all_t_axis));
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
for cc = 1:n_probes
    
    fprintf('Computing trig avgs for MUA %d of %d\n',cc,n_probes);
    %trial averages (USE ALL OF TRIAL HERE)
    [mua_data(cc).simsac_trial_avg,trial_lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),simsac_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0);
    [mua_data(cc).grayback_trial_avg,trial_lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),grayback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0);
    [mua_data(cc).imback_trial_avg,trial_lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),imback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0);
    
    %general averages
    [mua_data(cc).msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).gsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).simsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),all_sim_sacs(used_simsacs),backlag,forwardlag,nboot,used_trialvec,0);
    
    [mua_data(cc).msac_end_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_stop_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).gsac_end_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_stop_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0);
    
    %background dependent
    [mua_data(cc).msac_gray_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(intersect(used_msacs,used_gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).gsac_gray_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(intersect(used_gsacs,used_gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).msac_im_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(intersect(used_msacs,used_iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).gsac_im_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(intersect(used_gsacs,used_iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0);

    %sac-location dependent
    [mua_data(cc).gsac_left_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(gsac_to_left),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).gsac_right_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(gsac_to_right),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).gsac_top_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(gsac_to_top),backlag,forwardlag,nboot,used_trialvec,0);
    [mua_data(cc).gsac_bottom_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(gsac_to_bottom),backlag,forwardlag,nboot,used_trialvec,0);

    if Expt_name(1) == 'G'
        %saccade direction dependent
        [mua_data(cc).msac_ver_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(intersect(used_msacs,ver_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data(cc).gsac_ver_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(intersect(used_gsacs,ver_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data(cc).simsac_ver_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),all_sim_sacs(intersect(used_simsacs,ver_sim_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data(cc).msac_hor_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(intersect(used_msacs,hor_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data(cc).gsac_hor_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),sac_start_inds(intersect(used_gsacs,hor_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
        [mua_data(cc).simsac_hor_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),all_sim_sacs(intersect(used_simsacs,hor_sim_sacs)),backlag,forwardlag,nboot,used_trialvec,0);
    end
    
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
        
        %trial averages (USE ALL OF TRIAL HERE)
        [sua_data(ss).simsac_trial_avg,trial_lags,sua_data(ss).simsac_trial_sem,sua_data(ss).nused_simsac_trials] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),simsac_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0,cur_inds);
        [sua_data(ss).grayback_trial_avg,trial_lags,sua_data(ss).grayback_trial_sem,sua_data(ss).nused_grayback_trials] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),grayback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0,cur_inds);
        [sua_data(ss).imback_trial_avg,trial_lags,sua_data(ss).imback_trial_sem,sua_data(ss).nused_imback_trials] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),imback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0,cur_inds);
        
        %general averages
        [sua_data(ss).msac_avg,lags,sua_data(ss).msac_sem,sua_data(ss).nused_msac] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).gsac_avg,lags,sua_data(ss).gsac_sem,sua_data(ss).nused_gsac] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).simsac_avg,lags,sua_data(ss).simsac_sem,sua_data(ss).nused_simsac] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),all_sim_sacs(used_simsacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        
        [sua_data(ss).msac_end_avg,lags,sua_data(ss).msac_end_sem,sua_data(ss).nused_msac] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_stop_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).gsac_end_avg,lags,sua_data(ss).gsac_end_sem,sua_data(ss).nused_gsac] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_stop_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);

        [sua_data(ss).gsac_left_avg,lags,sua_data(ss).gsac_left_sem] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(gsac_to_left),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).gsac_right_avg,lags,sua_data(ss).gsac_right_sem] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(gsac_to_right),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).gsac_top_avg,lags,sua_data(ss).gsac_top_sem] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(gsac_to_top),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).gsac_bottom_avg,lags,sua_data(ss).gsac_bottom_sem] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(gsac_to_bottom),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).msac_left_avg,lags,sua_data(ss).msac_left_sem] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(msac_left),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).msac_right_avg,lags,sua_data(ss).msac_right_sem] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(msac_right),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).msac_top_avg,lags,sua_data(ss).msac_top_sem] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(msac_top),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).msac_bottom_avg,lags,sua_data(ss).msac_bottom_sem] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(msac_bottom),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);

        %background dependent
        [sua_data(ss).msac_gray_avg,lags,sua_data(ss).msac_gray_sem,sua_data(ss).nused_msac_gray] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(intersect(used_msacs,used_gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).gsac_gray_avg,lags,sua_data(ss).gsac_gray_sem,sua_data(ss).nused_gsac_gray] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(intersect(used_gsacs,used_gback_sacs)),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).msac_im_avg,lags,sua_data(ss).msac_im_sem,sua_data(ss).nused_msac_im] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(intersect(used_msacs,used_iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        [sua_data(ss).gsac_im_avg,lags,sua_data(ss).gsac_im_sem,sua_data(ss).nused_gsac_im] = ...
            get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(intersect(used_gsacs,used_iback_sacs)),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
        
        %sassade direction dependent
        if Expt_name(1) == 'G'
            [sua_data(ss).msac_ver_avg,lags] = ...
                get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(intersect(used_msacs,ver_sacs)),backlag,forwardlag,[],used_trialvec,0,cur_use_inds);
            [sua_data(ss).gsac_ver_avg,lags] = ...
                get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(intersect(used_gsacs,ver_sacs)),backlag,forwardlag,[],used_trialvec,0,cur_use_inds);
            [sua_data(ss).simsac_ver_avg,lags] = ...
                get_event_trig_avg(all_sua_rate_norm(:,ss),all_sim_sacs(intersect(used_simsacs,ver_sim_sacs)),backlag,forwardlag,[],used_trialvec,0,cur_use_inds);
            [sua_data(ss).msac_hor_avg,lags] = ...
                get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(intersect(used_msacs,hor_sacs)),backlag,forwardlag,[],used_trialvec,0,cur_use_inds);
            [sua_data(ss).gsac_hor_avg,lags] = ...
                get_event_trig_avg(all_sua_rate_norm(:,ss),sac_start_inds(intersect(used_gsacs,hor_sacs)),backlag,forwardlag,[],used_trialvec,0,cur_use_inds);
            [sua_data(ss).simsac_hor_avg,lags] = ...
                get_event_trig_avg(all_sua_rate_norm(:,ss),all_sim_sacs(intersect(used_simsacs,hor_sim_sacs)),backlag,forwardlag,[],used_trialvec,0,cur_use_inds);
        end
    
        sua_data(ss).avg_rate = mean(all_binned_sua(cur_use_inds,ss));
        sua_data(ss).tot_nspikes = sum(all_binned_sua(cur_use_inds,ss));
        %start and stop trig avgs
        %micro sac burst and non-burst
    else
        sua_data(ss).used = false;
    end
end

%% FOR MODEL-BASED SAC-MOD ANALYSIS
% fprintf('Computing model-based MUA analysis\n');
% sac_bin_width = 1;
% sac_binspace = sac_bin_width*dt;
% sac_bin_edges = -(backlag*dt-sac_binspace/2):sac_binspace:(forwardlag*dt+sac_binspace/2);
% sac_bin_cents = 0.5*sac_bin_edges(1:end-1) + 0.5*sac_bin_edges(2:end);
% n_sac_bins = length(sac_bin_cents);
% 
% L2_params = create_L2_params([],[1 n_sac_bins],n_sac_bins);
% L2_params = create_L2_params(L2_params,n_sac_bins + [1 n_sac_bins],n_sac_bins);
% L2_params = create_L2_params(L2_params,2*n_sac_bins + [1 n_sac_bins],n_sac_bins);
% 
% gsac_inds = interp_sac_start_inds(used_gsacs);
% msac_inds = interp_sac_start_inds(used_msacs);
% trial_simsac_mat = zeros(length(all_t_axis),n_sac_bins);
% trial_gsac_mat = zeros(length(all_t_axis),n_sac_bins);
% trial_msac_mat = zeros(length(all_t_axis),n_sac_bins);
% for i = 1:n_sac_bins
%     for cc = 1:sac_bin_width
%         cur_inds = all_sim_sacs + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
%         cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
%         trial_simsac_mat(cur_inds,i) = 1;
%         cur_inds = gsac_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
%         cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
%         trial_gsac_mat(cur_inds,i) = 1;
%         cur_inds = msac_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
%         cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
%         trial_msac_mat(cur_inds,i) = 1;
%     end
% end
% 
% Xexpt = zeros(length(all_t_axis),length(cur_block_set)-1);
% for i = 1:length(cur_block_set)-1
%     cur_set = find(all_blockvec==i);
%     Xexpt(cur_set,i) = 1;
% end
% 
% cur_Xmat = [trial_gsac_mat trial_msac_mat trial_simsac_mat Xexpt];
% clear trial_gsac_mat trial_msac_mat trial_simsac_mat Xexpt
% 
% exclude_mua_probes = [su_probes 16]; %16 is a consistently bad probe
% include_mua_probes = setdiff(1:96,exclude_mua_probes);
% ov_mua_rate = mean(all_spike_rate_norm(:,include_mua_probes),2);
% [gen_data.mua_msac_avg,lags] = get_event_trig_avg(ov_mua_rate,interp_sac_start_inds(used_msacs),backlag,forwardlag,0,all_trialvec,0,used_inds);
% [gen_data.mua_gsac_avg,lags] = get_event_trig_avg(ov_mua_rate,interp_sac_start_inds(used_gsacs),backlag,forwardlag,0,all_trialvec,0,used_inds);
% [gen_data.mua_simsac_avg,lags] = get_event_trig_avg(ov_mua_rate,all_sim_sacs,backlag,forwardlag,0,all_trialvec,0,used_inds);
% [gen_data.mua_nb_msac,lags] = get_event_trig_avg(ov_mua_rate,interp_sac_start_inds(non_burst_msacs),backlag,forwardlag,0,all_trialvec,0,used_inds);
% [gen_data.mua_b_msac,lags] = get_event_trig_avg(ov_mua_rate,interp_sac_start_inds(burst_msacs),backlag,forwardlag,0,all_trialvec,0,used_inds);
% 
% Robs = sum(all_binned_spikes(used_inds,include_mua_probes),2);
% [fitmod] = regGLM_fit(cur_Xmat(used_inds,:),Robs,L2_params,ones(length(L2_params),1)*200,[],[],1);
% gen_data.mua_gsac_kern = fitmod.K((1:n_sac_bins));
% gen_data.mua_msac_kern = fitmod.K((1:n_sac_bins) + n_sac_bins);
% gen_data.mua_simsac_kern = fitmod.K((1:n_sac_bins) + 2*n_sac_bins);
% 
% clear cur_Xmat

