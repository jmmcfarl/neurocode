clear all
close all

% Expt_name = 'G086';
Expt_name = 'M275';

Expt_num = str2num(Expt_name(2:end));
data_dir = ['~/Data/bruce/' Expt_name];

save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end

is_random_ori = false;

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

%dont fit stim models using these blocks
if Expt_num == 270
    ignore_blocks = [5];
elseif Expt_num == 275
    ignore_blocks = [10];
elseif Expt_num == 277
    ignore_blocks = [26];
else
   ignore_blocks = []; 
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
is_sim_msac_expt = false;
if is_random_ori
    include_expts = {'rls.orXFaRC'};
    if any(strcmp(Expt_name,TBT_expts))
        is_TBT_expt = true;
    end
else
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
    else
        include_expts = {'rls.Fa', 'rls.FaXimi'};
    end
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
cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

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
if is_random_ori
    gsac_amp = unique(expt_sac_amp(cur_block_set));
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
all_trial_Fa = [];
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
    
%         fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
%     load(fname,'Clusters');
%     for cc = 1:n_probes
%         all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
%         all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
%         all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
%     end

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
    if isfield(Expts{cur_block}.Trials(1),'Ff')
        trial_Ff = [Expts{cur_block}.Trials(:).Ff];
    else
        trial_Ff = nan(1,length(trial_durs));
    end
%     if isfield(Expts{cur_block}.Trials(1),'Fa')
%         trial_Fa = [Expts{cur_block}.Trials(:).Fa];
%     else
%         trial_Fa = nan(1,length(trial_durs));
%     end

    if is_TBT_expt
        if isfield(Expts{cur_block}.Trials(1),'Bs')
            trial_back = strcmp('image',{Expts{cur_block}.Trials(:).Bs});
            trial_back = trial_back(id_inds);
            all_trial_back = cat(1,all_trial_back,trial_back(use_trials)');
        end
        
        if isfield(Expts{cur_block}.Trials(1),'Ff')
        trial_Ff = [Expts{cur_block}.Trials(:).Ff];
        trial_Ff = trial_Ff(id_inds);
        all_trial_Ff = cat(1,all_trial_Ff,trial_Ff(use_trials)');
        end
        
%         if isfield(Expts{cur_block}.Trials(1),'Fa')
%             trial_Fa = [Expts{cur_block}.Trials(:).Fa];
%             trial_Fa = trial_Fa(id_inds);
%             all_trial_Fa = cat(1,all_trial_Fa,trial_Fa(use_trials)');
%         end
    end
    
    trial_se = [Expts{cur_block}.Trials(:).se];
    trial_se = trial_se(id_inds);
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
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);

%compute corrected eye data in bar-oriented frame
if is_random_ori
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,zeros(size(cur_block_set)),used_inds);    
else
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);
end

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);

if is_random_ori
par_thresh = 5;
orth_thresh = 5;
else
par_thresh = 4;
orth_thresh = 1.2;
end
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
if is_random_ori
    delta_sacpar = sqrt(sac_deltaX.^2 + sac_deltaY.^2);
else
    delta_sacpar = abs(sac_deltaX);
end
is_gsac = delta_sacpar' >= gsac_thresh;
sac_angle = atan(sac_postY./sac_postX);

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

%% PICK OUT SACCADES FOR ANALYSIS

%define which saccades to use
used_msacs = find(is_micro & ~is_blink & ismember(sac_start_inds,used_inds));
used_gsacs = find(is_gsac' & ~is_blink & ismember(sac_start_inds,used_inds));
used_simsacs = find(ismember(all_sim_sacs,used_inds));

if is_TBT_expt
    grayback_gs_trials = find(all_trial_back == 0);
    imback_gs_trials = find(all_trial_back  == 1);
    
    used_gback_sacs = find(ismember(all_trialvec(sac_start_inds),grayback_gs_trials));
    used_iback_sacs = find(ismember(all_trialvec(sac_start_inds),imback_gs_trials));
else
    used_gback_sacs = find(ismember(all_blockvec(sac_start_inds),grayback_gs_expts));
    used_iback_sacs = find(ismember(all_blockvec(sac_start_inds),imback_gs_expts));
end

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

%%
lfp_trial_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_trial_start_times));
lfp_trial_stop_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_trial_end_times));
lfp_sac_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),sac_start_times));
lfp_simsac_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_t_axis(all_sim_sacs)));
lfp_trial_inds = [lfp_trial_start_inds(:) lfp_trial_stop_inds(:)];
lfp_trial_start_inds(isnan(lfp_trial_start_inds)) = [];
lfp_trial_stop_inds(isnan(lfp_trial_stop_inds)) = [];

%% COMPUTE TRIG AVGS FOR LFPs
nboot = [];
%set trial numbers to Inf so they don't get included in trig averaging
used_trialvec = ones(size(all_trialvec))*Inf;
used_trialvec(used_inds) = all_trialvec(used_inds);
used_trialvec_interp = interp1(all_t_axis,used_trialvec,full_lfp_taxis);
trialvec_interp = round(interp1(all_t_axis,all_trialvec,full_lfp_taxis));
tsince_start_interp = interp1(all_t_axis,all_tsince_start,full_lfp_taxis);

csd_params.Fs = Fsd; %sample freq
csd_params.BrainBound = 1; %first channel that is in the brain
csd_params.ChanSep = 0.05; %channel sep in mm
csd_params.diam = 2; %current disc diameter (in mm)

[lfp_data.trial_onset_csd,lags] = get_event_trig_csd(full_lfps,lfp_trial_start_inds,backlag,forwardlag,csd_params);
[lfp_data.trial_offset_csd,lags] = get_event_trig_csd(full_lfps,lfp_trial_stop_inds,backlag,forwardlag,csd_params);

%general averages
[lfp_data.msac_avg,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec_interp,0);
[lfp_data.gsac_avg,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec_interp,0);
[lfp_data.simsac_avg,lags] = get_event_trig_avg(full_lfps,lfp_simsac_start_inds(used_simsacs),backlag,forwardlag,nboot,used_trialvec_interp,0);

[lfp_data.gsac_csd,lags] = get_event_trig_csd(full_lfps,lfp_sac_start_inds(used_gsacs),backlag,forwardlag,csd_params,used_trialvec_interp,0);
[lfp_data.msac_csd,lags] = get_event_trig_csd(full_lfps,lfp_sac_start_inds(used_msacs),backlag,forwardlag,csd_params,used_trialvec_interp,0);
[lfp_data.simsac_csd,lags] = get_event_trig_csd(full_lfps,lfp_simsac_start_inds,backlag,forwardlag,csd_params,used_trialvec_interp,0);

%%
% 
% trialvec_interp = round(interp1(all_t_axis,all_trialvec,full_lfp_taxis));
% sac_angle = atan2(sac_deltaY,sac_deltaX);
% sac_edges = [-pi -2.9 -2.1 -1.2 -0.25 0.25 1.2 2.1 2.9 pi];
% [b,I] = histc(sac_angle,sac_edges);
% I(I==9) = 1;
% 
% clear g_avg
% for ii = 1:8
%     [g_avg(ii,:,:),tlags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(used_gsacs(I(used_gsacs) == ii)),backlag,forwardlag,[],trialvec_interp,1);
% end
% 
% close all
% for ii = 1:8
% imagesc(tlags/Fsd,1:24,squeeze(g_avg(ii,:,:))')
% caxis([-0.5 0.5])
% pause
% clf
% end
%%
% % [tonset_avg,tlags,ev_std,n_events] = get_event_trig_avg(full_lfps,lfp_trial_start_inds,backlag,full_forwardlag,5,trialvec_interp,1);
% [tonset_avg,tlags] = get_event_trig_avg(full_lfps,lfp_trial_inds(:,1),backlag,full_forwardlag,[],trialvec_interp,1);
% [tonset_csd,tlags] = get_event_trig_csd(full_lfps,lfp_trial_inds(:,1),backlag,full_forwardlag,csd_params);
% 
% %%
% [g1_avg,tlags] = get_event_trig_avg(full_lfps,lfp_trial_inds(g1,1),backlag,full_forwardlag,[],trialvec_interp,1);
% [g2_avg,tlags] = get_event_trig_avg(full_lfps,lfp_trial_inds(g2,1),backlag,full_forwardlag,[],trialvec_interp,1);
% % [g3_avg,tlags] = get_event_trig_avg(full_lfps,lfp_trial_inds(g3,1),backlag,full_forwardlag,[],trialvec_interp,1);
% % [g4_avg,tlags] = get_event_trig_avg(full_lfps,lfp_trial_inds(g4,1),backlag,full_forwardlag,[],trialvec_interp,1);
% 
% % [g4_csd,tlags] = get_event_trig_csd(full_lfps,lfp_trial_inds(g4,1),backlag,full_forwardlag,csd_params);

%%
% cd(save_dir)
% sname = 'lfp_trig_avg_data';
% save(sname,'lfp_data','lags','dt');
