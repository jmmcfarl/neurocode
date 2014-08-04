clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

fig_dir = '/home/james/Analysis/bruce/saccade_modulation/ori_tuning2/';

%run on 86 [0, 90];
Expt_num = 296;
Expt_name = sprintf('M%.3d',Expt_num);
if Expt_num > 280
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end

cd(data_dir);

% load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
% load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

good_coils = [1 1];

has_data = cellfun(@(x) length(x),Expts) > 0;
expt_names(has_data) = cellfun(@(x) x.Header.expname,Expts(has_data),'UniformOutput',false);

cur_block_set = find(strcmp('rls.orXFaRC',expt_names));

if Expt_num == 277
    cur_block_set(cur_block_set == 26) = []; %different set of orientations used in this block
end

n_blocks = length(cur_block_set);
%%
stim_fs = 100; %in Hz
dt = 0.01;
Fr = 1;

min_trial_dur = 0.75;

beg_buffer = 0.15;
end_buffer = 0.05;
if Expt_num < 266
    trial_dur = 2;
else
    trial_dur = 4;
end
n_probes = 24;

use_right_eye = false;

n_use_blocks = Inf;

sac_backlag = round(0.15/dt);
sac_forlag = round(0.35/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_stim_or = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
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
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        cur_stim_or = Expts{cur_block}.Trials(use_trials(tt)).or;
        cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        n_frames = length(cur_stim_or);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times+cur_toffset < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times(:) + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis(:) + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_stim_or = [all_stim_or; cur_stim_or(:)];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
    trial_toffset(ee) = all_t_bin_edges(end);
    cur_toffset = trial_toffset(ee);
    cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
end
orig_stim_or = all_stim_or;
%%
poss_oris = unique(all_stim_or);
poss_oris(poss_oris < -90) = [];

n_coarse_bins = 9;
if Expt_num == 277
    bin_edges = linspace(-80,90,n_coarse_bins+1);
else
    bin_edges = linspace(0,180,n_coarse_bins+1);
end
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);

coarse_stim_or = all_stim_or;
for ii = 1:length(bin_cents)
    cur_set = find(all_stim_or >= bin_edges(ii) & all_stim_or < bin_edges(ii+1));
    coarse_stim_or(cur_set) = bin_cents(ii);
end
coarse_stim_or(all_stim_or >= bin_edges(end)) = bin_cents(end);

% all_stim_or = coarse_stim_or;
% poss_oris = bin_cents;

%add blank back into ori vec
poss_oris = cat(1,poss_oris(:),-1009);

full_nPix = length(poss_oris);
all_stim_mat = zeros(length(all_t_axis),full_nPix);
for ii = 1:full_nPix
    cur_set = find(all_stim_or == poss_oris(ii));
    all_stim_mat(cur_set,ii) = 1;
end
all_blank_vec = zeros(length(all_t_axis),1);
all_blank_vec(all_stim_or == -1009) = 1;
all_junk_vec = zeros(length(all_t_axis),1);
all_junk_vec(all_stim_or == -1005) = 1;

%% BIN SPIKES FOR MU AND SU
rec_type = 'LP';
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


%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_stim_times),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% PROCESS EYE TRACKING DATA
expt_bar_ori = zeros(size(cur_block_set));

% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,good_coils);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori,used_inds);

[saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params);

is_blink = detect_blinks(all_eye_ts,all_eye_vals,saccades,et_params);

[saccades,is_blink] = merge_blinks(saccades,is_blink);

sac_start_times = [saccades(:).start_time];
sac_stop_times = [saccades(:).stop_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
sac_stop_times(bad_sacs) = [];
saccades(bad_sacs) = [];
is_blink(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);

used_is_blink = is_blink(used_saccade_set);

sac_amps = [saccades(used_saccade_set).amplitude];
micro_sacs = find(sac_amps < 1 & ~used_is_blink');
big_sacs = find(sac_amps > 1.5 & ~used_is_blink');

sac_durs = [saccades(used_saccade_set).duration];
sac_prepos = reshape([saccades(used_saccade_set).pre_pos],[],length(used_saccade_set));
sac_postpos = reshape([saccades(used_saccade_set).post_pos],[],length(used_saccade_set));
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

sac_burst_isi = 0.15;
sacburst_set = find([saccades(used_saccade_set).isi] < sac_burst_isi | [saccades(used_saccade_set).next_isi] < sac_burst_isi);
micro_sacs(ismember(micro_sacs,sacburst_set)) = [];

%%
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

if Expt_num == 277 %problem with fixation point being in RFs with certain saccade directions
    trial_sacpos = nan(length(trial_start_inds),4);
    for ii = 1:length(trial_start_inds)
        cur_set = find(all_eye_ts >= all_trial_start_times(ii) & all_eye_ts <= all_trial_end_times(ii));
        if ~isempty(cur_set)
            trial_sacpos(ii,:) = nanmean(all_eye_vals(cur_set,:),1);
        end
    end
    trial_sac_angle = atan2(trial_sacpos(:,2),trial_sacpos(:,1))*180/pi;
    g1 = find(trial_sac_angle < -160 | trial_sac_angle > 160);
    g2 = find(trial_sac_angle > 100 & trial_sac_angle < 160);
    g3 = find(trial_sac_angle < 100 & trial_sac_angle > 40);
    g4 = find(trial_sac_angle > -160 & trial_sac_angle < -100);
    
    temp = all_tsince_start(used_inds(saccade_start_inds(big_sacs)));
    outward_sacs = big_sacs((temp >= 0.5 & temp < 1.3) | ...
        (temp >= 2.0 & temp < 2.6 | temp >= 3.5));
    inward_sacs = big_sacs((temp >= 1.3 & temp < 2.0) | ...
        (temp >= 2.6 & temp < 3.5));
    
    % exclude_ts = [g2 g4];
    exclude_ts = [g4];
    big_sacs = big_sacs(~ismember(all_trialvec(used_inds(saccade_start_inds(big_sacs))),exclude_ts));
end
%%
saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));
% saccade_trial_inds = all_trialvec(used_inds(saccade_stop_inds));

Xsac = zeros(NT,n_sac_bins);
Xmsac = zeros(NT,n_sac_bins);
for ii = 1:n_sac_bins
        cur_sac_target = saccade_start_inds(big_sacs) + sac_bincents(ii);
%     cur_sac_target = saccade_stop_inds(big_sacs) + sac_bincents(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(big_sacs(uu))) = [];
    Xsac(cur_sac_target,ii) = 1;
    
    cur_sac_target = saccade_start_inds(micro_sacs) + sac_bincents(ii);
%      cur_sac_target = saccade_stop_inds(micro_sacs) + sac_bincents(ii);
   uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(micro_sacs(uu))) = [];
    Xmsac(cur_sac_target,ii) = 1;
end



%%
all_norm_binned_mua = bsxfun(@rdivide,all_binned_mua,nanmean(all_binned_mua));
backlag = round(0.1/dt);
forwardlag = round(0.3/dt);
[all_sac_avgs,lags] = get_event_trig_avg(all_norm_binned_mua(used_inds,:),saccade_start_inds(big_sacs),backlag,forwardlag);
[all_msac_avgs,lags] = get_event_trig_avg(all_norm_binned_mua(used_inds,:),saccade_start_inds(micro_sacs),backlag,forwardlag);

blank_locs = find(all_stim_or(used_inds) == -1009);
[all_blank_avgs,lags] = get_event_trig_avg(all_norm_binned_mua(used_inds,:),blank_locs,backlag,forwardlag);

%%
stim_probs = nan(length(poss_oris),1);
for ii = 1:length(poss_oris)
    stim_probs(ii) = sum(all_stim_or(used_inds) == poss_oris(ii));
end
stim_probs = stim_probs/sum(stim_probs);
%%
blag = 2;
flag = 20;
ori_smwin = 1;

all_SU_tavgs = nan(length(su_probes),length(poss_oris),blag+flag+1);
for ss = 1:length(su_probes)
    fprintf('SU %d of %d\n',ss,length(su_probes));
    cur_inds = used_inds(~isnan(all_binned_sua(used_inds,ss)));
    cur_Robs = all_binned_sua(cur_inds,ss);
    if ~isempty(cur_Robs)
        su_used(ss) = true;
        for ii = 1:length(poss_oris)
            [all_SU_tavgs(ss,ii,:),ori_lags] = get_event_trig_avg_v3(cur_Robs,find(all_stim_or(cur_inds) == poss_oris(ii)),blag,flag);
        end
        
        for ii = 1:length(ori_lags)
            all_SU_tavgs(ss,1:end-1,ii) = jmm_smooth_circbc(all_SU_tavgs(ss,1:end-1,ii),ori_smwin);
        end
        
    else
        su_used(ss) = false;
    end
    su_avgrate(ss) = mean(cur_Robs);
end

%%
all_SU_Navgs = bsxfun(@minus,all_SU_tavgs,su_avgrate');
all_SU_NRavgs = bsxfun(@rdivide,all_SU_Navgs,std(all_SU_Navgs,[],3));
all_SU_NTavgs = bsxfun(@rdivide,all_SU_Navgs,std(all_SU_Navgs,[],2));
all_SU_NBavgs = bsxfun(@rdivide,all_SU_Navgs,all_SU_Navgs(:,end,:));
all_SU_powlags = squeeze(std(all_SU_Navgs,[],2));
[~,SU_bestlags] = max(all_SU_powlags,[],2);


SU_oritune_slice = nan(length(su_probes),length(poss_oris));
for ss = 1:length(su_probes)
    SU_oritune_slice(ss,:) = all_SU_tavgs(ss,:,SU_bestlags(ss));
end


backgnd = mean(all_SU_powlags(:,1:3),2);
for ss = 1:length(su_probes)
    if su_used(ss)
%         for ii = 1:length(poss_oris)
%             all_SU_NTavgs(ss,1:end-1,ii) = jmm_smooth_circbc(all_SU_NTavgs(ss,1:end-1,ii),ori_smwin);
%         end
        nogoodlags = find(all_SU_powlags(ss,:) <= 3*backgnd(ss));
        all_SU_NTavgs(ss,:,nogoodlags) = nan;
    end
end

all_su_tkern = nan(length(su_probes),length(ori_lags),2);
all_su_okern = nan(length(su_probes),length(poss_oris),2);
for ss = 1:length(su_probes)
    if su_used(ss)
        temp = squeeze(all_SU_Navgs(ss,:,:));
        [C,O] = princomp(temp);
        all_su_tkern(ss,:,1) = C(:,1);
        all_su_okern(ss,:,1) = O(:,1);
        all_su_tkern(ss,:,2) = C(:,2);
        all_su_okern(ss,:,2) = O(:,2);
    end
end


% SU_bestlags = SU_bestlags -1;

%% FOR MU
%%
all_MU_tavgs = nan(n_probes,length(poss_oris),blag+flag+1);
for ss = 1:n_probes
    fprintf('MU %d of %d\n',ss,n_probes);
    cur_inds = used_inds(~isnan(all_binned_mua(used_inds,ss)));
    cur_Robs = all_binned_mua(cur_inds,ss);
    if ~isempty(cur_Robs)
        for ii = 1:length(poss_oris)
            [all_MU_tavgs(ss,ii,:),ori_lags] = get_event_trig_avg_v3(cur_Robs,find(all_stim_or(cur_inds) == poss_oris(ii)),blag,flag);
        end
        
        for ii = 1:length(ori_lags)
            all_MU_tavgs(ss,1:end-1,ii) = jmm_smooth_circbc(all_MU_tavgs(ss,1:end-1,ii),ori_smwin);
        end
        
    end
    mu_avgrate(ss) = mean(cur_Robs);
end
all_MU_Navgs = bsxfun(@minus,all_MU_tavgs,mu_avgrate');
all_MU_NRavgs = bsxfun(@rdivide,all_MU_Navgs,std(all_MU_Navgs,[],3));
all_MU_NTavgs = bsxfun(@rdivide,all_MU_Navgs,std(all_MU_Navgs,[],2));
all_MU_powlags = squeeze(std(all_MU_Navgs,[],2));
[~,MU_bestlags] = max(all_MU_powlags,[],2);

MU_oritune_slice = nan(n_probes,length(poss_oris));
for ss = 1:n_probes
    MU_oritune_slice(ss,:) = all_MU_tavgs(ss,:,MU_bestlags(ss));
end

backgnd = mean(all_MU_powlags(:,1:3),2);
for ss = 1:n_probes
    %         for ii = 1:length(poss_oris)
    %             all_MU_NTavgs(ss,1:end-1,ii) = jmm_smooth_circbc(all_MU_NTavgs(ss,1:end-1,ii),ori_smwin);
    %         end
    nogoodlags = find(all_MU_powlags(ss,:) <= 3*backgnd(ss));
    all_MU_NTavgs(ss,:,nogoodlags) = nan;
end

all_mu_tkern = nan(length(su_probes),length(ori_lags),2);
all_mu_okern = nan(length(su_probes),length(poss_oris),2);
for ss = 1:24
    temp = squeeze(all_MU_Navgs(ss,:,:));
    [C,O] = princomp(temp);
    all_mu_tkern(ss,:,1) = C(:,1);
    all_mu_okern(ss,:,1) = O(:,1);
    all_mu_tkern(ss,:,2) = C(:,2);
    all_mu_okern(ss,:,2) = O(:,2);
end



%%
all_su_stainfo = nan(length(su_probes),length(sac_bincents));
all_su_ori_avg = nan(length(su_probes),length(sac_bincents),length(poss_oris));
all_su_sacrate = nan(length(su_probes),length(sac_bincents));
for ss = 1:length(su_probes)
    ss
    %%
    cur_Robs = all_binned_sua(used_inds,ss);
    cc_uinds = find(~isnan(cur_Robs));
    cur_inds = used_inds(cc_uinds);
    cur_Robs = all_binned_sua(cur_inds,ss);
    
    if ~isempty(cur_inds)
        
        cur_bestlag = ori_lags(SU_bestlags(ss));
        lagvec = 1:length(cc_uinds);
        lagvec = lagvec - cur_bestlag;
        
        bad_lags = find(lagvec < 1);
        lagvec(bad_lags) = 1;
        bad_lags = cat(1,bad_lags',find(all_trialvec(cur_inds(lagvec)) ~= all_trialvec(cur_inds)));
        
        lag_ori = all_stim_or(cur_inds(lagvec));
        lag_stim_mat = all_stim_mat(cur_inds(lagvec),:);
        lag_ori(bad_lags) = nan;
        lag_stim_mat(bad_lags,:) = nan;
        
        cur_Xsac = Xsac(cc_uinds,:);
%         cur_Xsac = Xmsac(cc_uinds,:);
        any_sac_inds = find(any(cur_Xsac > 0,2));
                %%
%         has_bar = find(lag_ori >= -90);
        has_bar = find(lag_ori >= -90 | lag_ori == -1009);
        ori_sta = nan(length(sac_bincents),length(poss_oris));
        ori_avg = nan(length(sac_bincents),length(poss_oris));
        sac_arate = nan(1,length(sac_bincents));
        for ii = 1:length(sac_bincents)
            cur_set = has_bar(cur_Xsac(has_bar,ii) == 1);
            sac_avg = sum(bsxfun(@times,lag_stim_mat(cur_set,:),cur_Robs(cur_set)))/sum(cur_Robs(cur_set));
            ov_avg = mean(lag_stim_mat(cur_set,:));
            ori_sta(ii,:) = sac_avg-ov_avg;
            
            ori_avg(ii,:) = mean(bsxfun(@times,lag_stim_mat(cur_set,:),cur_Robs(cur_set)))./mean(lag_stim_mat(cur_set,:));
            sac_arate(ii) = mean(cur_Robs(cur_set));
        end
        [~,mmloc] = minmax(sac_arate);
        mmloc = sort(mmloc);
        mmloc(1) = mmloc(1)-1;
                
        sta_ssinfo = nanmean(ori_avg.*log2(bsxfun(@rdivide,ori_avg,sac_arate')),2)./sac_arate';
        
        ov_avg = nansum(bsxfun(@times,lag_stim_mat,cur_Robs))/sum(cur_Robs);
        ov_sta = ov_avg - nanmean(lag_stim_mat);
        
        ori_navg = bsxfun(@rdivide,ori_avg,mean(ori_avg));        
        all_su_stainfo(ss,:) = sta_ssinfo;
        all_su_ori_avg(ss,:,:) = ori_avg;
        all_su_sacrate(ss,:) = sac_arate;
    end
end

sm_win = 0.5;
sm_su_ori_avg = all_su_ori_avg;
sm_su_sacrate = all_su_sacrate;
for ss = 1:length(su_probes)
    for ii = 1:length(poss_oris)
        sm_su_ori_avg(ss,:,ii) = jmm_smooth_1d_cor(sm_su_ori_avg(ss,:,ii),sm_win);
    end
    
    for ii = 1:length(sac_bincents)
       sm_su_ori_avg(ss,ii,1:end-1) = jmm_smooth_circbc(squeeze(sm_su_ori_avg(ss,ii,1:end-1)),ori_smwin); 
    end
    
    sm_su_sacrate(ss,:) = jmm_smooth_1d_cor(sm_su_sacrate(ss,:),sm_win);
end

sm_su_nori_avg = bsxfun(@rdivide,sm_su_ori_avg,reshape(SU_oritune_slice,length(su_probes),1,[]));
% sm_su_nori_avg = bsxfun(@rdivide,sm_su_ori_avg,sm_su_ori_avg(:,:,end));
% sm_su_nori_avg = bsxfun(@rdivide,sm_su_ori_avg,mean(sm_su_ori_avg,2));
sm_su_lori_avg = bsxfun(@rdivide,sm_su_ori_avg,mean(sm_su_ori_avg,3));

sm_su_info = nanmean(bsxfun(@times,sm_su_ori_avg.*log2(bsxfun(@rdivide,sm_su_ori_avg,sm_su_sacrate)),reshape(stim_probs,1,1,[])),3)./sm_su_sacrate;
sm_su_info_rate = bsxfun(@times,sm_su_info,sm_su_sacrate);

%%

for ss = 1:n_probes
    ss
    %%
    cur_Robs = all_binned_mua(used_inds,ss);
    cc_uinds = find(~isnan(cur_Robs));
    cur_inds = used_inds(cc_uinds);
    cur_Robs = all_binned_mua(cur_inds,ss);
    
    if ~isempty(cur_inds)
        
        cur_bestlag = ori_lags(MU_bestlags(ss));
        lagvec = 1:length(cc_uinds);
        lagvec = lagvec - cur_bestlag;
        
        bad_lags = find(lagvec < 1);
        lagvec(bad_lags) = 1;
        bad_lags = cat(1,bad_lags',find(all_trialvec(cur_inds(lagvec)) ~= all_trialvec(cur_inds)));
        
        lag_ori = all_stim_or(cur_inds(lagvec));
        lag_stim_mat = all_stim_mat(cur_inds(lagvec),:);
        lag_ori(bad_lags) = nan;
        lag_stim_mat(bad_lags,:) = nan;
        
        cur_Xsac = Xsac(cc_uinds,:);
%         cur_Xsac = Xmsac(cc_uinds,:);
        any_sac_inds = find(any(cur_Xsac > 0,2));
                %%
%         has_bar = find(lag_ori >= -90);
        has_bar = find(lag_ori >= -90 | lag_ori == -1009);
        ori_sta = nan(length(sac_bincents),length(poss_oris));
        ori_avg = nan(length(sac_bincents),length(poss_oris));
        sac_arate = nan(1,length(sac_bincents));
        for ii = 1:length(sac_bincents)
            cur_set = has_bar(cur_Xsac(has_bar,ii) == 1);
            sac_avg = sum(bsxfun(@times,lag_stim_mat(cur_set,:),cur_Robs(cur_set)))/sum(cur_Robs(cur_set));
            ov_avg = mean(lag_stim_mat(cur_set,:));
            ori_sta(ii,:) = sac_avg-ov_avg;
            
            ori_avg(ii,:) = mean(bsxfun(@times,lag_stim_mat(cur_set,:),cur_Robs(cur_set)))./mean(lag_stim_mat(cur_set,:));
            sac_arate(ii) = mean(cur_Robs(cur_set));
        end
        [~,mmloc] = minmax(sac_arate);
        mmloc = sort(mmloc);
        mmloc(1) = mmloc(1)-1;
                
        sta_ssinfo = nanmean(ori_avg.*log2(bsxfun(@rdivide,ori_avg,sac_arate')),2)./sac_arate';
        
        ov_avg = nansum(bsxfun(@times,lag_stim_mat,cur_Robs))/sum(cur_Robs);
        ov_sta = ov_avg - nanmean(lag_stim_mat);
        
%         ori_navg = bsxfun(@rdivide,ori_avg,mean(ori_avg));        
        all_mu_stainfo(ss,:) = sta_ssinfo;
        all_mu_ori_avg(ss,:,:) = ori_avg;
        all_mu_sacrate(ss,:) = sac_arate;
        
    end
end

sm_mu_ori_avg = all_mu_ori_avg;
sm_mu_sacrate = all_mu_sacrate;
for ss = 1:n_probes
    for ii = 1:length(poss_oris)
        sm_mu_ori_avg(ss,:,ii) = jmm_smooth_1d_cor(sm_mu_ori_avg(ss,:,ii),sm_win);
    end
    
    for ii = 1:length(sac_bincents)
        sm_mu_ori_avg(ss,ii,1:end-1) = jmm_smooth_circbc(squeeze(sm_mu_ori_avg(ss,ii,1:end-1)),ori_smwin);
    end

    sm_mu_sacrate(ss,:) = jmm_smooth_1d_cor(sm_mu_sacrate(ss,:),sm_win);
end

sm_mu_nori_avg = bsxfun(@rdivide,sm_mu_ori_avg,reshape(MU_oritune_slice,n_probes,1,[]));
% sm_mu_nori_avg = bsxfun(@rdivide,sm_mu_ori_avg,mean(sm_mu_ori_avg,2));
sm_mu_lori_avg = bsxfun(@rdivide,sm_mu_ori_avg,mean(sm_mu_ori_avg,3));

sm_mu_info = nanmean(bsxfun(@times,sm_mu_ori_avg.*log2(bsxfun(@rdivide,sm_mu_ori_avg,sm_mu_sacrate)),reshape(stim_probs,1,1,[])),3)./sm_mu_sacrate;
sm_mu_info_rate = bsxfun(@times,sm_mu_info,sm_mu_sacrate);


%%
xl = [-0.1 0.25];
close all
for ii = 1:length(su_probes)
    cur_inds = all_binned_sua(used_inds,ii);
    if su_used(ii)
    ii
    
%     f1 = figure();
    
    subplot(2,3,1)
    imagesc(ori_lags*dt,1:length(poss_oris),squeeze(all_SU_Navgs(ii,:,:)))
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    xlabel('Lag (s)');
    ylabel('Ori');
    cur_xl = xlim();
    line(cur_xl,[length(poss_oris) length(poss_oris)]-0.5,'color','w')
    title('Trig avgs');
    
    subplot(2,3,4)
    imagesc(ori_lags*dt,1:length(poss_oris),squeeze(all_SU_NRavgs(ii,:,:)))
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]*0.5);
    xlabel('Lag (s)');
    ylabel('Ori');
    cur_xl = xlim();
    line(cur_xl,[length(poss_oris) length(poss_oris)]-0.5,'color','w')
    title('Normalized trig avgs');

    subplot(2,3,2)
    imagesc(sac_bincents*dt,1:length(poss_oris),squeeze(sm_su_ori_avg(ii,:,:))');
    xlim(xl);
    line(xl,[length(poss_oris) length(poss_oris)]-0.5,'color','w')
    xlabel('Time since saccade (s)');
    ylabel('Ori');
    title('Avg rates');
    
    subplot(2,3,5)
    imagesc(sac_bincents*dt,1:length(poss_oris),squeeze(sm_su_nori_avg(ii,:,:))');
    xlim(xl);
    line(xl,[length(poss_oris) length(poss_oris)]-0.5,'color','w')
    xlabel('Time since saccade (s)');
    ylabel('Ori');
    title('Normalized avg rates');

    [~,mmloc] = minmax(sm_su_sacrate(ii,:));
    mmloc = sort(mmloc);
    mmloc(1) = mmloc(1) - 2;
    mmloc(2) = mmloc(1) + 4;
    
    cmap = jet(6);
    
    subplot(2,3,6)
    ax = plotyy(sac_bincents*dt,sm_su_info_rate(ii,:),sac_bincents*dt,sm_su_sacrate(ii,:)/dt);
    hold on
%     for jj = 1:6
%         plot(sac_bincents(jj+mmloc(1))*dt,sm_info_rate(ii,jj+mmloc(1)),'o','color',cmap(jj,:),'linewidth',2);
%     end
    xlim(ax(1),xl);
    xlim(ax(2),xl);
    xlabel('Time since saccade (s)');
    ylabel(ax(1),'Info rate (bits/sec)');
    ylabel(ax(2),'Rate (Hz)');
    title('Info rate');

    
%     hold on
%     for jj = 1:6
%         plot(squeeze(sm_su_lori_avg(ii,jj+mmloc(1),:)),'color',cmap(jj,:),'linewidth',2);
%     end
    
    subplot(2,3,3)
    ax = plotyy(sac_bincents*dt,sm_su_info(ii,:),sac_bincents*dt,sm_su_sacrate(ii,:)/dt);
    hold on
%     for jj = 1:6
%         plot(sac_bincents(jj+mmloc(1))*dt,sm_su_info(ii,jj+mmloc(1)),'o','color',cmap(jj,:),'linewidth',2);
%     end
    xlim(ax(1),xl);
    xlim(ax(2),xl);
    xlabel('Time since saccade (s)');
    ylabel(ax(1),'Info (bits/spk)');
    ylabel(ax(2),'Rate (Hz)');
    title('SS-info');
    
    pause
    clf
    
%     fig_width = 15;
%     rel_height = 0.8*2/3;
%     figufy(f1);
%     fname = [fig_dir sprintf('Gsac_orimod_E%d_SU%d.pdf',Expt_num,ii)];
% %      fname = [fig_dir sprintf('Msac_orimod_E%d_SU%d.pdf',Expt_num,ii)];
%    exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width);
%     fname
%     close(f1);

    end
end

%%
xl = [-0.05 0.2];
close all
for ii = 1:n_probes
    cur_inds = all_binned_mua(used_inds,ii);
    if ~isempty(cur_inds)
    ii
%     f1 = figure();
    subplot(2,3,1)
    imagesc(ori_lags*dt,1:length(poss_oris),squeeze(all_MU_Navgs(ii,:,:)))
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
     cur_xl = xlim();
    line(cur_xl,[length(poss_oris) length(poss_oris)]-0.5,'color','w')
   xlabel('Lag (s)');
    ylabel('Ori');
    title('Trig avgs');
    
    subplot(2,3,4)
    imagesc(ori_lags*dt,1:length(poss_oris),squeeze(all_MU_NRavgs(ii,:,:)))
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    cur_xl = xlim();
    line(cur_xl,[length(poss_oris) length(poss_oris)]-0.5,'color','w')
    xlabel('Lag (s)');
    ylabel('Ori');
    title('Normalized trig avgs');

    subplot(2,3,2)
    imagesc(sac_bincents*dt,1:length(poss_oris),squeeze(sm_mu_ori_avg(ii,:,:))');
    xlim(xl);
    line(xl,[length(poss_oris) length(poss_oris)]-0.5,'color','w')
    xlabel('Time since saccade (s)');
    ylabel('Ori');
    title('Avg rates');
    
    subplot(2,3,5)
    imagesc(sac_bincents*dt,1:length(poss_oris),squeeze(sm_mu_nori_avg(ii,:,:))');
    xlim(xl);
    line(xl,[length(poss_oris) length(poss_oris)]-0.5,'color','w')
    xlabel('Time since saccade (s)');
    ylabel('Ori');
    title('Normalized avg rates');

    [~,mmloc] = minmax(sm_mu_sacrate(ii,:));
    mmloc = sort(mmloc);
    mmloc(1) = mmloc(1) - 2;
    mmloc(2) = mmloc(1) + 4;
    
    cmap = jet(6);
    
    subplot(2,3,6)
    ax = plotyy(sac_bincents*dt,sm_mu_info_rate(ii,:),sac_bincents*dt,sm_mu_sacrate(ii,:)/dt);
    hold on
%     for jj = 1:6
%         plot(sac_bincents(jj+mmloc(1))*dt,sm_info_rate(ii,jj+mmloc(1)),'o','color',cmap(jj,:),'linewidth',2);
%     end
    xlim(ax(1),xl);
    xlim(ax(2),xl);
    xlabel('Time since saccade (s)');
    ylabel(ax(1),'Info rate (bits/sec)');
    ylabel(ax(2),'Rate (Hz)');
    title('Info rate');

    
%     hold on
%     for jj = 1:6
%         plot(squeeze(sm_mu_lori_avg(ii,jj+mmloc(1),:)),'color',cmap(jj,:),'linewidth',2);
%     end
    
    subplot(2,3,3)
    ax = plotyy(sac_bincents*dt,sm_mu_info(ii,:),sac_bincents*dt,sm_mu_sacrate(ii,:)/dt);
    hold on
%     for jj = 1:6
%         plot(sac_bincents(jj+mmloc(1))*dt,sm_mu_info(ii,jj+mmloc(1)),'o','color',cmap(jj,:),'linewidth',2);
%     end
    xlim(ax(1),xl);
    xlim(ax(2),xl);
    xlabel('Time since saccade (s)');
    ylabel(ax(1),'Info (bits/spk)');
    ylabel(ax(2),'Rate (Hz)');
    title('SS-info');
    
    pause
    clf

    
%     fig_width = 15;
%     rel_height = 0.8*2/3;
%     figufy(f1);
%     fname = [fig_dir sprintf('Gsac_orimod_E%d_MU%d.pdf',Expt_num,ii)];
% %     fname = [fig_dir sprintf('Msac_orimod_E%d_MU%d.pdf',Expt_num,ii)];
%     exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width);
%     close(f1);

    end
end
