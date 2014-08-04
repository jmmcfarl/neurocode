% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G094';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

%%
dsf = 1; %lfps originally sampled at 400Hz
use_lfps = 1:2:96;
lcf = 0.5;

cur_expt_set = [35];

%%
fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_tsince_start = [];
all_exptvec = [];
all_trialvec = [];
all_trial_durs = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_V = [];
% all_trial_nph = [];
% all_trial_or = [];
% all_trial_tf = [];
all_trial_sf = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
%     trial_or = [Expts{cur_expt}.Trials(:).or];
%     trial_nph = [Expts{cur_expt}.Trials(:).nph];
%     trial_tf = [Expts{cur_expt}.Trials(:).tf];
    trial_sf = [Expts{cur_expt}.Trials(:).sf];
    
    [un_ids,id_inds] = unique(trial_ids);
    use_trials = id_inds(trial_durs(id_inds) >= 0.5);
    
    trial_start_times = trial_start_times(use_trials);
    trial_end_times = trial_end_times(use_trials);
    trial_durs = trial_durs(use_trials);
%     trial_or = trial_or(use_trials);
%     trial_tf = trial_tf(use_trials);
    trial_sf = trial_sf(use_trials);
%     trial_nph = trial_nph(use_trials);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times');
    all_trial_durs = cat(1,all_trial_durs,trial_durs');
%     all_trial_nph = cat(1,all_trial_nph,trial_nph');
%     all_trial_tf = cat(1,all_trial_tf,trial_tf');
    all_trial_sf = cat(1,all_trial_sf,trial_sf');
%     all_trial_or = cat(1,all_trial_or,trial_or');
    
    n_trials = length(use_trials);
    
    %load lfps
    lfp_fname = sprintf('Expt%d_LFP.mat',cur_expt);
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
    all_t_axis = [all_t_axis; cur_lfp_t'];
    
    expt_tsince_start = nan(length(cur_lfp_t),1);
    expt_exptvec = nan(length(cur_lfp_t),1);
    expt_trialvec = nan(length(cur_lfp_t),1);
    for tt = 1:n_trials
        cur_samples = find(cur_lfp_t >= trial_start_times(tt) & cur_lfp_t <= trial_end_times(tt));
        
        expt_tsince_start(cur_samples) = cur_lfp_t(cur_samples) - trial_start_times(tt);
        expt_exptvec(cur_samples) = ee;
        expt_trialvec(cur_samples) = tt + trial_cnt;
    end
    trial_cnt = trial_cnt + n_trials;
    
    all_exptvec = cat(1,all_exptvec,expt_exptvec);
    all_tsince_start = cat(1,all_tsince_start,expt_tsince_start);
    all_trialvec = cat(1,all_trialvec,expt_trialvec);
end
Fsd = Fs/dsf;

%% CLASSIFY TRIALS
start_buffer = round(Fsd*0.5);
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));

% all_drifting_grate_trials = find(all_trial_nph == 0);
% all_rp_grate_trials = find(all_trial_nph > 0);

%% TRIAL-AVERAGED POWER SPECTRA
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));

params.Fs = Fsd;
params.tapers = [3 5];
movingwin = [1 1];
sMarkers = [all_trial_start_inds(:)+start_buffer all_trial_end_inds(:)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];


clear S_rp_or S_dg_or
un_sfs = unique(all_trial_sf);
for ii = 1:length(un_sfs)
    fprintf('SF %d of %d\n',ii,length(un_sfs));
    cur_sf_trials = find(all_trial_sf == un_sfs(ii) & used_trial_durs > movingwin(1));
    uset = randperm(length(cur_sf_trials)); if length(uset) > 3; uset(4:end) = []; end;
    cur_sf_trials = cur_sf_trials(uset);
%     cur_dg_trials = all_drifting_grate_trials(ismember(all_drifting_grate_trials,cur_tf_trials));
%     cur_rp_trials = all_rp_grate_trials(ismember(all_rp_grate_trials,cur_tf_trials));
    
    for ll = 1:length(use_lfps)
        
        [S_or(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_sf_trials,:));
%         [S_rp_or(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_sp_trials,:));
        
    end
end