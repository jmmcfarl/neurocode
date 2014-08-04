clear all
% close all

addpath('~/James_scripts/bruce/G081/')

cd ~/Data/bruce/G081/
load jbeG081Expts.mat
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));

stim_fs = 100; %in Hz
use_units = 1:96;

Fs = 200;
new_dt = .005;

beg_buffer = round(Fs*0.2); %don't use the first X data after start of a trial.
end_buffer = round(Fs*0.2);

%% PARSE TRIAL DATA STRUCTURES
for i = 1:length(Expts)
    if strcmp(Expts{i}.Header.expname,'grating.OpXseRC') | strcmp(Expts{i}.Header.expname,'grating.OpRC')
        is_bar_expt(i) = 1;
    else
        is_bar_expt(i) = 0;
    end
    
    if strcmp(Expts{i}.Stimvals.Bs,'image')
        expt_image_back(i) = 1;
    else
        expt_image_back(i) = 0;
    end
    
    expt_sim_sacs(i) = Expts{i}.Stimvals.ijump;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
    
end
expt_bar_ori(expt_bar_ori == -45) = 135;

load ./all_un_bar_pos
n_bar_pos = size(all_un_bar_pos,1);


%% DETERMINE SET OF EXPERIMENTS TO USE
% flen = 12;
% bar_oris = [0];
% un_bar_pos = all_un_bar_pos(:,1);

% fprintf('Analyzing %d ori expts\n',bar_oris);

%     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE));

% %expts with X deg bars and gray back (sim sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris);

%expts with X deg bars and any back (including sim sacs)
cur_expt_set = find(is_bar_expt==1);

cur_expt_set(cur_expt_set<=8) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
% cur_expt_set(cur_expt_set >= 46 & cur_expt_set <= 51) = []; %problem with image background
cur_expt_set(ismember(cur_expt_set,[46 48 49 51])) = []; %problem with image background
cur_expt_set(cur_expt_set > 60) = []; %no rect

%% COMPUTE TRIAL DATA
all_stim_times = [];
all_phase = [];
all_Op = [];
all_bar_mat = [];

all_t_axis = [];
all_rel_stimes = [];
all_rel_etimes = [];
all_binned_spks = [];
all_used_inds = [];
all_exptvec = [];
all_trialvec = [];
all_se = [];
all_trial_start_times = [];
all_trial_end_times = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    n_trials = length(Expts{cur_expt}.Trials);
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
    all_trial_start_times = [all_trial_start_times; trial_start_times(:)];
    all_trial_end_times = [all_trial_end_times; trial_end_times(:)];
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(tt).Start/1e4;
        cur_Op = Expts{cur_expt}.Trials(tt).Op;
        cur_phase = Expts{cur_expt}.Trials(tt).ph;
        if isfield(Expts{cur_expt}.Trials(tt),'se')
            cur_se = Expts{cur_expt}.Trials(tt).se;
        else
            cur_se = nan;
        end
        cur_t_edges = [(Expts{cur_expt}.Trials(tt).Start(1)/1e4):1/Fs:(Expts{cur_expt}.Trials(tt).End(end)/1e4)];
        cur_t_cents = cur_t_edges(1:end-1) + 2/Fs;
        cur_binned_spks = nan(length(cur_t_cents),length(use_units));
        for cc = 1:length(use_units)
            cur_hist = histc(Clusters{use_units(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        cur_used_inds = ones(length(cur_t_cents),1);
        cur_used_inds(1:beg_buffer) = 0;
        cur_used_inds(end-end_buffer+1:end) = 0;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_t_axis = [all_t_axis; cur_t_cents(:)];
        all_rel_stimes = [all_rel_stimes; cur_t_cents(:)- trial_start_times(tt)];
        all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_t_cents(:)];
        all_Op = [all_Op; cur_Op];
        all_phase = [all_phase; cur_phase'];
        all_used_inds = [all_used_inds; cur_used_inds];
        %         all_bar_mat = [all_bar_mat; bar_Xmat];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(length(cur_t_cents),1)*ee];
        all_se = [all_se; ones(length(cur_t_cents),1)*cur_se];
        all_trialvec = [all_trialvec; ones(length(cur_t_cents),1)*tt];
    end
    
end



%%
lin_X = zeros(length(all_t_axis),length(cur_expt_set)-1);
for i = 1:length(cur_expt_set)-1
    cur_set = find(all_exptvec==i);
    lin_X(cur_set,i) = 1;
end

[trial_set,ia,ic] = unique([all_trialvec all_exptvec],'rows');
n_trials = length(ia);

% [trial_set,ia,ic] = unique(all_trialvec);
% n_trials = length(trial_set);

rp = randperm(n_trials);
rand_trial_vec = rp(ic);
[~,ind_shuff] = sort(rand_trial_vec);

buf_time = 0.15;

bar135_expts = find(expt_bar_ori(cur_expt_set)==135);

n_xv_folds = 10;
xv_frac = 0.1;
n_xv_trials = floor(n_trials*xv_frac);
rp_tset = randperm(n_trials);
all_tr_inds = [];
for ii = 1:n_xv_folds
    cur = (ii-1)*n_xv_trials + (1:n_xv_trials);
    cur_tset = rp_tset(cur);
    xv_inds{ii} = find(ismember(ic,cur_tset));
    xv_inds{ii}(all_rel_stimes(xv_inds{ii}) < buf_time | all_rel_etimes(xv_inds{ii}) < buf_time) = [];
%     xv_inds{ii}(all_rel_stimes(xv_inds{ii}) < buf_time) = [];
    xv_inds{ii}(ismember(all_exptvec(xv_inds{ii}),bar135_expts)) = [];
    all_tr_inds = [all_tr_inds; xv_inds{ii}];
end
all_tr_inds = sort(all_tr_inds);

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
load ./jbeG081.em.mat
all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
for ee = 1:length(cur_expt_set);
    cur_set = find(all_exptvec==ee);
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));

    eye_dt = median(diff(eye_ts_interp));
    eye_fs = 1/eye_dt;
    lEyeXY = eye_vals_interp(:,1:2);
    rEyeXY = eye_vals_interp(:,3:4);
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),3);
    sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;

    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);

    all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
    all_eye_speed = [all_eye_speed; eye_speed];
    all_eye_ts = [all_eye_ts; eye_ts_interp'];
end
back_pts = 1 + find(diff(all_eye_ts) <= 0);
double_samples = [];
for i = 1:length(back_pts)
    next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
    double_samples = [double_samples back_pts(i):next_forward];
end
all_eye_ts(double_samples) = [];
all_eye_speed(double_samples) = [];
all_eye_vals(double_samples,:) = [];
eye_fs = 1/median(diff(all_eye_ts));

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);


all_eye_intrial = zeros(size(all_eye_ts));
for i = 1:length(all_trial_start_times)
    curset = find(all_eye_ts >= all_trial_start_times(i) & all_eye_ts <= all_trial_end_times(i));
    all_eye_intrial(curset) = 1;
end

%%
%compute times saccades
sac_thresh = 10;
min_isi = 0.05;

orig_saccade_inds = 1 + find(all_eye_speed(1:end-1) < sac_thresh & all_eye_speed(2:end) > sac_thresh);
orig_saccade_inds = unique(orig_saccade_inds);
isis = [Inf; diff(orig_saccade_inds)]/eye_fs;

% double_sacs = find(isis < min_isi);
% orig_saccade_inds(double_sacs) = [];
isis(all_eye_intrial(orig_saccade_inds) == 0) = [];
orig_saccade_inds(all_eye_intrial(orig_saccade_inds) == 0) = [];
interp_saccade_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_eye_ts(orig_saccade_inds)));
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

peri_thresh = 3;
sac_starts = nan(size(interp_saccade_inds));
sac_stops = nan(size(interp_saccade_inds));
for i = 1:length(interp_saccade_inds)
    cur_up = find(interp_eye_speed(1:interp_saccade_inds(i)) < peri_thresh,1,'last');
    cur_down = find(interp_eye_speed(interp_saccade_inds(i):end) < peri_thresh,1,'first');
    if ~isempty(cur_up) & ~isempty(cur_down)
        sac_starts(i) = cur_up;
        sac_stops(i) = interp_saccade_inds(i) + cur_down;
    end
end
bad = find(isnan(sac_starts));
interp_saccade_inds(bad) = []; sac_starts(bad) = []; sac_stops(bad) = []; isis(bad) = [];

sac_starts = sac_starts - round(Fs*0.01);
sac_stops = sac_stops + round(Fs*0.01);

%%
all_msacs = [];
all_out_sacs = [];
all_ret_sacs = [];
bar_oris = [0 45 90 135];
for bb = 1:4
    
%     orth_eye_pos = interp_eye_vals(:,2);
    par_eye_pos = interp_eye_vals(:,1)*cos(bar_oris(bb)*pi/180) - interp_eye_vals(:,2)*sin(bar_oris(bb)*pi/180);
    bad_pts = find(abs(par_eye_pos) > 1); %cheap way of detecting blinks because the right eye signal was unreliable even for that at times

    cur_bar_expts = find(expt_bar_ori(cur_expt_set) == bar_oris(bb));
    cur_inds = find(ismember(all_exptvec,cur_bar_expts));
    cur_sac_set = find(ismember(interp_saccade_inds,cur_inds));

    pre_pos = par_eye_pos(sac_starts(cur_sac_set));
    post_pos = par_eye_pos(sac_stops(cur_sac_set));

    sac_dist = abs(post_pos-pre_pos);
    if bar_oris(bb) == 0
        outgoing_sacs = find(pre_pos > -0.9 & post_pos < -0.9);
        returning_sacs = find(pre_pos < -1 & post_pos > -1);
    elseif bar_oris(bb) == 135
        outgoing_sacs = find(pre_pos > -1 & post_pos < -1);
        returning_sacs = find(pre_pos < -1.1 & post_pos > -1.2);
    elseif bar_oris(bb) == 45
        outgoing_sacs = find(pre_pos > -1 & post_pos < -1);
        returning_sacs = find(pre_pos < -1.1 & post_pos > -1);
    elseif bar_oris(bb) == 90
        outgoing_sacs = find(pre_pos > -1 & post_pos < -1);
        returning_sacs = find(pre_pos < -1 & post_pos > -1);
    end
    msacs = find(sac_dist' < 1); msacs(ismember(msacs,outgoing_sacs)) = []; msacs(ismember(msacs,returning_sacs)) = [];
    
    
    all_out_sacs = [all_out_sacs; cur_sac_set(outgoing_sacs)];
    all_ret_sacs = [all_ret_sacs; cur_sac_set(returning_sacs)];
    all_msacs = [all_msacs; (cur_sac_set(msacs))];
    
    
end
all_msacs = sort(all_msacs);
all_ret_sacs = sort(all_ret_sacs);
all_out_sacs = sort(all_out_sacs);

%%
min_isi = 0.05;
all_out_sacs(isis(all_out_sacs) < min_isi) = [];
all_ret_sacs(isis(all_ret_sacs) < min_isi) = [];

bad = find(isis(all_msacs) < min_isi);
% bad = [bad; 1+ find(isis(all_msacs(2:end)-1) < min_isi)];
all_msacs(bad) = [];


max_isi = 0.2;
% bad = find(isis(all_msacs) > max_isi);
% bad = find(isis(all_msacs(1:end-1)+1) > max_isi);
% all_msacs(bad) = [];
short_isi_msacs = all_msacs(isis(all_msacs) <= max_isi);
long_isi_msacs = all_msacs(isis(all_msacs) > max_isi);

%%
all_out_sacs = interp_saccade_inds(all_out_sacs);
all_ret_sacs = interp_saccade_inds(all_ret_sacs);

all_msacs = interp_saccade_inds(all_msacs);
short_isi_msacs = interp_saccade_inds(short_isi_msacs);
long_isi_msacs = interp_saccade_inds(long_isi_msacs);

all_msacs(all_rel_stimes(all_msacs) < 0.2 | all_rel_stimes(all_msacs) > 1.8) = [];
short_isi_msacs(all_rel_stimes(short_isi_msacs) < 0.2 | all_rel_stimes(short_isi_msacs) > 1.8) = [];
long_isi_msacs(all_rel_stimes(long_isi_msacs) < 0.2 | all_rel_stimes(long_isi_msacs) > 1.8) = [];

use_out = find(all_rel_stimes(all_out_sacs) > 0.8 & all_rel_stimes(all_out_sacs) < 1.3);
use_ret = find(all_rel_stimes(all_ret_sacs) > 1.4 & all_rel_stimes(all_ret_sacs) < 1.8);
all_out_sacs = all_out_sacs(use_out);
all_ret_sacs = all_ret_sacs(use_ret);

all_on_stim = 1+find(all_rel_stimes(1:end-1) <= 0.7 & all_rel_stimes(2:end) > 0.7);
all_off_stim = 1+find(all_rel_stimes(1:end-1) <= 1.4 & all_rel_stimes(2:end) > 1.4);
sim_expt_set = find(expt_sim_sacs(cur_expt_set) > 0);
all_on_stim(~ismember(all_exptvec(all_on_stim),sim_expt_set)) = [];
all_off_stim(~ismember(all_exptvec(all_off_stim),sim_expt_set)) = [];
all_sim_sacs = sort([all_on_stim; all_off_stim]);

% all_don_stim = 1+find(all_rel_stimes(1:end-1) <= 0.7 & all_rel_stimes(2:end) > 0.7);
% all_doff_stim = 1+find(all_rel_stimes(1:end-1) <= 1.4 & all_rel_stimes(2:end) > 1.4);
% nsim_expt_set = find(expt_sim_sacs(cur_expt_set) == 0);
% all_don_stim(~ismember(all_exptvec(all_don_stim),nsim_expt_set)) = [];
% all_doff_stim(~ismember(all_exptvec(all_doff_stim),nsim_expt_set)) = [];
% all_sim_sacs = sort([all_don_stim all_doff_stim]);

dtimes = [-Inf; diff(all_rel_stimes)];
gray_expts = find(expt_image_back(cur_expt_set) == 0);
all_gray_tstarts = find(dtimes < -0.5 & ismember(all_exptvec,gray_expts));
all_gray_tstarts(ismember(all_exptvec(all_gray_tstarts),bar135_expts)) = [];

all_im_expts = find(expt_image_back(cur_expt_set) == 1);
im_expts = find(expt_image_back(cur_expt_set) == 1 & expt_sim_sacs(cur_expt_set) == 0);
all_im_tstarts = find(dtimes < -0.5 & ismember(all_exptvec,im_expts));
all_im_tstarts(ismember(all_exptvec(all_im_tstarts),bar135_expts)) = [];

sim_expts = find(expt_image_back(cur_expt_set) == 1 & expt_sim_sacs(cur_expt_set) > 0);
all_sim_tstarts = find(dtimes < -0.5 & ismember(all_exptvec,sim_expts));
all_sim_tstarts(ismember(all_exptvec(all_sim_tstarts),bar135_expts)) = [];

% %restrict to expts witp unique stim sequences per image back
% all_im_tstarts(~ismember(cur_expt_set(all_exptvec(all_im_tstarts)),52:60)) = [];
% all_sim_tstarts(~ismember(cur_expt_set(all_exptvec(all_sim_tstarts)),52:60)) = [];
% all_sim_sacs(~ismember(cur_expt_set(all_exptvec(all_sim_sacs)),52:60)) = [];
% all_off_stim(~ismember(cur_expt_set(all_exptvec(all_off_stim)),52:60)) = [];
% all_on_stim(~ismember(cur_expt_set(all_exptvec(all_on_stim)),52:60)) = [];

%%
trial_sac_inds = zeros(size(all_t_axis));
trial_sac_inds(interp_saccade_inds) = 1;

trial_msac_inds = zeros(size(all_t_axis));
trial_msac_inds(all_msacs) = 1;

trial_fsac_inds = zeros(size(all_t_axis));
trial_fsac_inds(all_out_sacs) = 1;

trial_ssac_inds = zeros(size(all_t_axis));
trial_ssac_inds(all_ret_sacs) = 1;

big_sacs = unique([all_out_sacs(:); all_ret_sacs(:)]);
trial_bsac_inds = zeros(size(all_t_axis));
trial_bsac_inds(big_sacs) = 1;

trial_simsac_inds = zeros(size(all_t_axis));
trial_simsac_inds(all_sim_sacs) = 1;

trial_onsac_inds = zeros(size(all_t_axis));
trial_onsac_inds(all_on_stim) = 1;

trial_offsac_inds = zeros(size(all_t_axis));
trial_offsac_inds(all_off_stim) = 1;

trial_graystart_inds = zeros(size(all_t_axis));
trial_graystart_inds(all_gray_tstarts) = 1;
trial_imstart_inds = zeros(size(all_t_axis));
trial_imstart_inds(all_im_tstarts) = 1;
trial_simstart_inds = zeros(size(all_t_axis));
trial_simstart_inds(all_sim_tstarts) = 1;

cur_dt = 0.005;
tent_centers = [0:cur_dt:0.9];
tent_centers = round(tent_centers/cur_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

shift = round(0.4/cur_dt);
tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
tent_centers = tent_centers-shift;

trial_sac_inds = zeros(size(all_t_axis));
trial_sac_inds(big_sacs) = 1;
trial_sac_mat = zeros(length(all_t_axis),ntents);
for i = 1:ntents
    trial_sac_mat(:,i) = conv(trial_sac_inds,tbmat(i,:),'same');
end

trial_simsac_mat = zeros(length(all_t_axis),ntents);
for i = 1:ntents
    trial_simsac_mat(:,i) = conv(trial_simsac_inds,tbmat(i,:),'same');
end

trial_msac_inds = zeros(size(all_t_axis));
trial_msac_inds(all_msacs) = 1;
trial_msac_mat = zeros(length(all_t_axis),ntents);
for i = 1:ntents
    trial_msac_mat(:,i) = conv(trial_msac_inds,tbmat(i,:),'same');
end

trial_mssac_inds = zeros(size(all_t_axis));
trial_mssac_inds(short_isi_msacs) = 1;
trial_mssac_mat = zeros(length(all_t_axis),ntents);
for i = 1:ntents
    trial_mssac_mat(:,i) = conv(trial_mssac_inds,tbmat(i,:),'same');
end

trial_mlsac_inds = zeros(size(all_t_axis));
trial_mlsac_inds(long_isi_msacs) = 1;
trial_mlsac_mat = zeros(length(all_t_axis),ntents);
for i = 1:ntents
    trial_mlsac_mat(:,i) = conv(trial_mlsac_inds,tbmat(i,:),'same');
end

%%
all_binned_spks_norm = bsxfun(@rdivide,all_binned_spks,mean(all_binned_spks));
sm_sig = (Fs*0.005);
for cc = 1:96
    all_binned_spks_norm(:,cc) = jmm_smooth_1d_cor(all_binned_spks_norm(:,cc),sm_sig);
end

tr_msac_inds = find(trial_msac_inds(all_tr_inds)==1);
tr_bsac_inds = find(trial_bsac_inds(all_tr_inds)==1);
tr_fsac_inds = find(trial_fsac_inds(all_tr_inds)==1);
tr_ssac_inds = find(trial_ssac_inds(all_tr_inds)==1);
tr_mssac_inds = find(trial_mssac_inds(all_tr_inds)==1);
tr_mlsac_inds = find(trial_mlsac_inds(all_tr_inds)==1);

tr_msac_gray_inds = tr_msac_inds(ismember(all_exptvec(all_tr_inds(tr_msac_inds)),gray_expts));
tr_msac_im_inds = tr_msac_inds(ismember(all_exptvec(all_tr_inds(tr_msac_inds)),all_im_expts));
tr_bsac_gray_inds = tr_bsac_inds(ismember(all_exptvec(all_tr_inds(tr_bsac_inds)),gray_expts));
tr_bsac_im_inds = tr_bsac_inds(ismember(all_exptvec(all_tr_inds(tr_bsac_inds)),all_im_expts));
tr_fsac_gray_inds = tr_fsac_inds(ismember(all_exptvec(all_tr_inds(tr_fsac_inds)),gray_expts));
tr_fsac_im_inds = tr_fsac_inds(ismember(all_exptvec(all_tr_inds(tr_fsac_inds)),all_im_expts));
tr_ssac_gray_inds = tr_ssac_inds(ismember(all_exptvec(all_tr_inds(tr_ssac_inds)),gray_expts));
tr_ssac_im_inds = tr_ssac_inds(ismember(all_exptvec(all_tr_inds(tr_ssac_inds)),all_im_expts));

tr_simsac_inds = find(trial_simsac_inds(all_tr_inds)==1);
tr_onsac_inds = find(trial_onsac_inds(all_tr_inds)==1);
tr_offsac_inds = find(trial_offsac_inds(all_tr_inds)==1);

tr_graystart_inds = find(trial_graystart_inds==1);
tr_imstart_inds = find(trial_imstart_inds==1);
tr_simstart_inds = find(trial_simstart_inds==1);

% for cc = 1:96
%     
%     [avg_sac_trig(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(all_tr_inds,cc),tr_bsac_inds,round(0.4*Fs),round(0.5*Fs));    
%     [avg_msac_trig(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(all_tr_inds,cc),tr_msac_inds,round(0.4*Fs),round(0.5*Fs));
%     [avg_mssac_trig(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(all_tr_inds,cc),tr_mssac_inds,round(0.4*Fs),round(0.5*Fs));
%     [avg_mlsac_trig(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(all_tr_inds,cc),tr_mlsac_inds,round(0.4*Fs),round(0.5*Fs));
%     
% end

disp('Computing big sac stats')
n_boot = 100;
[sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_bsac_inds,round(0.4*Fs),round(0.5*Fs));
sac_trg_mat = reshape(sac_trg_mat,length(tr_bsac_inds),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_bsac_trig_rate = squeeze(mean(boot_avgs));
std_bsac_trig_rate = squeeze(std(boot_avgs));

% disp('Computing big sac gray stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_bsac_gray_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_bsac_gray_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_bsac_gray_trig_rate = squeeze(mean(boot_avgs));
% std_bsac_gray_trig_rate = squeeze(std(boot_avgs));
% 
% disp('Computing big sac im stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_bsac_im_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_bsac_im_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_bsac_im_trig_rate = squeeze(mean(boot_avgs));
% std_bsac_im_trig_rate = squeeze(std(boot_avgs));

disp('Computing micro sac stats')
[sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_msac_inds,round(0.4*Fs),round(0.5*Fs));
sac_trg_mat = reshape(sac_trg_mat,length(tr_msac_inds),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_msac_trig_rate = squeeze(mean(boot_avgs));
std_msac_trig_rate = squeeze(std(boot_avgs));
% 
% disp('Computing micro sac gray stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_msac_gray_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_msac_gray_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_msac_gray_trig_rate = squeeze(mean(boot_avgs));
% std_msac_gray_trig_rate = squeeze(std(boot_avgs));
% 
% disp('Computing micro sac im stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_msac_im_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_msac_im_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_msac_im_trig_rate = squeeze(mean(boot_avgs));
% std_msac_im_trig_rate = squeeze(std(boot_avgs));
% 
% disp('Computing l-ISI micro sac stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_mlsac_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_mlsac_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_mlsac_trig_rate = squeeze(mean(boot_avgs));
% std_mlsac_trig_rate = squeeze(std(boot_avgs));
% 
% disp('Computing s-ISI micro sac stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_mssac_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_mssac_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_mssac_trig_rate = squeeze(mean(boot_avgs));
% std_mssac_trig_rate = squeeze(std(boot_avgs));
% % 
% disp('Computing first sac stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_fsac_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_fsac_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_fsac_trig_rate = squeeze(mean(boot_avgs));
% std_fsac_trig_rate = squeeze(std(boot_avgs));
% disp('Computing first sac gray stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_fsac_gray_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_fsac_gray_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_fsac_gray_trig_rate = squeeze(mean(boot_avgs));
% std_fsac_gray_trig_rate = squeeze(std(boot_avgs));
% disp('Computing first sac im stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_fsac_im_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_fsac_im_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_fsac_im_trig_rate = squeeze(mean(boot_avgs));
% std_fsac_im_trig_rate = squeeze(std(boot_avgs));

% disp('Computing second sac stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_ssac_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_ssac_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_ssac_trig_rate = squeeze(mean(boot_avgs));
% std_ssac_trig_rate = squeeze(std(boot_avgs));
% disp('Computing second sac gray stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_ssac_gray_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_ssac_gray_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_ssac_gray_trig_rate = squeeze(mean(boot_avgs));
% std_ssac_gray_trig_rate = squeeze(std(boot_avgs));
% disp('Computing second sac im stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_ssac_im_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_ssac_im_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_ssac_im_trig_rate = squeeze(mean(boot_avgs));
% std_ssac_im_trig_rate = squeeze(std(boot_avgs));

disp('Computing sim sac stats')
[sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_simsac_inds,round(0.4*Fs),round(0.5*Fs));
sac_trg_mat = reshape(sac_trg_mat,length(tr_simsac_inds),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_simsac_trig_rate = squeeze(mean(boot_avgs));
std_simsac_trig_rate = squeeze(std(boot_avgs));

% disp('Computing on sac stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_onsac_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_onsac_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_onsac_trig_rate = squeeze(mean(boot_avgs));
% std_onsac_trig_rate = squeeze(std(boot_avgs));
% 
% disp('Computing off sac stats')
% [sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_offsac_inds,round(0.4*Fs),round(0.5*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_offsac_inds),96*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
% avg_offsac_trig_rate = squeeze(mean(boot_avgs));
% std_offsac_trig_rate = squeeze(std(boot_avgs));
% 
% disp('Computing tstart')
% [sac_trg_mat,cur_lagsf] = get_event_trig_mat(all_binned_spks_norm,tr_graystart_inds,round(0*Fs),round(2*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_graystart_inds),96*length(cur_lagsf));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lagsf) 96]);
% avg_graystart_trig_rate = squeeze(mean(boot_avgs));
% std_graystart_trig_rate = squeeze(std(boot_avgs));
% 
% [sac_trg_mat,cur_lagsf] = get_event_trig_mat(all_binned_spks_norm,tr_imstart_inds,round(0*Fs),round(2*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_imstart_inds),96*length(cur_lagsf));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lagsf) 96]);
% avg_imstart_trig_rate = squeeze(mean(boot_avgs));
% std_imstart_trig_rate = squeeze(std(boot_avgs));
% 
% [sac_trg_mat,cur_lagsf] = get_event_trig_mat(all_binned_spks_norm,tr_simstart_inds,round(0*Fs),round(2*Fs));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_simstart_inds),96*length(cur_lagsf));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lagsf) 96]);
% avg_simstart_trig_rate = squeeze(mean(boot_avgs));
% std_simstart_trig_rate = squeeze(std(boot_avgs));

%%
figure; hold on
plot(cur_lags*new_dt,mean(avg_bsac_trig_rate,2),'k','linewidth',1)
plot(cur_lags*new_dt,mean(avg_msac_trig_rate,2),'r','linewidth',1)
plot(cur_lags*new_dt,mean(avg_simsac_trig_rate,2),'b','linewidth',1)
legend('Guided','Micro','Simulated')
shadedErrorBar(cur_lags*new_dt,mean(avg_bsac_trig_rate,2),std(avg_bsac_trig_rate,[],2)/sqrt(96),{'k'})
shadedErrorBar(cur_lags*new_dt,mean(avg_msac_trig_rate,2),std(avg_msac_trig_rate,[],2)/sqrt(96),{'r'})
shadedErrorBar(cur_lags*new_dt,mean(avg_simsac_trig_rate,2),std(avg_simsac_trig_rate,[],2)/sqrt(96),{'b'})
xlabel('Time since saccade onset (s)','fontsize',16)
ylabel('Relative rate','fontsize',16)
xlim([-0.3 0.4])
xl = xlim();
line(xl,[1 1],'color','k')
yl = ylim();
line([0 0],yl,'color','k')

figure; hold on
plot(cur_lags*new_dt,mean(avg_bsac_gray_trig_rate,2),'k','linewidth',1)
plot(cur_lags*new_dt,mean(avg_bsac_im_trig_rate,2),'r','linewidth',1)
plot(cur_lags*new_dt,mean(avg_simsac_trig_rate,2),'b','linewidth',1)
legend('Gray back','Image back','Simulated')
shadedErrorBar(cur_lags*new_dt,mean(avg_bsac_gray_trig_rate,2),std(avg_bsac_gray_trig_rate,[],2)/sqrt(96),{'k'})
shadedErrorBar(cur_lags*new_dt,mean(avg_bsac_im_trig_rate,2),std(avg_bsac_im_trig_rate,[],2)/sqrt(96),{'r'})
shadedErrorBar(cur_lags*new_dt,mean(avg_simsac_trig_rate,2),std(avg_simsac_trig_rate,[],2)/sqrt(96),{'b'})
xlabel('Time since saccade onset (s)','fontsize',16)
ylabel('Relative rate','fontsize',16)
xlim([-0.3 0.4])
xl = xlim();
line(xl,[1 1],'color','k')
yl = ylim();
line([0 0],yl,'color','k')

figure; hold on
plot(cur_lags*new_dt,mean(avg_msac_gray_trig_rate,2),'k','linewidth',1)
plot(cur_lags*new_dt,mean(avg_msac_im_trig_rate,2),'r','linewidth',1)
legend('Micro gray back','Micro image back')
shadedErrorBar(cur_lags*new_dt,mean(avg_msac_gray_trig_rate,2),std(avg_msac_gray_trig_rate,[],2)/sqrt(96),{'k'})
shadedErrorBar(cur_lags*new_dt,mean(avg_msac_im_trig_rate,2),std(avg_msac_im_trig_rate,[],2)/sqrt(96),{'r'})
xlabel('Time since saccade onset (s)','fontsize',16)
ylabel('Relative rate','fontsize',16)
xlim([-0.3 0.4])
xl = xlim();
line(xl,[1 1],'color','k')
yl = ylim();
line([0 0],yl,'color','k')

figure; hold on
plot(cur_lags*new_dt,mean(avg_fsac_trig_rate,2),'k','linewidth',1)
plot(cur_lags*new_dt,mean(avg_ssac_trig_rate,2),'r','linewidth',1)
legend('Outgoing','Returning')
shadedErrorBar(cur_lags*new_dt,mean(avg_fsac_trig_rate,2),std(avg_fsac_trig_rate,[],2)/sqrt(96),{'k'})
shadedErrorBar(cur_lags*new_dt,mean(avg_ssac_trig_rate,2),std(avg_ssac_trig_rate,[],2)/sqrt(96),{'r'})
xlabel('Time since saccade onset (s)','fontsize',16)
ylabel('Relative rate','fontsize',16)
xlim([-0.3 0.4])
xl = xlim();
line(xl,[1 1],'color','k')
yl = ylim();
line([0 0],yl,'color','k')

figure; hold on
plot(cur_lags*new_dt,mean(avg_onsac_trig_rate,2),'k','linewidth',1)
plot(cur_lags*new_dt,mean(avg_offsac_trig_rate,2),'r','linewidth',1)
legend('Outgoing sim','Returning sim')
shadedErrorBar(cur_lags*new_dt,mean(avg_onsac_trig_rate,2),std(avg_onsac_trig_rate,[],2)/sqrt(96),{'k'})
shadedErrorBar(cur_lags*new_dt,mean(avg_offsac_trig_rate,2),std(avg_offsac_trig_rate,[],2)/sqrt(96),{'r'})
xlabel('Time since saccade onset (s)','fontsize',16)
ylabel('Relative rate','fontsize',16)
xlim([-0.3 0.4])
xl = xlim();
line(xl,[1 1],'color','k')
yl = ylim();
line([0 0],yl,'color','k')

%%
% poss_ses = [...
%     403004
%     419002
%     433008
%     434005
%     435002
%     503002
%     505006
%     511007
%     562001
%     572001];
% 
% for ee = 1:length(poss_ses);
%     ee
%     all_im_expts = find(expt_image_back(cur_expt_set) == 1);
%     all_im_tstarts = find(dtimes < -0.5 & ...
%         ismember(all_exptvec,all_im_expts) & all_se == poss_ses(ee));
%     
%     n_stims(ee) = length(all_im_tstarts);
%     [sac_trg_mat,cur_lagsf] = get_event_trig_mat(all_binned_spks_norm,all_im_tstarts,round(0*Fs),round(2*Fs));
%     sac_trg_mat = reshape(sac_trg_mat,length(all_im_tstarts),96*length(cur_lagsf));
%     boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
%     boot_avgs = reshape(boot_avgs,[n_boot length(cur_lagsf) 96]);
%     avg_stim_dep_trig_rate(ee,:,:) = squeeze(mean(boot_avgs));
%     std_stim_dep_trig_rate(ee,:,:) = squeeze(std(boot_avgs));
%     
% end
% 
% for cc =1 :96
%     cc
%     for ee = 1:length(poss_ses)
%         subplot(5,2,ee)
%         shadedErrorBar(cur_lagsf*new_dt,squeeze(avg_stim_dep_trig_rate(ee,:,cc)),squeeze(std_stim_dep_trig_rate(ee,:,cc)))
%     end
%     pause
%     clf
% end
%%
cd ~/Analysis/bruce/G081/sac_trg_avgs/

close all
% for cc = 1:96
%     cc
%     subplot(2,1,1)
%     h=errorbar(cur_lags/Fs,squeeze(mean(sac_trig_avg_rate(cc,:,:),2)),squeeze(std(sac_trig_avg_rate(cc,:,:),[],2))/sqrt(10));
%     errorbar_tick(h,.001,'units')
%     xlim([-0.3 0.4])
%     subplot(2,1,2)
%     h=errorbar(cur_lags/Fs,squeeze(mean(msac_trig_avg_rate(cc,:,:),2)),squeeze(std(msac_trig_avg_rate(cc,:,:),[],2))/sqrt(10));
%     errorbar_tick(h,.001,'units')
%     xlim([-0.3 0.4])
%     
%     pause
%     clf
% end

for cc = 1:96
    cc
    f = figure('visible','off');
%     subplot(1,3,1)
    hold on
    shadedErrorBar(cur_lags*new_dt,avg_msac_trig_rate(:,cc),std_msac_trig_rate(:,cc),{'r'},1);

    shadedErrorBar(cur_lags*new_dt,avg_bsac_trig_rate(:,cc),std_bsac_trig_rate(:,cc),{'b'},1);
    shadedErrorBar(cur_lags*new_dt,avg_simsac_trig_rate(:,cc),std_simsac_trig_rate(:,cc),{'k'},1);
    xlim([-0.3 0.4]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k')
    line([0 0],yl,'color','k')    
    
%      subplot(1,3,2)
%     hold on
%     shadedErrorBar(cur_lags*new_dt,avg_msac_gray_trig_rate(:,cc),std_msac_gray_trig_rate(:,cc),{'r'},1);
% 
%     shadedErrorBar(cur_lags*new_dt,avg_bsac_gray_trig_rate(:,cc),std_bsac_gray_trig_rate(:,cc),{'b'},1);
% %     shadedErrorBar(cur_lags*new_dt,avg_simsac_trig_rate(:,cc),std_simsac_trig_rate(:,cc),{'k'},1);
%     xlim([-0.3 0.4]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k')
%     line([0 0],yl,'color','k')    
%    
%      subplot(1,3,3)
%     hold on
%     shadedErrorBar(cur_lags*new_dt,avg_msac_im_trig_rate(:,cc),std_msac_im_trig_rate(:,cc),{'r'},1);
% 
%     shadedErrorBar(cur_lags*new_dt,avg_bsac_im_trig_rate(:,cc),std_bsac_im_trig_rate(:,cc),{'b'},1);
%     shadedErrorBar(cur_lags*new_dt,avg_simsac_trig_rate(:,cc),std_simsac_trig_rate(:,cc),{'k'},1);
%     xlim([-0.3 0.4]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k')
%     line([0 0],yl,'color','k')    

fname = sprintf('Unit%d',cc);
print(fname,'-dpng');
close


%     pause
%     clf
end

% for cc = 1:96
%     cc
% %     subplot(1,3,1)
%     hold on
%     shadedErrorBar(cur_lags*new_dt,avg_msac_trig_rate(:,cc),std_msac_trig_rate(:,cc),{'r'},1);
% 
%     shadedErrorBar(cur_lags*new_dt,avg_bsac_trig_rate(:,cc),std_bsac_trig_rate(:,cc),{'b'},1);
%     shadedErrorBar(cur_lags*new_dt,avg_simsac_trig_rate(:,cc),std_simsac_trig_rate(:,cc),{'k'},1);
%     xlim([-0.3 0.4]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k')
%     line([0 0],yl,'color','k')    
%     
% %      subplot(1,3,2)
% %     hold on
% % %     shadedErrorBar(cur_lags*new_dt,avg_msac_trig_rate(:,cc),std_msac_trig_rate(:,cc),{'r'},1);
% % 
% % %     shadedErrorBar(cur_lags*new_dt,avg_fsac_trig_rate(:,cc),std_fsac_trig_rate(:,cc),{'b'},1);
% %     shadedErrorBar(cur_lags*new_dt,avg_fsac_gray_trig_rate(:,cc),std_fsac_gray_trig_rate(:,cc),{'b'},1);
% %     shadedErrorBar(cur_lags*new_dt,avg_fsac_im_trig_rate(:,cc),std_fsac_im_trig_rate(:,cc),{'r'},1);
% %     shadedErrorBar(cur_lags*new_dt,avg_onsac_trig_rate(:,cc),std_onsac_trig_rate(:,cc),{'k'},1);
% %     xlim([-0.3 0.4]);
% %     xl = xlim(); yl = ylim();
% %     line(xl,[1 1],'color','k')
% %     line([0 0],yl,'color','k')    
% %    
% %      subplot(1,3,3)
% %     hold on
% % %     shadedErrorBar(cur_lags*new_dt,avg_msac_trig_rate(:,cc),std_msac_trig_rate(:,cc),{'r'},1);
% % 
% % %     shadedErrorBar(cur_lags*new_dt,avg_ssac_trig_rate(:,cc),std_ssac_trig_rate(:,cc),{'b'},1);
% %     shadedErrorBar(cur_lags*new_dt,avg_ssac_gray_trig_rate(:,cc),std_ssac_gray_trig_rate(:,cc),{'b'},1);
% %     shadedErrorBar(cur_lags*new_dt,avg_ssac_im_trig_rate(:,cc),std_ssac_im_trig_rate(:,cc),{'r'},1);
% %     shadedErrorBar(cur_lags*new_dt,avg_offsac_trig_rate(:,cc),std_offsac_trig_rate(:,cc),{'k'},1);
% %     xlim([-0.3 0.4]);
% %     xl = xlim(); yl = ylim();
% %     line(xl,[1 1],'color','k')
% %     line([0 0],yl,'color','k')    
% 
%     pause
%     clf
% end

%%
close all
for cc = 1:96
    cc
    hold on
    subplot(3,1,1)
    shadedErrorBar(cur_lagsf*new_dt,avg_graystart_trig_rate(:,cc),std_graystart_trig_rate(:,cc),{'r'},1);
    xlim([0 2]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k')
    line([0.7 0.7],yl,'color','k')    
    line([1.4 1.4],yl,'color','k')    

    
    subplot(3,1,2)
    shadedErrorBar(cur_lagsf*new_dt,avg_imstart_trig_rate(:,cc),std_imstart_trig_rate(:,cc),{'b'},1);
    xlim([0 2]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k')
    line([0.7 0.7],yl,'color','k')    
    line([1.4 1.4],yl,'color','k')    

    subplot(3,1,3)
    shadedErrorBar(cur_lagsf*new_dt,avg_simstart_trig_rate(:,cc),std_simstart_trig_rate(:,cc),{'k'},1);
    xlim([0 2]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k')
    line([0.7 0.7],yl,'color','k')    
    line([1.4 1.4],yl,'color','k')    
    
    pause
    clf
end

%%
NL_type = 0;
dl_time = 10;
dl2_time = 1000;
silent = 1;

cur_msac_inds = find(trial_msac_inds(all_tr_inds)==1);
cur_bsac_inds = find(trial_bsac_inds(all_tr_inds)==1);
cur_simsac_inds = find(trial_simsac_inds(all_tr_inds)==1);
cur_mssac_inds = find(trial_mssac_inds(all_tr_inds)==1);
cur_mlsac_inds = find(trial_mlsac_inds(all_tr_inds)==1);

for cc = 1:96
    fprintf('Cell %d of %d\n',cc,96);
    
    Robs = all_binned_spks(all_tr_inds,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    cur_Xmat = [lin_X(all_tr_inds,:)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_null,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);

%     cur_Xmat = [trial_sac_mat(all_tr_inds,:) lin_X(all_tr_inds,:)];
%     lamrange = [dl_time 1 ntents 0];
%     lamrange2 = [dl2_time 1 ntents 0];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_sac,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, lamrange, lamrange2,[], [], [], [], NL_type);
%     sac_kern(cc,:) = fitp_sac.k(1:ntents);
% %     sac_kern_se(cc,:) = se(1:ntents);
% %     predrate_out = cur_Xmat*fitp_sac.k(1:end-1) + fitp_sac.k(end);
% %     [st_avg_sacmod(cc,:),cur_lags] = get_event_trig_avg(exp(predrate_out),tr_sac_inds,round(0.2/new_dt),round(0.5/new_dt));
    
%     cur_Xmat = [trial_msac_mat(all_tr_inds,:) trial_sac_mat(all_tr_inds,:) lin_X(all_tr_inds,:)];
%     lamrange = [dl_time 1 ntents 0; dl_time ntents+1 2*ntents 0];
%     lamrange2 = [dl2_time 1 ntents 0; dl2_time ntents+1 2*ntents 0];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_mulsac,grad,se] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, lamrange, lamrange2,[], [], [], [], NL_type);
%     co_msac_kern(cc,:) = fitp_mulsac.k(1:ntents);
%     co_bsac_kern(cc,:) = fitp_mulsac.k(ntents+1:2*ntents);
%     co_msac_kern_se(cc,:) = se(1:ntents);
%     co_bsac_kern_se(cc,:) = se(ntents+1:2*ntents);
%     predrate_out = cur_Xmat*fitp_mulsac.k(1:end-1) + fitp_mulsac.k(end);
%     predrate_out = log(1+exp(predrate_out));
%     [co_avg_bsacmod(cc,:),cur_lags] = get_event_trig_avg(predrate_out,cur_bsac_inds,round(0.4/cur_dt),round(0.5/cur_dt));
%     [co_avg_msacmod(cc,:),cur_lags] = get_event_trig_avg(predrate_out,cur_msac_inds,round(0.4/cur_dt),round(0.5/cur_dt));
%     [co_avg_mssacmod(cc,:),cur_lags] = get_event_trig_avg(predrate_out,cur_mssac_inds,round(0.4/cur_dt),round(0.5/cur_dt));
%     [co_avg_mlsacmod(cc,:),cur_lags] = get_event_trig_avg(predrate_out,cur_mlsac_inds,round(0.4/cur_dt),round(0.5/cur_dt));
    
    cur_Xmat = [trial_msac_mat(all_tr_inds,:) trial_sac_mat(all_tr_inds,:) trial_simsac_mat(all_tr_inds,:) lin_X(all_tr_inds,:)];
    lamrange = [dl_time 1 ntents 0; dl_time ntents+1 2*ntents 0; dl_time 2*ntents+1 3*ntents 0];
    lamrange2 = [dl2_time 1 ntents 0; dl2_time ntents+1 2*ntents 0; dl2_time 2*ntents+1 3*ntents 0];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_mulsac,grad,se] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, lamrange, lamrange2,[], [], [], [], NL_type);
    co_msac_kern(cc,:) = fitp_mulsac.k(1:ntents);
    co_bsac_kern(cc,:) = fitp_mulsac.k(ntents+1:2*ntents);
    co_simsac_kern(cc,:) = fitp_mulsac.k(2*ntents+1:3*ntents);
    co_msac_kern_se(cc,:) = se(1:ntents);
    co_bsac_kern_se(cc,:) = se(ntents+1:2*ntents);
    co_simsac_kern_se(cc,:) = se(2*ntents+1:3*ntents);
    predrate_out = cur_Xmat*fitp_mulsac.k(1:end-1) + fitp_mulsac.k(end);
    predrate_out = log(1+exp(predrate_out));
    [co_avg_bsacmod(cc,:),cur_lags] = get_event_trig_avg(predrate_out,cur_bsac_inds,round(0.4/cur_dt),round(0.5/cur_dt));
    [co_avg_simsacmod(cc,:),cur_lags] = get_event_trig_avg(predrate_out,cur_simsac_inds,round(0.4/cur_dt),round(0.5/cur_dt));
    [co_avg_msacmod(cc,:),cur_lags] = get_event_trig_avg(predrate_out,cur_msac_inds,round(0.4/cur_dt),round(0.5/cur_dt));
    [co_avg_mssacmod(cc,:),cur_lags] = get_event_trig_avg(predrate_out,cur_mssac_inds,round(0.4/cur_dt),round(0.5/cur_dt));
    [co_avg_mlsacmod(cc,:),cur_lags] = get_event_trig_avg(predrate_out,cur_mlsac_inds,round(0.4/cur_dt),round(0.5/cur_dt));

%     
%      cur_Xmat = [trial_msac_mat(all_tr_inds,:) trial_sac_mat(all_tr_inds,:) lin_X(all_tr_inds,:)];
%     lamrange = [dl_time 1 ntents 0; dl_time ntents+1 2*ntents 0];
%     lamrange2 = [dl2_time 1 ntents 0; dl2_time ntents+1 2*ntents 0];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_mulsac,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, lamrange, lamrange2,[], [], [], [], NL_type);
%     co_msac_kern(cc,:) = fitp_mulsac.k(1:ntents);
%     co_bsac_kern(cc,:) = fitp_mulsac.k(ntents+1:2*ntents);
%     predrate_out = cur_Xmat*fitp_mulsac.k(1:end-1) + fitp_mulsac.k(end);
%     [co_avg_msacmod(cc,:),cur_lags] = get_event_trig_avg(exp(predrate_out),cur_msac_inds,round(0.4/cur_dt),round(0.5/cur_dt));
    
end

%%
cd /home/james/Desktop

close all
for cc = 1:96
    fprintf('Cell %d of %d\n',cc,96);
    
    cur_avg_rate = mean(all_binned_spks(:,cc));
    subplot(3,1,1)
    hold on
    shadedErrorBar(cur_lags*new_dt,avg_msac_trig_rate(:,cc),std_msac_trig_rate(:,cc),{'r'},1);
%     shadedErrorBar(cur_lags*new_dt,avg_mlsac_trig_rate(:,cc),std_mlsac_trig_rate(:,cc),{'k'},1);
%     shadedErrorBar(cur_lags*new_dt,avg_msac_trig_rate(:,cc),std_msac_trig_rate(:,cc),{'b'},1);
    
axis tight
    xlim([-0.3 0.35]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k')
    line([0 0],yl,'color','k')
%     plot(cur_lags*new_dt,co_avg_msacmod(cc,:)/cur_avg_rate,'k')
    xlabel('Time since saccade onset (s)','fontsize',16)
    ylabel('Relative rate','fontsize',16)
    title('Microsaccades','fontsize',16)

    subplot(3,1,2)
    hold on
    shadedErrorBar(cur_lags*new_dt,avg_bsac_trig_rate(:,cc),std_bsac_trig_rate(:,cc),{'r'},1);
%     shadedErrorBar(cur_lags*new_dt,avg_fsac_trig_rate(:,cc),std_fsac_trig_rate(:,cc),{'b'},1);
%     shadedErrorBar(cur_lags*new_dt,avg_ssac_trig_rate(:,cc),std_ssac_trig_rate(:,cc),{'k'},1);
    plot(cur_lags*new_dt,avg_fsac_trig_rate(:,cc),'b','linewidth',1);
    plot(cur_lags*new_dt,avg_ssac_trig_rate(:,cc),'k','linewidth',1);
axis tight
    xlim([-0.3 0.35]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k')
    line([0 0],yl,'color','k')
    hold on
%     plot(cur_lags*new_dt,co_avg_bsacmod(cc,:)/cur_avg_rate,'k')
    xlabel('Time since saccade onset (s)','fontsize',16)
    ylabel('Relative rate','fontsize',16)
    title('Guided saccades','fontsize',16)

    subplot(3,1,3)
    hold on
    shadedErrorBar(cur_lags*new_dt,avg_simsac_trig_rate(:,cc),std_simsac_trig_rate(:,cc),{'r'},1);
%     shadedErrorBar(cur_lags*new_dt,avg_onsac_trig_rate(:,cc),std_onsac_trig_rate(:,cc),{'b'},1);
%     shadedErrorBar(cur_lags*new_dt,avg_offsac_trig_rate(:,cc),std_offsac_trig_rate(:,cc),{'k'},1);
    plot(cur_lags*new_dt,avg_onsac_trig_rate(:,cc),'b','linewidth',1);
    plot(cur_lags*new_dt,avg_offsac_trig_rate(:,cc),'k','linewidth',1);
axis tight
    xlim([-0.3 0.35]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k')
    line([0 0],yl,'color','k')
    hold on
%     plot(cur_lags*new_dt,co_avg_simsacmod(cc,:)/cur_avg_rate,'k')
    xlabel('Time since saccade onset (s)','fontsize',16)
    ylabel('Relative rate','fontsize',16)
    title('Simulated saccades','fontsize',16)

%     subplot(2,3,4)
% %     shadedErrorBar(tent_centers*new_dt,co_msac_kern(cc,:),co_msac_kern_se(cc,:));
%     plot(tent_centers*new_dt,co_msac_kern(cc,:));
%     xlim([-0.4 0.4]);
%        xl = xlim(); yl = ylim();
%     line(xl,[0 0],'color','k')
%     line([0 0],yl,'color','k')
%  
%     subplot(2,3,5)
% %     shadedErrorBar(tent_centers*new_dt,co_bsac_kern(cc,:),co_bsac_kern_se(cc,:));
%     plot(tent_centers*new_dt,co_bsac_kern(cc,:));
%     xlim([-0.4 0.4]);
%        xl = xlim(); yl = ylim();
%     line(xl,[0 0],'color','k')
%     line([0 0],yl,'color','k')
% 
%     subplot(2,3,6)
% %     shadedErrorBar(tent_centers*new_dt,co_bsac_kern(cc,:),co_bsac_kern_se(cc,:));
%     plot(tent_centers*new_dt,co_simsac_kern(cc,:));
%     xlim([-0.4 0.4]);
%        xl = xlim(); yl = ylim();
%     line(xl,[0 0],'color','k')
%     line([0 0],yl,'color','k')
% 


    fillPage(gcf,'Papersize',[5 10]);
    pname = sprintf('Array_sacmod_unit%d',cc);
    print(pname,'-dpdf');
    close
    
    
%     pause
%     clf
end

%%
cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans
gabor_params_used = gabor_params_f{2};
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

cd ~/Data/bruce/G081

%fit smoothed retinotopic surface
interp_x = nan(96,1);
interp_y = nan(96,1);
tempinds = zeros(10,10);
for i = 1:96
    tempx(Y_pos(i),X_pos(i)) = gabor_params_used(i,1);
    tempy(Y_pos(i),X_pos(i)) = gabor_params_used(i,2);
    tempinds(Y_pos(i),X_histpos(i)) = i;
end
weights = ones(10,10);
weights(tempx==0) = 0;
xpos_interp = smoothn(tempx,weights,'robust');
ypos_interp = smoothn(tempy,weights,'robust');
used_inds = find(weights == 1);
tempinds = tempinds(used_inds);
interp_x(tempinds) = xpos_interp(used_inds);
interp_y(tempinds) = ypos_interp(used_inds);
xi = linspace(min(interp_x),max(interp_x),50);
yi = linspace(min(interp_y),max(interp_y),50);
[Xi,Yi] = meshgrid(xi,yi);

id_mat = nan(10,10);
for i = 1:10
    for j = 1:10
        cur = find(X_pos==j&Y_pos==i);
        if ~isempty(cur)
            id_mat(i,j) = cur;
        end
    end
end
use_ids = find(~isnan(id_mat));

el_fdist = sqrt(interp_x.^2 + interp_y.^2);
[temp,f_ord] = sort(el_fdist);

cent_point = [mean(interp_x) mean(interp_y)];
el_edgedist = sqrt((interp_x-cent_point(1)).^2 + (interp_y-mean(cent_point(2))).^2);
[temp,e_ord] = sort(el_edgedist);

%%
bar_oris = [0 45 90 135];
for bb = 1:4
    cur_bar_expts = find(expt_bar_ori(cur_expt_set) == bar_oris(bb));
    cur_inds = all_tr_inds(ismember(all_exptvec(all_tr_inds),cur_bar_expts));
    bar_msac_inds = find(trial_msac_inds(cur_inds)==1);
    bar_bsac_inds = find(trial_bsac_inds(cur_inds)==1);
    bar_fsac_inds = find(trial_fsac_inds(cur_inds)==1);
    bar_ssac_inds = find(trial_ssac_inds(cur_inds)==1);
        for cc = 1:96
        [bar_sac_trig_avg_rate(cc,bb,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(cur_inds,cc),bar_bsac_inds,round(0.4*Fs),round(0.5*Fs));       
        [bar_msac_trig_avg_rate(cc,bb,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(cur_inds,cc),bar_msac_inds,round(0.4*Fs),round(0.5*Fs));
        [bar_fsac_trig_avg_rate(cc,bb,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(cur_inds,cc),bar_fsac_inds,round(0.4*Fs),round(0.5*Fs));
        [bar_ssac_trig_avg_rate(cc,bb,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(cur_inds,cc),bar_ssac_inds,round(0.4*Fs),round(0.5*Fs));
    end
end

%%
cur_expts = find(expt_image_back(cur_expt_set) == 1);
cur_inds = all_tr_inds(ismember(all_exptvec(all_tr_inds),cur_expts));
cur_msac_inds = find(trial_msac_inds(cur_inds)==1);
cur_bsac_inds = find(trial_bsac_inds(cur_inds)==1);
for cc = 1:96
    [im_sac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(cur_inds,cc),cur_bsac_inds,round(0.4*Fs),round(0.5*Fs));
    
    [im_msac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(cur_inds,cc),cur_msac_inds,round(0.4*Fs),round(0.5*Fs));
end

cur_expts = find(expt_image_back(cur_expt_set) == 0);
cur_inds = all_tr_inds(ismember(all_exptvec(all_tr_inds),cur_expts));
cur_msac_inds = find(trial_msac_inds(cur_inds)==1);
cur_bsac_inds = find(trial_bsac_inds(cur_inds)==1);
for cc = 1:96
    [nim_sac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(cur_inds,cc),cur_bsac_inds,round(0.4*Fs),round(0.5*Fs));
    
    [nim_msac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(cur_inds,cc),cur_msac_inds,round(0.4*Fs),round(0.5*Fs));
end

cur_simsac_inds = find(trial_simsac_inds(all_tr_inds)==1);
cur_onsac_inds = find(trial_onsac_inds(all_tr_inds)==1);
cur_offsac_inds = find(trial_offsac_inds(all_tr_inds)==1);
for cc = 1:96
    [sim_sac_trig_avg_rate(cc,:),cur_lags2] = get_event_trig_avg(all_binned_spks_norm(all_tr_inds,cc),cur_simsac_inds,round(0.5*Fs),round(0.75*Fs));
    [on_sac_trig_avg_rate(cc,:),cur_lags2] = get_event_trig_avg(all_binned_spks_norm(all_tr_inds,cc),cur_onsac_inds,round(0.5*Fs),round(0.75*Fs));
    [off_sac_trig_avg_rate(cc,:),cur_lags2] = get_event_trig_avg(all_binned_spks_norm(all_tr_inds,cc),cur_offsac_inds,round(0.5*Fs),round(0.75*Fs));
end
%%
close all
for cc = 1:96
    cc
    plot(cur_lags/Fs,avg_sac_trig(cc,:))
    hold on
    plot(cur_lags/Fs,avg_msac_trig(cc,:),'r')
    plot(cur_lags2/Fs,sim_sac_trig_avg_rate(cc,:),'k')
    
    xlim([-0.3 0.4])

    %     subplot(2,1,1)
%     plot(cur_lags/Fs,squeeze(bar_sac_trig_avg_rate(cc,:,:)))
%     xlim([-0.3 0.4])
%     subplot(2,1,2)
%     plot(cur_lags/Fs,squeeze(bar_msac_trig_avg_rate(cc,:,:)))
%     xlim([-0.3 0.4])
    
    pause
    clf
end



%%
close all
for i = 1:96
subplot(2,1,1)
shadedErrorBar(cur_lags/Fs,squeeze(mean(msac_trig_avg_rate(i,:,:),2)),squeeze(std(msac_trig_avg_rate(i,:,:),[],2))/sqrt(10))
subplot(2,1,2)
plot(tent_centers*cur_dt,co_msac_kern(i,:))
pause
clf
end