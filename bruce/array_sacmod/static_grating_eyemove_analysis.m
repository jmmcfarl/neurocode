clear all
% close all

Expt_name = 'G087';
data_dir = ['~/Data/bruce/' Expt_name];
addpath('~/James_scripts/bruce/G081/')

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = '_rph_gratings';

%%
stim_fs = 100; %in Hz
dt = 0.005;
beg_buffer = round(0.15/dt);
end_buffer = round(0.15/dt);

%%
include_expts = {'square.or', 'square.orXnph'};
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    included_type(ii) = any(strcmp(expt_names{ii},include_expts));
end

use_expts = find(included_type);
use_expts(use_expts==36) = []; %no orientation stored
%%
cd(data_dir);

trial_cnter = 0;

all_stim_times = [];
all_used_inds = false(0);
all_phase_Xmat = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
all_trial_ors = [];
all_trial_start_times = [];
all_trial_end_times = [];
for ee = 1:length(use_expts);
    fprintf('Expt %d of %d\n',ee,length(use_expts));
    cur_expt = use_expts(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
        
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    
    [un_ids,id_inds] = unique(trial_ids);
    trial_durs = trial_durs(id_inds);
    trial_start_times = trial_start_times(id_inds);
    trial_end_times = trial_end_times(id_inds);
    
    trial_ors = [Expts{cur_expt}.Trials(id_inds).or];
    if isfield(Expts{cur_expt}.Trials(1),'nph')
        trial_nph = [Expts{cur_expt}.Trials(id_inds).nph];
    else
       trial_nph = zeros(size(trial_ors)); 
    end
    use_trials = find(trial_durs >= 0.5 & trial_nph == 0);
    
    all_trial_ors = [all_trial_ors; trial_ors(use_trials)'];
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start/1e4;
        cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_binned_spks = nan(length(cur_t_edges),96);
        for cc = 1:96
            cur_hist = histc(Clusters{cc}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist;
            %             temp = convert_to_spikebins(cur_hist(1:end-1));
        end
        
        
        cur_used_inds = true(length(cur_t_edges),1);
        cur_used_inds(1:beg_buffer) = false;
        cur_used_inds((end-end_buffer+1):end) = false;
        
        all_stim_times = [all_stim_times; cur_t_edges];
        all_used_inds = [all_used_inds; cur_used_inds];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_t_edges))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_t_edges))*(tt + trial_cnter)];
    end
    trial_cnter = trial_cnter + n_trials;
end

%%
[c,ia,ic] = unique(all_trialvec);
n_trials = length(ia);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_set))';

tr_inds(~all_used_inds(tr_inds)) = [];
xv_inds(~all_used_inds(xv_inds)) = [];

Xexpt = zeros(length(all_stim_times),length(use_expts)-1);
for i = 1:length(use_expts)-1
    cur_set = find(all_exptvec==i);
    Xexpt(cur_set,i) = 1;
end

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
emfile = ['jbe' Expt_name '.em.mat'];
load(emfile);

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
for ee = 1:length(use_expts);
    fprintf('Loading eye data for expt %d of %d\n',ee,length(use_expts));
    cur_set = find(all_exptvec==ee);
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_stim_times(cur_set([1 end])));

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

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_stim_times);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_stim_times);

all_eye_intrial = zeros(size(all_eye_ts));
for i = 1:length(all_trial_start_times)
    curset = find(all_eye_ts >= all_trial_start_times(i) & all_eye_ts <= all_trial_end_times(i));
    all_eye_intrial(curset) = 1;
end

%%
%compute times saccades
sac_thresh = 10;
% min_isi = 0.05;
min_isi = 0.05;
max_isi = Inf;
% max_isi = 0.2;

orig_saccade_inds = 1 + find(all_eye_speed(1:end-1) < sac_thresh & all_eye_speed(2:end) > sac_thresh);
orig_saccade_inds = unique(orig_saccade_inds);
isis = [Inf; diff(orig_saccade_inds)]/eye_fs;
bad_isis = find(isis < min_isi | isis > max_isi);
orig_saccade_inds(bad_isis) = [];
isis(bad_isis) = [];

burst_isi_inds = find(isis(1:end-1) < 0.15 | isis(2:end) < 0.15);
% orig_saccade_inds = orig_saccade_inds(burst_isi_inds);
% isis = isis(burst_isi_inds);
orig_saccade_inds(burst_isi_inds) = [];
isis(burst_isi_inds) = [];


isis(all_eye_intrial(orig_saccade_inds) == 0) = [];
orig_saccade_inds(all_eye_intrial(orig_saccade_inds) == 0) = [];
interp_saccade_inds = round(interp1(all_stim_times,1:length(all_stim_times),all_eye_ts(orig_saccade_inds)));

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

sac_starts = sac_starts - round(0.01/dt);
sac_stops = sac_stops + round(0.01/dt);

pre_pos = interp_eye_vals(sac_starts,1:2);
post_pos = interp_eye_vals(sac_stops,1:2);
sac_delta_pos = post_pos - pre_pos;
sac_amp = sqrt(sum(sac_delta_pos.^2,2));
sac_dir = atan(sac_delta_pos(:,2)./sac_delta_pos(:,1));

micro_inds = find(sac_amp <= 1);
sac_stim_ori = all_trial_ors(all_trialvec(interp_saccade_inds));

hor_stim_micros = micro_inds(sac_stim_ori(micro_inds) == 0);
ver_stim_micros = micro_inds(sac_stim_ori(micro_inds) == 90);

hor_sac_micros = micro_inds(sac_dir(micro_inds) > -pi/8 & sac_dir(micro_inds) < pi/8);
ver_sac_micros = micro_inds(sac_dir(micro_inds) > 3*pi/8 | sac_dir(micro_inds) < -3*pi/8);

par_micros = [intersect(hor_sac_micros,hor_stim_micros); intersect(ver_sac_micros,ver_stim_micros)];
orth_micros = [intersect(hor_sac_micros,ver_stim_micros); intersect(ver_sac_micros,hor_stim_micros)];

hsac_hstim_micros = intersect(hor_sac_micros,hor_stim_micros);
hsac_vstim_micros = intersect(hor_sac_micros,ver_stim_micros);
vsac_hstim_micros = intersect(ver_sac_micros,hor_stim_micros);
vsac_vstim_micros = intersect(ver_sac_micros,ver_stim_micros);

%%
sm_sig = (0.01/dt);
for cc = 1:96
    all_binned_spks_sm(:,cc) = jmm_smooth_1d_cor(all_binned_spks(:,cc),sm_sig);
end

all_binned_spks_norm = nan(size(all_binned_spks));
for ee = 1:length(use_expts)
    cur_expt_inds = find(all_exptvec == ee);
    expt_mean_rates(ee,:) = mean(all_binned_spks_sm(cur_expt_inds,:))/dt;
    all_binned_spks_norm(cur_expt_inds,:) = bsxfun(@rdivide,all_binned_spks_sm(cur_expt_inds,:),mean(all_binned_spks_sm(cur_expt_inds,:)));
end

ov_avg_rates = mean(expt_mean_rates);
ov_cv_rates = std(expt_mean_rates)./mean(expt_mean_rates);
use_units = find(ov_avg_rates >= 1 & ov_cv_rates <= 0.1);

%%
close all
backlag = round(0.25/dt);
forwardlag = round(0.45/dt);
nboot = 100;

use_trial_ids = all_trialvec; use_trial_ids(~all_used_inds) = 0;
% use_trial_ids = [];
[avg_vstim,lags,std_vstim] = get_event_trig_avg(all_binned_spks_sm,interp_saccade_inds(ver_stim_micros),backlag,forwardlag,nboot,use_trial_ids);
[avg_hstim,lags,std_hstim] = get_event_trig_avg(all_binned_spks_sm,interp_saccade_inds(hor_stim_micros),backlag,forwardlag,nboot,use_trial_ids);

[avg_vstim_vsac,lags,std_vstim_vsac] = get_event_trig_avg(all_binned_spks_sm,interp_saccade_inds(vsac_vstim_micros),backlag,forwardlag,nboot,use_trial_ids);
[avg_vstim_hsac,lags,std_vstim_hsac] = get_event_trig_avg(all_binned_spks_sm,interp_saccade_inds(hsac_vstim_micros),backlag,forwardlag,nboot,use_trial_ids);
[avg_hstim_vsac,lags,std_hstim_vsac] = get_event_trig_avg(all_binned_spks_sm,interp_saccade_inds(vsac_hstim_micros),backlag,forwardlag,nboot,use_trial_ids);
[avg_hstim_hsac,lags,std_hstim_hsac] = get_event_trig_avg(all_binned_spks_sm,interp_saccade_inds(hsac_hstim_micros),backlag,forwardlag,nboot,use_trial_ids);

%%
close all
for cc = 1:length(use_units)
    fprintf('Unit %d %d/%d\n',use_units(cc),cc,length(use_units));
    subplot(2,1,1)
    shadedErrorBar(lags*dt,avg_vstim(:,use_units(cc)),std_vstim(:,use_units(cc)));
    hold on
    shadedErrorBar(lags*dt,avg_hstim(:,use_units(cc)),std_hstim(:,use_units(cc)),{'color','r'});

%     subplot(2,1,2)
%     shadedErrorBar(lags*dt,avg_par(:,use_units(cc)),std_par(:,use_units(cc)))
%     hold on
%     shadedErrorBar(lags*dt,avg_orth(:,use_units(cc)),std_orth(:,use_units(cc)),{'color','r'})
    subplot(2,1,2)
    shadedErrorBar(lags*dt,avg_hstim_hsac(:,use_units(cc)),std_hstim_hsac(:,use_units(cc)))
    hold on
    shadedErrorBar(lags*dt,avg_vstim_hsac(:,use_units(cc)),std_vstim_hsac(:,use_units(cc)),{'color','r'})
    shadedErrorBar(lags*dt,avg_vstim_vsac(:,use_units(cc)),std_vstim_vsac(:,use_units(cc)),{'color','g'})
    shadedErrorBar(lags*dt,avg_hstim_vsac(:,use_units(cc)),std_hstim_vsac(:,use_units(cc)),{'color','b'})

    pause
    clf
end
%%
back_shift = round(0.2/dt);
forw_shift = round(0.3/dt);
flen = back_shift + forw_shift+1;
stim_params = NIMcreate_stim_params(flen,dt,1,1,length(use_expts)-1);

micro_Xvec = zeros(size(all_stim_times));
micro_Xvec(interp_saccade_inds(micro_inds) - back_shift) = 1;
micro_Xmat = create_time_embedding(micro_Xvec,stim_params);

reg_params = NIMcreate_reg_params('lambda_d2T',500);


silent = 1;
cc = 1;
Robs = all_binned_spks(all_used_inds,cc);
fit0(cc) = NIMinitialize_model(stim_params,1,{'lin'},reg_params); %initialize NIM
fit0(cc) = NIMfit_filters(fit0(cc),Robs,micro_Xmat(all_used_inds,:),Xexpt(all_used_inds,:),[],silent); %fit stimulus filters

