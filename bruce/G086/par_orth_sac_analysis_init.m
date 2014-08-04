clear all
% close all

addpath('~/James_scripts/bruce/G081/')

cd ~/Data/bruce/G086/
load jbeG086Expts.mat

stim_fs = 100; %in Hz
use_units = 1:96;

Fs = 200;
new_dt = .005;

%% PARSE TRIAL DATA STRUCTURES
target_expt_name = ['rls.seXFaRC'];
is_target_type = false(length(Expts),1);
expt_frame_dur = nan(length(Expts),1);
expt_sac_dir = nan(length(Expts),1);
for i = 1:length(Expts)
    if strcmp(Expts{i}.Header.expname,target_expt_name)
        is_target_type(i) = true;
    end
    expt_frame_dur(i) = Expts{i}.Stimvals.Fr;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
    expt_sac_dir(i) = mod(Expts{i}.Stimvals.Fa,180);
end



%% DETERMINE SET OF EXPERIMENTS TO USE
%expts with X deg bars and any back (including sim sacs)
cur_expt_set = find(is_target_type);

%% COMPUTE TRIAL DATA
all_stim_times = [];
all_t_axis = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_rel_stimes = [];
all_rel_etimes = [];
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

        cur_t_edges = [(Expts{cur_expt}.Trials(tt).Start(1)/1e4):1/Fs:(Expts{cur_expt}.Trials(tt).End(end)/1e4)];
        cur_t_cents = cur_t_edges(1:end-1) + 2/Fs;
        cur_binned_spks = nan(length(cur_t_cents),length(use_units));
        for cc = 1:length(use_units)
            cur_hist = histc(Clusters{use_units(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_t_axis = [all_t_axis; cur_t_cents(:)];
        all_rel_stimes = [all_rel_stimes; cur_t_cents(:)- trial_start_times(tt)];
        all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_t_cents(:)];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(length(cur_t_cents),1)*ee];
        all_trialvec = [all_trialvec; ones(length(cur_t_cents),1)*tt];
    end
    
end


%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
load ./jbeG086.em.mat
all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
for ee = 1:length(cur_expt_set);
    fprintf('Loading eye data for expt %d of %d\n',ee,length(cur_expt_set));
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
all_big_sacs = [];
sac_oris = [0 90];
for bb = 1:2
    
    par_eye_pos = interp_eye_vals(:,1)*cos(sac_oris(bb)*pi/180) - interp_eye_vals(:,2)*sin(sac_oris(bb)*pi/180);
    bad_pts = find(abs(par_eye_pos) > 1); %cheap way of detecting blinks because the right eye signal was unreliable even for that at times

    cur_bar_expts = find(expt_sac_dir(cur_expt_set) == sac_oris(bb));
    cur_inds = find(ismember(all_exptvec,cur_bar_expts));
    cur_sac_set = find(ismember(interp_saccade_inds,cur_inds));

    pre_pos = par_eye_pos(sac_starts(cur_sac_set));
    post_pos = par_eye_pos(sac_stops(cur_sac_set));

    sac_dist = abs(post_pos-pre_pos);
    big_sacs = find(sac_dist' >= 2);
    msacs = find(sac_dist' < 1); 
    
    all_big_sacs = [all_big_sacs; cur_sac_set(big_sacs)];
    all_msacs = [all_msacs; (cur_sac_set(msacs))];    
    
end
all_msacs = sort(all_msacs);
all_big_sacs = sort(all_big_sacs);

%%
min_isi = 0.05;
all_big_sacs(isis(all_big_sacs) < min_isi) = [];
bad = find(isis(all_msacs) < min_isi);
all_msacs(bad) = [];

all_msac_inds = interp_saccade_inds(all_msacs);
all_bigsac_inds = interp_saccade_inds(all_big_sacs);


%%
sm_sig = (Fs*0.005);
for cc = 1:96
    all_binned_spks_sm(:,cc) = jmm_smooth_1d_cor(all_binned_spks(:,cc),sm_sig);
end

%%
all_binned_spks_norm = nan(size(all_binned_spks));
for ee = 1:length(cur_expt_set)
    cur_expt_inds = find(all_exptvec == ee);
    expt_mean_rates(ee,:) = mean(all_binned_spks_sm(cur_expt_inds,:))/new_dt;
    all_binned_spks_norm(cur_expt_inds,:) = bsxfun(@rdivide,all_binned_spks_sm(cur_expt_inds,:),mean(all_binned_spks_sm(cur_expt_inds,:)));
end

%%
ov_avg_rates = mean(expt_mean_rates);
ov_cv_rates = std(expt_mean_rates)./mean(expt_mean_rates);
use_units = find(ov_avg_rates >= 1 & ov_cv_rates <= 0.1);
% use_units = 1:96;

%%
backlag = round(0.1*Fs);
forwardlag = round(0.3*Fs);

par_100 = find(expt_frame_dur(cur_expt_set) == 1 & expt_sac_dir(cur_expt_set) == 90);
par_30 = find(expt_frame_dur(cur_expt_set) == 3 & expt_sac_dir(cur_expt_set) == 90);
orth_100 = find(expt_frame_dur(cur_expt_set) == 1 & expt_sac_dir(cur_expt_set) == 0);
orth_30 = find(expt_frame_dur(cur_expt_set) == 3 & expt_sac_dir(cur_expt_set) == 0);

parallel_bsacs_100 = all_bigsac_inds(ismember(all_exptvec(all_bigsac_inds),par_100));
parallel_bsacs_30 = all_bigsac_inds(ismember(all_exptvec(all_bigsac_inds),par_30));
orth_bsacs_100 = all_bigsac_inds(ismember(all_exptvec(all_bigsac_inds),orth_100));
orth_bsacs_30 = all_bigsac_inds(ismember(all_exptvec(all_bigsac_inds),orth_30));

parallel_msacs_100 = all_msac_inds(ismember(all_exptvec(all_msac_inds),par_100));
parallel_msacs_30 = all_msac_inds(ismember(all_exptvec(all_msac_inds),par_30));
orth_msacs_100 = all_msac_inds(ismember(all_exptvec(all_msac_inds),orth_100));
orth_msacs_30 = all_msac_inds(ismember(all_exptvec(all_msac_inds),orth_30));

%%
[parallel_bsac100_trig_avgs,lags] = get_event_trig_avg(all_binned_spks_norm,parallel_bsacs_100,backlag,forwardlag);
[parallel_bsac30_trig_avgs,lags] = get_event_trig_avg(all_binned_spks_norm,parallel_bsacs_30,backlag,forwardlag);
[orth_bsac100_trig_avgs,lags] = get_event_trig_avg(all_binned_spks_norm,orth_bsacs_100,backlag,forwardlag);
[orth_bsac30_trig_avgs,lags] = get_event_trig_avg(all_binned_spks_norm,orth_bsacs_30,backlag,forwardlag);

[parallel_msac100_trig_avgs,lags] = get_event_trig_avg(all_binned_spks_norm,parallel_msacs_100,backlag,forwardlag);
[parallel_msac30_trig_avgs,lags] = get_event_trig_avg(all_binned_spks_norm,parallel_msacs_30,backlag,forwardlag);
[orth_msac100_trig_avgs,lags] = get_event_trig_avg(all_binned_spks_norm,orth_msacs_100,backlag,forwardlag);
[orth_msac30_trig_avgs,lags] = get_event_trig_avg(all_binned_spks_norm,orth_msacs_30,backlag,forwardlag);

%% BSACS
disp('Computing p100 stats')
n_boot = 100;
% [p100_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm,parallel_bsacs_100,backlag,forwardlag);
[p100_mat,cur_lags] = get_event_trig_mat(all_binned_spks_sm,parallel_bsacs_100,backlag,forwardlag);
p100_mat = reshape(p100_mat,length(parallel_bsacs_100),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,p100_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_p100 = squeeze(mean(boot_avgs))/new_dt;
std_p100 = squeeze(std(boot_avgs))/new_dt;

disp('Computing p30 stats')
n_boot = 100;
% [p30_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm,parallel_bsacs_30,backlag,forwardlag);
[p30_mat,cur_lags] = get_event_trig_mat(all_binned_spks_sm,parallel_bsacs_30,backlag,forwardlag);
p30_mat = reshape(p30_mat,length(parallel_bsacs_30),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,p30_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_p30 = squeeze(mean(boot_avgs))/new_dt;
std_p30 = squeeze(std(boot_avgs))/new_dt;

disp('Computing o100 stats')
n_boot = 100;
% [o100_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm,orth_bsacs_100,backlag,forwardlag);
[o100_mat,cur_lags] = get_event_trig_mat(all_binned_spks_sm,orth_bsacs_100,backlag,forwardlag);
o100_mat = reshape(o100_mat,length(orth_bsacs_100),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,o100_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_o100 = squeeze(mean(boot_avgs))/new_dt;
std_o100 = squeeze(std(boot_avgs))/new_dt;

disp('Computing o30 stats')
n_boot = 100;
% [o30_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm,orth_bsacs_30,backlag,forwardlag);
[o30_mat,cur_lags] = get_event_trig_mat(all_binned_spks_sm,orth_bsacs_30,backlag,forwardlag);
o30_mat = reshape(o30_mat,length(orth_bsacs_30),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,o30_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_o30 = squeeze(mean(boot_avgs))/new_dt;
std_o30 = squeeze(std(boot_avgs))/new_dt;

%% MSACS
disp('Computing p100 stats')
n_boot = 100;
[p100m_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm,parallel_msacs_100,backlag,forwardlag);
p100m_mat = reshape(p100m_mat,length(parallel_msacs_100),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,p100m_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_p100m = squeeze(mean(boot_avgs));
std_p100m = squeeze(std(boot_avgs));

disp('Computing p30 stats')
n_boot = 100;
[p30m_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm,parallel_msacs_30,backlag,forwardlag);
p30m_mat = reshape(p30m_mat,length(parallel_msacs_30),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,p30m_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_p30m = squeeze(mean(boot_avgs));
std_p30m = squeeze(std(boot_avgs));

disp('Computing o100 stats')
n_boot = 100;
[o100m_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm,orth_msacs_100,backlag,forwardlag);
o100m_mat = reshape(o100m_mat,length(orth_msacs_100),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,o100m_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_o100m = squeeze(mean(boot_avgs));
std_o100m = squeeze(std(boot_avgs));

disp('Computing o30 stats')
n_boot = 100;
[o30m_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm,orth_msacs_30,backlag,forwardlag);
o30m_mat = reshape(o30m_mat,length(orth_msacs_30),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,o30m_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_o30m = squeeze(mean(boot_avgs));
std_o30m = squeeze(std(boot_avgs));

%%
out_dir = '~/Analysis/bruce/G086/par_orth_sacs/';
close all
for cc = use_units
    cc
    %     subplot(2,1,1)
    plot(0,ov_avg_rates(cc),'k-');     hold on
    plot(0,ov_avg_rates(cc),'r'); plot(0,ov_avg_rates(cc),'b-');plot(0,ov_avg_rates(cc),'g-');
    legend('Par 100Hz','Par 30Hz','Orth 100Hz','Orth 30Hz');
    
    shadedErrorBar(cur_lags/Fs*1000,avg_p100(:,cc),std_p100(:,cc),{'color','k'})
    shadedErrorBar(cur_lags/Fs*1000,avg_p30(:,cc),std_p30(:,cc),{'color','r'})
    shadedErrorBar(cur_lags/Fs*1000,avg_o100(:,cc),std_o100(:,cc),{'color','b'})
    shadedErrorBar(cur_lags/Fs*1000,avg_o30(:,cc),std_o30(:,cc),{'color','g'})
    xlim([-backlag forwardlag]/Fs*1000)
    grid on
    xlabel('Time from saccade onset (ms)','fontsize',14);
    ylabel('Average rate (Hz)','fontsize',14);
    
    fname = strcat(out_dir,sprintf('Unit%d',cc));
    print(fname,'-dpdf');
    close
    
%     pause
%     clf
end

%%
figure
plot(0,1,'k-');     hold on
plot(0,1,'r'); plot(0,1,'b-');plot(0,1,'g-');
legend('Par 100Hz','Par 30Hz','Orth 100Hz','Orth 30Hz');
shadedErrorBar(lags/Fs,mean(parallel_bsac100_trig_avgs(:,use_units),2),std(parallel_bsac100_trig_avgs(:,use_units),[],2)/sqrt(length(use_units)))
shadedErrorBar(lags/Fs,mean(parallel_bsac30_trig_avgs(:,use_units),2),std(parallel_bsac30_trig_avgs(:,use_units),[],2)/sqrt(length(use_units)),{'color','r'})
shadedErrorBar(lags/Fs,mean(orth_bsac100_trig_avgs(:,use_units),2),std(orth_bsac100_trig_avgs(:,use_units),[],2)/sqrt(length(use_units)),{'color','b'})
shadedErrorBar(lags/Fs,mean(orth_bsac30_trig_avgs(:,use_units),2),std(orth_bsac30_trig_avgs(:,use_units),[],2)/sqrt(length(use_units)),{'color','g'})
    xlim([-backlag forwardlag]/Fs)
    xlabel('Time from saccade onset (ms)','fontsize',14);
    ylabel('Average rate (Hz)','fontsize',14);
grid on
%%
figure
plot(0,1,'k-');     hold on
plot(0,1,'r'); plot(0,1,'b-');plot(0,1,'g-');
legend('Par 100Hz','Par 30Hz','Orth 100Hz','Orth 30Hz');
shadedErrorBar(lags/Fs,mean(parallel_msac100_trig_avgs(:,use_units),2),std(parallel_msac100_trig_avgs(:,use_units),[],2)/sqrt(length(use_units)))
hold on
shadedErrorBar(lags/Fs,mean(parallel_msac30_trig_avgs(:,use_units),2),std(parallel_msac30_trig_avgs(:,use_units),[],2)/sqrt(length(use_units)),{'color','r'})
shadedErrorBar(lags/Fs,mean(orth_msac100_trig_avgs(:,use_units),2),std(orth_msac100_trig_avgs(:,use_units),[],2)/sqrt(length(use_units)),{'color','b'})
shadedErrorBar(lags/Fs,mean(orth_msac30_trig_avgs(:,use_units),2),std(orth_msac30_trig_avgs(:,use_units),[],2)/sqrt(length(use_units)),{'color','g'})
    xlim([-backlag forwardlag]/Fs)
    xlabel('Time from saccade onset (ms)','fontsize',14);
    ylabel('Average rate (Hz)','fontsize',14);
grid on

%%
disp('Computing big sac stats')
n_boot = 100;
[sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_bsac_inds,round(0.4*Fs),round(0.5*Fs));
sac_trg_mat = reshape(sac_trg_mat,length(tr_bsac_inds),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_bsac_trig_rate = squeeze(mean(boot_avgs));
std_bsac_trig_rate = squeeze(std(boot_avgs));


disp('Computing micro sac stats')
[sac_trg_mat,cur_lags] = get_event_trig_mat(all_binned_spks_norm(all_tr_inds,:),tr_msac_inds,round(0.4*Fs),round(0.5*Fs));
sac_trg_mat = reshape(sac_trg_mat,length(tr_msac_inds),96*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 96]);
avg_msac_trig_rate = squeeze(mean(boot_avgs));
std_msac_trig_rate = squeeze(std(boot_avgs));
