clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G088';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/G081/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/par_orth_sacs'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

recompute_prepdata = 1;
prep_data_name = [anal_dir '/sac_prepdata.mat'];
recompute_sacdata = 1;
sac_data_name = [anal_dir '/sac_sacdata.mat'];

%% PARSE TRIAL DATA STRUCTURES

stim_fs = 100; %in Hz
dt = 0.005;
beg_buffer = round(0.15/dt);
end_buffer = round(0.15/dt);

%% PARSE TRIAL DATA STRUCTURES
target_expt_name = {'rls.seXFaRC','rls.seXFaXFrRC','rls.orXFaXsl'};
is_target_type = false(length(Expts),1);
expt_frame_dur = nan(length(Expts),1);
expt_sac_dir = nan(length(Expts),1);
expt_bar_ori = nan(length(Expts),1);
for i = 1:length(Expts)
    expt_names{i} = Expts{i}.Header.expname;
    if any(strcmp(expt_names{i},target_expt_name))
        is_target_type(i) = true;
    end
    expt_frame_dur(i) = Expts{i}.Stimvals.Fr;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
    expt_sac_dir(i) = mod(Expts{i}.Stimvals.Fa,180);
end

cur_expt_set = find(is_target_type);
if strcmp(Expt_name,'G088')
    cur_expt_set(cur_expt_set == 22) = []; %only 4 trials and causes problems
end

%% COMPUTE TRIAL DATA
cd(data_dir);

if recompute_prepdata == 1 || ~exist(prep_data_name,'file')
 %%
 fprintf('Computing prep data\n');
    trial_cnt = 0;
    
    all_stim_times = [];
    all_t_axis = [];
    all_tsince_start = [];
    all_used_inds = false(0);
    all_phase_Xmat = [];
    all_binned_spks = [];
    all_exptvec = [];
    all_bar_or = [];
    all_frame_dur = [];
    all_sac_dir = [];
    all_trialvec = [];
    all_trial_start_times = [];
    all_trial_end_times = [];
    for ee = 1:length(cur_expt_set);
        fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
        cur_expt = cur_expt_set(ee);
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
        
        use_trials = find(trial_durs >= 0.5);
        
        if isfield(Expts{cur_expt}.Trials(1),'Fr')
            trial_Fr = [Expts{cur_expt}.Trials(:).Fr];
            all_frame_dur = [all_frame_dur; trial_Fr(use_trials)'];
        elseif isfield(Expts{cur_expt}.Trials(1),'sl')
            trial_sl = [Expts{cur_expt}.Trials(:).sl];
            trial_sl(trial_sl == 0) = 1;
            all_frame_dur = [all_frame_dur; trial_sl(use_trials)'];            
        else
            all_frame_dur = [all_frame_dur; ones(length(use_trials),1)*expt_frame_dur(cur_expt)];
        end
        if isfield(Expts{cur_expt}.Trials(1),'or')
            trial_or = [Expts{cur_expt}.Trials(:).or];
            all_bar_or = [all_bar_or; trial_or(use_trials)'];
        else
            all_bar_or = [all_bar_or; ones(length(use_trials),1)*expt_bar_ori(cur_expt)];            
        end
        if isfield(Expts{cur_expt}.Trials(1),'Fa')
            trial_Fa = [Expts{cur_expt}.Trials(:).Fa];
            all_sac_dir = [all_sac_dir; trial_Fa(use_trials)'];
        else
            all_sac_dir = [all_sac_dir; ones(length(use_trials),1)*expt_sac_dir(cur_expt)];                        
        end
        
        all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
        all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
        n_trials = length(use_trials);
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start/1e4;
            if length(cur_stim_times) == 1
                cur_stim_times = trial_start_times(use_trials(tt)):1/stim_fs:trial_end_times(use_trials(tt));
            end
            cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
            cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
            cur_binned_spks = nan(length(cur_t_axis),96);
            for cc = 1:96
                cur_hist = histc(Clusters{cc}.times,cur_t_edges);
                cur_binned_spks(:,cc) = cur_hist(1:end-1);
                %             temp = convert_to_spikebins(cur_hist(1:end-1));
            end
            
            cur_used_inds = true(length(cur_t_axis),1);
            cur_used_inds(1:beg_buffer) = false;
            cur_used_inds((end-end_buffer+1):end) = false;
            
            cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
            
            all_stim_times = [all_stim_times; cur_stim_times(:)];
            all_t_axis = [all_t_axis; cur_t_axis];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_used_inds = [all_used_inds; cur_used_inds];
            all_binned_spks = [all_binned_spks; cur_binned_spks];
            all_exptvec = [all_exptvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        end
        trial_cnt = trial_cnt + n_trials;
    end
    save(prep_data_name,'all_*');
else
    fprintf('Loading pre-computed prep data\n');
    load(prep_data_name);
end

%%
[c,ia,ic] = unique(all_trialvec);
n_trials = length(ia);

xv_frac = 0;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_set))';

tr_inds(~all_used_inds(tr_inds)) = [];
xv_inds(~all_used_inds(xv_inds)) = [];

Xexpt = zeros(length(all_stim_times),length(cur_expt_set)-1);
for i = 1:length(cur_expt_set)-1
    cur_set = find(all_exptvec==i);
    Xexpt(cur_set,i) = 1;
end

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
if recompute_sacdata == 1 || ~exist(sac_data_name,'file')
    
    emfile = ['jbe' Expt_name '.em.mat'];
    load(emfile);
    
    all_eye_vals = [];
    all_eye_speed = [];
    all_eye_ts = [];
    all_eye_exptvec = [];
    eye_smooth = 3;
    for ee = 1:length(cur_expt_set);
        fprintf('Loading eye data for expt %d of %d\n',ee,length(cur_expt_set));
        
        cur_set = find(all_exptvec==ee);
        if ~isempty(cur_set)
        [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));
        
        eye_dt = Expt.Header.CRrates(1);
        eye_fs = 1/eye_dt;
        lEyeXY = eye_vals_interp(:,1:2);
        rEyeXY = eye_vals_interp(:,3:4);
        
        %slight smoothing before computing speed
        sm_avg_eyepos = lEyeXY; eye_vel = lEyeXY; %initialization
        sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),eye_smooth);
        sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),eye_smooth);
        eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
        eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
        
        eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
        all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
        all_eye_speed = [all_eye_speed; eye_speed];
        all_eye_ts = [all_eye_ts; eye_ts_interp'];
        all_eye_exptvec = [all_eye_exptvec; ee*ones(size(eye_speed))];
        end
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
    
    interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis);
    
    %%
    sac_thresh = 10;
    peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
    saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);
    
    peri_thresh = 3; %threshold eye speed for defining saccade boundary inds
    thresh_cross_up = 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh);
    thresh_cross_down = 1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh);
    sac_start_inds = nan(size(saccade_inds));
    sac_stop_inds = nan(size(saccade_inds));
    for ii = 1:length(saccade_inds)
        next_tc = find(thresh_cross_down > saccade_inds(ii),1,'first');
        if ~isempty(next_tc)
            sac_stop_inds(ii) = thresh_cross_down(next_tc);
        end
        prev_tc = find(thresh_cross_up < saccade_inds(ii),1,'last');
        if ~isempty(prev_tc)
            sac_start_inds(ii) = thresh_cross_up(prev_tc);
        end
        
    end
    
    %get rid of double-peaks
    min_isi = 0.05; max_isi = Inf;
    isis = [Inf; diff(sac_start_inds)]/eye_fs;
    bad_isis = (isis < min_isi | isis > max_isi);
    bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
    saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];
    
    saccade_times = all_eye_ts(saccade_inds);
    sac_start_times = all_eye_ts(sac_start_inds);
    sac_stop_times = all_eye_ts(sac_stop_inds);
    sac_durs = sac_stop_times - sac_start_times;
    
    sac_dbuff = round(0.005/eye_dt);
    pre_inds = sac_start_inds - sac_dbuff;
    pre_inds(pre_inds < 1) = 1;
    sac_pre_pos = all_eye_vals(pre_inds,:);
    post_inds = sac_stop_inds + sac_dbuff;
    post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
    sac_post_pos = all_eye_vals(post_inds,:);
    
    %use only left-eye signal here
    sac_delta_pos = sac_post_pos(:,1:2) - sac_pre_pos(:,1:2);
    sac_amps = sqrt(sum(sac_delta_pos.^2,2));
    sac_dirs = atan2(sac_delta_pos(:,2),sac_delta_pos(:,1));
    
    temp = ones(length(saccade_times),1);
    saccades = struct('peak_time',mat2cell(saccade_times,temp),'start_time',mat2cell(sac_start_times,temp),...
        'stop_time',mat2cell(sac_stop_times,temp),'isi',mat2cell(isis,temp),...
        'duration',mat2cell(sac_durs,temp),'amplitude',mat2cell(sac_amps,temp),'direction',mat2cell(sac_dirs,temp),...
        'pre_Lx',mat2cell(sac_pre_pos(:,1),temp),'post_Lx',mat2cell(sac_post_pos(:,1),temp),...
        'pre_Ly',mat2cell(sac_pre_pos(:,2),temp),'post_Ly',mat2cell(sac_post_pos(:,2),temp),...
        'pre_Rx',mat2cell(sac_pre_pos(:,3),temp),'post_Rx',mat2cell(sac_post_pos(:,3),temp),...
        'pre_Ry',mat2cell(sac_pre_pos(:,4),temp),'post_Ry',mat2cell(sac_post_pos(:,4),temp));
    
    save(sac_data_name,'saccades');
else
    fprintf('Loading pre-computed sac data\n');
    load(sac_data_name);
end

%%

%% PICK OUT SACCADES FOR ANALYSIS
sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
saccades(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];

max_dur = 0.1;
sac_durs = [saccades(:).duration];
is_blink = sac_durs > max_dur;

isis = [saccades(:).isi];
is_sacburst = (isis(1:end-1) < 0.15 | isis(2:end) < 0.15);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps < 1;

sac_post_Lx = [saccades(:).post_Lx];
sac_post_Ly = [saccades(:).post_Ly];
sac_pre_Lx = [saccades(:).pre_Lx];
sac_pre_Ly = [saccades(:).pre_Ly];
sac_deltaX = sac_post_Lx - sac_pre_Lx;
sac_deltaY = sac_post_Ly - sac_pre_Ly;

all_sac_dir = mod(all_sac_dir,180);

sac_dirs = [0 90];
sac_expts = all_exptvec(interp_sac_start_inds);
sac_trials = all_trialvec(interp_sac_start_inds);
delta_sacpar = nan(size(saccades));
delta_sacorth = nan(size(saccades));
for bb = 1:2    
    cur_sac_trials = find(all_sac_dir == sac_dirs(bb));
    cur_sac_set = find(ismember(sac_trials,cur_sac_trials));
    
    cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));    
    delta_sacpar(cur_sac_set) = cur_delta_sac_par;
    cur_delta_sac_orth = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180+90) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180+90));    
    delta_sacorth(cur_sac_set) = cur_delta_sac_orth;
end
is_gsac = delta_sacpar' > 2 & delta_sacpar' < 4;

%define which saccades to use
used_msacs = find(is_micro & ~is_blink);
used_gsacs = find(is_gsac & ~is_blink);
msac_inds = interp_sac_start_inds(used_msacs);
gsac_inds = interp_sac_start_inds(used_gsacs);

p100_trials = find(all_frame_dur(all_trialvec) == 1 & all_sac_dir(all_trialvec) == all_bar_or(all_trialvec));
p30_trials = find(all_frame_dur(all_trialvec) == 3 & all_sac_dir(all_trialvec) == all_bar_or(all_trialvec));
o100_trials = find(all_frame_dur(all_trialvec) == 1 & all_sac_dir(all_trialvec) ~= all_bar_or(all_trialvec));
o30_trials = find(all_frame_dur(all_trialvec) == 3 & all_sac_dir(all_trialvec) ~= all_bar_or(all_trialvec));

p100_gsacs = gsac_inds(ismember(gsac_inds,p100_trials));
p30_gsacs = gsac_inds(ismember(gsac_inds,p30_trials));
o100_gsacs = gsac_inds(ismember(gsac_inds,o100_trials));
o30_gsacs = gsac_inds(ismember(gsac_inds,o30_trials));

p100_msacs = msac_inds(ismember(msac_inds,p100_trials));
p30_msacs = msac_inds(ismember(msac_inds,p30_trials));
o100_msacs = msac_inds(ismember(msac_inds,o100_trials));
o30_msacs = msac_inds(ismember(msac_inds,o30_trials));


%%
all_binned_spks_sm = nan(size(all_binned_spks));
sm_sig = (0.005/dt);
for cc = 1:96
    all_binned_spks_sm(:,cc) = jmm_smooth_1d_cor(all_binned_spks(:,cc),sm_sig);
end
all_binned_spks_norm = nan(size(all_binned_spks));
for ee = 1:length(cur_expt_set)
    cur_expt_inds = find(all_exptvec == ee);
    expt_mean_rates(ee,:) = mean(all_binned_spks_sm(cur_expt_inds,:))/dt;
    all_binned_spks_norm(cur_expt_inds,:) = bsxfun(@rdivide,all_binned_spks_sm(cur_expt_inds,:),mean(all_binned_spks_sm(cur_expt_inds,:)));
end
ov_avg_rates = mean(expt_mean_rates);

%%
nboot = [];
backlag = round(0.2/dt);
forwardlag = round(0.5/dt);

[p100_gsac_pavg,lags] = get_event_trig_avg(all_binned_spks_norm,p100_gsacs,backlag,forwardlag,nboot);
[p30_gsac_pavg,lags] = get_event_trig_avg(all_binned_spks_norm,p30_gsacs,backlag,forwardlag,nboot);
[o100_gsac_pavg,lags] = get_event_trig_avg(all_binned_spks_norm,o100_gsacs,backlag,forwardlag,nboot);
[o30_gsac_pavg,lags] = get_event_trig_avg(all_binned_spks_norm,o30_gsacs,backlag,forwardlag,nboot);

[p100_msac_pavg,lags] = get_event_trig_avg(all_binned_spks_norm,p100_msacs,backlag,forwardlag,nboot);
[p30_msac_pavg,lags] = get_event_trig_avg(all_binned_spks_norm,p30_msacs,backlag,forwardlag,nboot);
[o100_msac_pavg,lags] = get_event_trig_avg(all_binned_spks_norm,o100_msacs,backlag,forwardlag,nboot);
[o30_msac_pavg,lags] = get_event_trig_avg(all_binned_spks_norm,o30_msacs,backlag,forwardlag,nboot);

%%
close all
to_print = 1;
fname = [anal_dir '/parorth_all_units'];

if to_print == 1
    figure('visible','off');
end
subplot(2,1,1)
shadedErrorBar(lags*dt,mean(p100_gsac_pavg,2),std(p100_gsac_pavg,[],2)/sqrt(96));
hold on
shadedErrorBar(lags*dt,mean(p30_gsac_pavg,2),std(p30_gsac_pavg,[],2)/sqrt(96),{'color','r'});
shadedErrorBar(lags*dt,mean(o100_gsac_pavg,2),std(o100_gsac_pavg,[],2)/sqrt(96),{'color','b'});
shadedErrorBar(lags*dt,mean(o30_gsac_pavg,2),std(o30_gsac_pavg,[],2)/sqrt(96),{'color','g'});
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
subplot(2,1,2)
shadedErrorBar(lags*dt,mean(p100_msac_pavg,2),std(p100_msac_pavg,[],2)/sqrt(96));
hold on
shadedErrorBar(lags*dt,mean(p30_msac_pavg,2),std(p30_msac_pavg,[],2)/sqrt(96),{'color','r'});
shadedErrorBar(lags*dt,mean(o100_msac_pavg,2),std(o100_msac_pavg,[],2)/sqrt(96),{'color','b'});
shadedErrorBar(lags*dt,mean(o30_msac_pavg,2),std(o30_msac_pavg,[],2)/sqrt(96),{'color','g'});
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
fillPage(gcf,'papersize',[6 10]);
if to_print == 1
    print(fname,'-dpsc');
    close all
end

%%
all_binned_spks_sm = nan(size(all_binned_spks));
sm_sig = (0.01/dt);
for cc = 1:96
    all_binned_spks_sm(:,cc) = jmm_smooth_1d_cor(all_binned_spks(:,cc),sm_sig);
end
all_binned_spks_norm = all_binned_spks_sm;

backlag = round(0.2/dt);
forwardlag = round(0.5/dt);
nboot = 100;

[p100_gsac_avg,lags,p100_gsac_std] = get_event_trig_avg(all_binned_spks_norm,p100_gsacs,backlag,forwardlag,nboot);
[p30_gsac_avg,lags,p30_gsac_std] = get_event_trig_avg(all_binned_spks_norm,p30_gsacs,backlag,forwardlag,nboot);
[o100_gsac_avg,lags,o100_gsac_std] = get_event_trig_avg(all_binned_spks_norm,o100_gsacs,backlag,forwardlag,nboot);
[o30_gsac_avg,lags,o30_gsac_std] = get_event_trig_avg(all_binned_spks_norm,o30_gsacs,backlag,forwardlag,nboot);

[p100_msac_avg,lags,p100_msac_std] = get_event_trig_avg(all_binned_spks_norm,p100_msacs,backlag,forwardlag,nboot);
[p30_msac_avg,lags,p30_msac_std] = get_event_trig_avg(all_binned_spks_norm,p30_msacs,backlag,forwardlag,nboot);
[o100_msac_avg,lags,o100_msac_std] = get_event_trig_avg(all_binned_spks_norm,o100_msacs,backlag,forwardlag,nboot);
[o30_msac_avg,lags,o30_msac_std] = get_event_trig_avg(all_binned_spks_norm,o30_msacs,backlag,forwardlag,nboot);

%%
sname = [anal_dir '/parorth_data.mat'];
save(sname,'p100*','p30*','lags','dt','o100*','o30*','ov_avg_rates');

%%
close all
for cc = 1:96
    cur_mean_rate = mean(all_binned_spks(tr_inds,cc))/dt;
    fprintf('Cell %d of %d\n',cc,96);
    if to_print == 1
        figure('visible','off');
    end
subplot(2,1,1)
shadedErrorBar(lags*dt,p100_gsac_avg(:,cc)/dt,p100_gsac_std(:,cc)/dt);
hold on
shadedErrorBar(lags*dt,p30_gsac_avg(:,cc)/dt,p30_gsac_std(:,cc)/dt,{'color','r'});
shadedErrorBar(lags*dt,o100_gsac_avg(:,cc)/dt,o100_gsac_std(:,cc)/dt,{'color','b'});
shadedErrorBar(lags*dt,o30_gsac_avg(:,cc)/dt,o30_gsac_std(:,cc)/dt,{'color','g'});
xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
    title(sprintf('Unit%d',cc));
subplot(2,1,2)
shadedErrorBar(lags*dt,p100_msac_avg(:,cc)/dt,p100_msac_std(:,cc)/dt);
hold on
shadedErrorBar(lags*dt,p30_msac_avg(:,cc)/dt,p30_msac_std(:,cc)/dt,{'color','r'});
shadedErrorBar(lags*dt,o100_msac_avg(:,cc)/dt,o100_msac_std(:,cc)/dt,{'color','b'});
shadedErrorBar(lags*dt,o30_msac_avg(:,cc)/dt,o30_msac_std(:,cc)/dt,{'color','g'});
xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
fillPage(gcf,'papersize',[6 10]);
    if to_print == 1
        print(fname,'-append','-dpsc');
        close all
    else
        pause
        clf
    end
end
