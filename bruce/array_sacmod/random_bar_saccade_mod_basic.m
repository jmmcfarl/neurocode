clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G081';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/G081/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
if ~strcmp(Expt_name,'G081')
    load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
end

anal_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

recompute_prepdata = 1;
prep_data_name = [anal_dir '/sac_prepdata.mat'];
recompute_sacdata = 0;
sac_data_name = [anal_dir '/sac_sacdata.mat'];

if exist('./CellList.mat','file')
    load('CellList.mat');
else
    disp('No CellList found.');
    CellList = [];
end

%% PARSE TRIAL DATA STRUCTURES

stim_fs = 100; %in Hz
dt = 0.005;
beg_buffer = round(0.15/dt);
end_buffer = round(0.15/dt);

if strcmp(Expt_name,'G081')
    trial_dur = 2;
else
    trial_dur = 4;
end

%%
if strcmp(Expt_name,'G081')
    for i = 1:length(Expts)
        if strcmp(Expts{i}.Header.expname,'grating.OpXseRC') | strcmp(Expts{i}.Header.expname,'grating.OpRC')
            is_bar_expt(i) = 1;
        else
            is_bar_expt(i) = 0;
        end
        
        if strcmp(Expts{i}.Stimvals.Bs,'image')
            expt_imback(i) = 1;
        else
            expt_imback(i) = 0;
        end
        
        expt_sim_sacs(i) = Expts{i}.Stimvals.ijump;
        expt_bar_ori(i) = Expts{i}.Stimvals.or;
    end
    expt_has_ds = (expt_sim_sacs == 0);
    expt_bar_ori(expt_bar_ori == -45) = 135;
    expt_binoc = zeros(size(expt_bar_ori));
    expt_imback = expt_imback';
    cur_expt_set = find(is_bar_expt);
else
    include_expts = {'rls.Fa', 'rls.FaXimi'};
    for ii = 1:length(Expts)
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    end
    expt_has_ds(isnan(expt_has_ds)) = 0;
    expt_has_ds(expt_has_ds == -1) = 0;
    expt_binoc(isnan(expt_binoc)) = 0;
    % cur_expt_set = find(included_type & expt_has_ds' == 1);
    cur_expt_set = find(included_type & ~expt_binoc');
    
    if strcmp(Expt_name,'G087')
        cur_expt_set(cur_expt_set == 15) = []; %only 6 trials and causes problems
    end
end

if strcmp(Expt_name,'G087')
    cur_expt_set(cur_expt_set == 15) = []; %only 6 trials and causes problems
end

%% LOAD ONE CLUSTER FILE TO IDENTIFY MARKED UNITS
fname = sprintf('Expt%dClusterTimes.mat',cur_expt_set(1));
load(fname);
for cc = 1:96
   if isfield(Clusters{cc},'marked')
       clust_mark(cc) = Clusters{cc}.marked;
   else
       clust_mark(cc) = 0;
   end
end

%% PARSE UNITS
single_units = find(any(CellList(:,:,1) > 0));
n_sus = length(single_units);
good_mus = find(clust_mark == 4);
bad_chs = find(clust_mark == 3);

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
    all_binned_mua = [];
    all_binned_sua = [];
    all_exptvec = [];
    all_trialvec = [];
    all_trial_start_times = [];
    all_trial_end_times = [];
    for ee = 1:length(cur_expt_set);
        fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
        cur_expt = cur_expt_set(ee);
        fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
        load(fname);
        fname = sprintf('Expt%dClusterTimesDetails.mat',cur_expt);
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
        
        all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
        all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
        n_trials = length(use_trials);
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start'/1e4;
            if length(cur_stim_times) == 1
                cur_stim_times = trial_start_times(use_trials(tt)):1/stim_fs:trial_end_times(use_trials(tt));
            end
            cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
            cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
            cur_binned_mua = nan(length(cur_t_axis),96);
            for cc = 1:96
                if ~ismember(cc,single_units)
%                     cur_hist = histc(Clusters{cc}.times,cur_t_edges);
%                     %use clustered MUA
                    cur_hist = histc(ClusterDetails{cc}.t,cur_t_edges); %use all thresholded MUA
                else
                    cur_hist = histc(ClusterDetails{cc}.t(ClusterDetails{cc}.clst==1),cur_t_edges);
                end
                cur_binned_mua(:,cc) = cur_hist(1:end-1);
            end
             cur_binned_sua = nan(length(cur_t_axis),n_sus);
            for cc = 1:n_sus
                cur_hist = histc(Clusters{cc}.times,cur_t_edges);
                cur_binned_sua(:,cc) = cur_hist(1:end-1);
            end
           
            cur_used_inds = true(length(cur_t_axis),1);
            cur_used_inds(1:beg_buffer) = false;
            cur_used_inds((end-end_buffer+1):end) = false;
            
            cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
            
            all_stim_times = [all_stim_times; cur_stim_times'];
            all_t_axis = [all_t_axis; cur_t_axis];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_used_inds = [all_used_inds; cur_used_inds];
            all_binned_mua = [all_binned_mua; cur_binned_mua];
            all_binned_sua = [all_binned_sua; cur_binned_sua];
            all_exptvec = [all_exptvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        end
        trial_cnt = trial_cnt + n_trials;
    end
    all_binned_mua = int16(all_binned_mua);
    all_binned_sua = int16(all_binned_sua);
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
    interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);
    
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
    pre_inds = saccade_inds - sac_dbuff;
    pre_inds(pre_inds < 1) = 1;
    sac_pre_pos = all_eye_vals(pre_inds,:);
    post_inds = saccade_inds + sac_dbuff;
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

%% PICK OUT SACCADES FOR ANALYSIS
sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
saccades(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];

sac_stop_times = [saccades(:).stop_time];
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;

sac_peak_times = [saccades(:).peak_time];
interp_sac_peak_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_peak_times));
interp_sac_peak_inds(isnan(interp_sac_peak_inds)) = 1;

max_dur = 0.1;
sac_durs = [saccades(:).duration];
is_blink = sac_durs > max_dur;

isis = [saccades(:).isi];
is_sacburst = false(length(saccades),1);
is_sacburst(isis(1:end-1) < 0.15 | isis(2:end) < 0.15) = true;

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps < 1;

sac_post_Lx = [saccades(:).post_Lx];
sac_post_Ly = [saccades(:).post_Ly];
sac_pre_Lx = [saccades(:).pre_Lx];
sac_pre_Ly = [saccades(:).pre_Ly];
sac_deltaX = sac_post_Lx - sac_pre_Lx;
sac_deltaY = sac_post_Ly - sac_pre_Ly;

if strcmp(Expt_name,'G081')
    sac_dirs = [0 45 90 135];
    sim_sac_times = [0.7 1.4];
    sac_thresh = 0.5;
else
    sac_dirs = [0 90];
    sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
    sac_thresh = 1;
end
sac_start_expts = all_exptvec(interp_sac_start_inds);
delta_sacpar = nan(size(saccades));
for bb = 1:length(sac_dirs)   
    cur_bar_expts = find(expt_bar_ori(cur_expt_set) == sac_dirs(bb));
    cur_sac_set = find(ismember(sac_start_expts,cur_bar_expts));
    
    cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));    
    delta_sacpar(cur_sac_set) = cur_delta_sac_par;
end
is_gsac = delta_sacpar' > sac_thresh;

sim_sac_expts = find(~expt_has_ds(cur_expt_set));

all_sim_sacs = [];
sim_expt_inds = find(ismember(all_exptvec,sim_sac_expts));
sim_sacs = cell(length(sim_sac_times),1);
for ii = 1:length(sim_sac_times)
    sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
        all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
    all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
end

%define which saccades to use
used_msacs = find(is_micro & ~is_blink);
used_gsacs = find(is_gsac & ~is_blink);

imback_gs_expts = find(expt_has_ds(cur_expt_set) & expt_imback(cur_expt_set)');
grayback_gs_expts = find(expt_has_ds(cur_expt_set) & ~expt_imback(cur_expt_set)' & ~expt_binoc(cur_expt_set));
hori_expts = find(expt_bar_ori(cur_expt_set) == 0);
ver_expts = find(expt_bar_ori(cur_expt_set) == 90);

trial_start_inds = find(all_tsince_start < dt);
simsac_trial_starts = trial_start_inds(ismember(all_exptvec(trial_start_inds),sim_sac_expts));
grayback_trial_starts = trial_start_inds(ismember(all_exptvec(trial_start_inds),grayback_gs_expts));
imback_trial_starts = trial_start_inds(ismember(all_exptvec(trial_start_inds),imback_gs_expts));

used_msacs_grayback = used_msacs(ismember(all_exptvec(interp_sac_start_inds(used_msacs)),grayback_gs_expts));
used_msacs_imback = used_msacs(ismember(all_exptvec(interp_sac_start_inds(used_msacs)),imback_gs_expts));
used_gsacs_grayback = used_gsacs(ismember(all_exptvec(interp_sac_start_inds(used_gsacs)),grayback_gs_expts));
used_gsacs_imback = used_gsacs(ismember(all_exptvec(interp_sac_start_inds(used_gsacs)),imback_gs_expts));

used_msacs_hori = used_msacs(ismember(all_exptvec(interp_sac_start_inds(used_msacs)),hori_expts));
used_msacs_ver = used_msacs(ismember(all_exptvec(interp_sac_start_inds(used_msacs)),ver_expts));
used_gsacs_hori = used_gsacs(ismember(all_exptvec(interp_sac_start_inds(used_gsacs)),hori_expts));
used_gsacs_ver = used_gsacs(ismember(all_exptvec(interp_sac_start_inds(used_gsacs)),ver_expts));

simsacs_hori = all_sim_sacs(ismember(all_exptvec(all_sim_sacs),hori_expts));
simsacs_ver = all_sim_sacs(ismember(all_exptvec(all_sim_sacs),ver_expts));

%%
all_binned_mua_sm = nan(size(all_binned_mua));
all_binned_sua_sm = nan(size(all_binned_sua));
sm_sig = (0.0075/dt);
for cc = 1:96
    all_binned_mua_sm(:,cc) = jmm_smooth_1d_cor(double(all_binned_mua(:,cc)),sm_sig);
end
for cc = 1:n_sus
    all_binned_sua_sm(:,cc) = jmm_smooth_1d_cor(double(all_binned_sua(:,cc)),sm_sig);
end
all_binned_mua_norm = nan(size(all_binned_mua));
all_binned_sua_norm = nan(size(all_binned_sua));
for ee = 1:length(cur_expt_set)
    cur_expt_inds = find(all_exptvec == ee);
    expt_mean_muarates(ee,:) = mean(all_binned_mua_sm(cur_expt_inds,:))/dt;
    all_binned_mua_norm(cur_expt_inds,:) = bsxfun(@rdivide,all_binned_mua_sm(cur_expt_inds,:),...
        mean(all_binned_mua_sm(cur_expt_inds,:)));
    expt_mean_suarates(ee,:) = mean(all_binned_sua_sm(cur_expt_inds,:))/dt;
    all_binned_sua_norm(cur_expt_inds,:) = bsxfun(@rdivide,all_binned_sua_sm(cur_expt_inds,:),...
        mean(all_binned_sua_sm(cur_expt_inds,:)));
end
ov_avg_mua_rates = mean(expt_mean_muarates);
ov_avg_sua_rates = mean(expt_mean_suarates);

binned_micro_sacs = hist(interp_sac_start_inds(used_msacs),1:length(all_t_axis));
binned_msac_sm = jmm_smooth_1d_cor(binned_micro_sacs,sm_sig);

%% WTIHOUT MEAN-RATE NORMALIZATION
all_binned_mua_sm = nan(size(all_binned_mua));
all_binned_sua_sm = nan(size(all_binned_sua));
sm_sig = (0.01/dt);
for cc = 1:96
    all_binned_mua_sm(:,cc) = jmm_smooth_1d_cor(double(all_binned_mua(:,cc)),sm_sig);
end
for cc = 1:n_sus
    all_binned_sua_sm(:,cc) = jmm_smooth_1d_cor(double(all_binned_sua(:,cc)),sm_sig);
end

%%
nboot = [];

%trial averages
% [trial_avgs,trial_lags] = get_event_trig_avg(all_binned_mua_norm,trial_start_inds,0,round(4/dt));
[simsac_trial_avgs,trial_lags] = get_event_trig_avg(all_binned_mua_norm,simsac_trial_starts,0,round(4/dt));
[grayback_trial_avgs,trial_lags] = get_event_trig_avg(all_binned_mua_norm,grayback_trial_starts,0,round(4/dt));
[imback_trial_avgs,trial_lags] = get_event_trig_avg(all_binned_mua_norm,imback_trial_starts,0,round(4/dt));

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

%general averages
[msac_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_msacs),backlag,forwardlag,nboot,all_trialvec);
[gsac_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_gsacs),backlag,forwardlag,nboot,all_trialvec);
[simsac_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,all_sim_sacs,backlag,forwardlag,nboot,all_trialvec);

%microsacs excluding bursts
non_burst_msacs = used_msacs(~is_sacburst(used_msacs));
burst_msacs = used_msacs(is_sacburst(used_msacs));
[nonburst_msac_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(non_burst_msacs),backlag,forwardlag,nboot,all_trialvec);
[burst_msac_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(burst_msacs),backlag,forwardlag,nboot,all_trialvec);


%stop- and peak-triggered averages
[msac_stop_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_stop_inds(used_msacs),backlag,forwardlag,nboot,all_trialvec);
[gsac_stop_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_stop_inds(used_gsacs),backlag,forwardlag,nboot,all_trialvec);
[msac_peak_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_peak_inds(used_msacs),backlag,forwardlag,nboot,all_trialvec);
[gsac_peak_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_peak_inds(used_gsacs),backlag,forwardlag,nboot,all_trialvec);

%background dependent
[msac_gray_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_msacs_grayback),backlag,forwardlag,nboot,all_trialvec);
[gsac_gray_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_gsacs_grayback),backlag,forwardlag,nboot,all_trialvec);
[msac_im_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_msacs_imback),backlag,forwardlag,nboot,all_trialvec);
[gsac_im_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_gsacs_imback),backlag,forwardlag,nboot,all_trialvec);

%saccade direction dependent
[msac_ver_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_msacs_ver),backlag,forwardlag,nboot,all_trialvec);
[gsac_ver_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_gsacs_ver),backlag,forwardlag,nboot,all_trialvec);
[simsac_ver_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,simsacs_ver,backlag,forwardlag,nboot,all_trialvec);
[msac_hor_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_msacs_hori),backlag,forwardlag,nboot,all_trialvec);
[gsac_hor_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,interp_sac_start_inds(used_gsacs_hori),backlag,forwardlag,nboot,all_trialvec);
[simsac_hor_avgs,lags] = get_event_trig_avg(all_binned_mua_norm,simsacs_hori,backlag,forwardlag,nboot,all_trialvec);

% for ii = 1:length(sim_sac_times)
%   [simsac_avg{ii},lags] = get_event_trig_avg(all_binned_spks_norm,sim_sacs{ii},backlag,forwardlag,nboot);
% end  

% eye_backlag = round(0.4/dt);
% eye_forwardlag = round(0.9/dt);
% [all_simsac_avg_msac,eye_lags] = get_event_trig_avg(binned_msac_sm(:),all_sim_sacs,eye_backlag,eye_forwardlag,nboot,all_trialvec);
% [all_gsac_avg_msac,eye_lags] = get_event_trig_avg(binned_msac_sm(:),interp_sac_start_inds(used_gsacs),eye_backlag,eye_forwardlag,nboot,all_trialvec);
% [all_msac_avg_msac,eye_lags] = get_event_trig_avg(binned_msac_sm(:),interp_sac_start_inds(used_msacs),eye_backlag,eye_forwardlag,nboot,all_trialvec);
% [all_simsac_avg_eyespeed,eye_lags] = get_event_trig_avg(interp_eye_speed,all_sim_sacs,eye_backlag,eye_forwardlag,nboot,all_trialvec);
% [all_gsac_avg_eyespeed,eye_lags] = get_event_trig_avg(interp_eye_speed,interp_sac_start_inds(used_gsacs),eye_backlag,eye_forwardlag,nboot,all_trialvec);
% [all_msac_avg_eyespeed,eye_lags] = get_event_trig_avg(interp_eye_speed,interp_sac_start_inds(used_msacs),eye_backlag,eye_forwardlag,nboot,all_trialvec);
% for ii = 1:length(sim_sac_times)
%   [simsac_avg_msac{ii},lags] = get_event_trig_avg(binned_msac_sm(:),sim_sacs{ii},backlag,forwardlag,nboot);
% end  

%% PRINT ENSEMBLE AVERAGES
to_print = 0;

% exclude_units = [bad_chs];
exclude_units = [bad_chs single_units];
used_units = setdiff(1:96,exclude_units);

fname1 = [anal_dir '/sacmod_all_units'];
if to_print == 1
    figure('visible','off');
else
    figure
end
shadedErrorBar(lags*dt,mean(simsac_avgs(:,used_units),2),std(simsac_avgs(:,used_units),[],2)/sqrt(length(used_units)));
hold on
shadedErrorBar(lags*dt,mean(gsac_avgs(:,used_units),2),std(gsac_avgs(:,used_units),[],2)/sqrt(length(used_units)),{'color','r'});
shadedErrorBar(lags*dt,mean(msac_avgs(:,used_units),2),std(msac_avgs(:,used_units),[],2)/sqrt(length(used_units)),{'color','b'});
plot(lags*dt,mean(gsac_stop_avgs(:,used_units),2),'r--');
plot(lags*dt,mean(msac_stop_avgs(:,used_units),2),'b--');
% plot(lags*dt,mean(gsac_peak_avgs,2),'r--');
% plot(lags*dt,mean(msac_peak_avgs,2),'b--');
xlim([-backlag*dt forwardlag*dt])
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
if to_print == 1
    print(fname1,'-dpsc');
    close all
end

fname2 = [anal_dir '/sacmod_all_units_dircompare'];
if to_print == 1
    figure('visible','off');
else
    figure
end
subplot(3,1,1)
shadedErrorBar(lags*dt,mean(gsac_ver_avgs(:,used_units),2),std(gsac_ver_avgs(:,used_units),[],2)/sqrt(length(used_units)));
hold on
shadedErrorBar(lags*dt,mean(gsac_hor_avgs(:,used_units),2),std(gsac_hor_avgs(:,used_units),[],2)/sqrt(length(used_units)),{'color','r'});
xlim([-backlag*dt forwardlag*dt])
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
subplot(3,1,2)
shadedErrorBar(lags*dt,mean(msac_ver_avgs(:,used_units),2),std(msac_ver_avgs(:,used_units),[],2)/sqrt(length(used_units)));
hold on
shadedErrorBar(lags*dt,mean(msac_hor_avgs(:,used_units),2),std(msac_hor_avgs(:,used_units),[],2)/sqrt(length(used_units)),{'color','r'});
xlim([-backlag*dt forwardlag*dt])
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
subplot(3,1,3)
shadedErrorBar(lags*dt,mean(simsac_hor_avgs(:,used_units),2),std(simsac_hor_avgs(:,used_units),[],2)/sqrt(length(used_units)));
hold on
shadedErrorBar(lags*dt,mean(simsac_ver_avgs(:,used_units),2),std(simsac_ver_avgs(:,used_units),[],2)/sqrt(length(used_units)),{'color','r'});
xlim([-backlag*dt forwardlag*dt])
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
if to_print == 1
    fillPage(gcf,'papersize',[8 15]);
    print(fname2,'-dpsc');
    close all
end

fname3 = [anal_dir '/sacmod_all_units_backcomp'];
if to_print == 1
    figure('visible','off');
else
    figure
end
subplot(2,1,1)
shadedErrorBar(lags*dt,mean(gsac_gray_avgs(:,used_units),2),std(gsac_gray_avgs(:,used_units),[],2)/sqrt(length(used_units)));
hold on
shadedErrorBar(lags*dt,mean(gsac_im_avgs(:,used_units),2),std(gsac_im_avgs(:,used_units),[],2)/sqrt(length(used_units)),{'color','r'});
xlim([-backlag*dt forwardlag*dt])
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
subplot(2,1,2)
shadedErrorBar(lags*dt,mean(msac_gray_avgs(:,used_units),2),std(msac_gray_avgs(:,used_units),[],2)/sqrt(length(used_units)));
hold on
shadedErrorBar(lags*dt,mean(msac_im_avgs(:,used_units),2),std(msac_im_avgs(:,used_units),[],2)/sqrt(length(used_units)),{'color','r'});
xlim([-backlag*dt forwardlag*dt])
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
if to_print == 1
    fillPage(gcf,'papersize',[8 12]);
    print(fname3,'-dpsc');
    close all
end

%%
backlag = round(0.25/dt);
forwardlag = round(0.55/dt);
nboot = 50;
[all_simsac_avgs,lags,all_simsac_std] = get_event_trig_avg(all_binned_mua_sm,all_sim_sacs,backlag,forwardlag,nboot,all_trialvec);
[all_gsac_avgs,lags,all_gsac_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_gsacs),backlag,forwardlag,nboot,all_trialvec);
[all_msac_avgs,lags,all_msac_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_msacs),backlag,forwardlag,nboot,all_trialvec);

[all_gsac_stop_avgs,lags] = get_event_trig_avg(all_binned_mua_sm,interp_sac_stop_inds(used_gsacs),backlag,forwardlag,0,all_trialvec);
[all_msac_stop_avgs,lags] = get_event_trig_avg(all_binned_mua_sm,interp_sac_stop_inds(used_msacs),backlag,forwardlag,0,all_trialvec);
[all_gsac_peak_avgs,lags] = get_event_trig_avg(all_binned_mua_sm,interp_sac_peak_inds(used_gsacs),backlag,forwardlag,0,all_trialvec);
[all_msac_peak_avgs,lags] = get_event_trig_avg(all_binned_mua_sm,interp_sac_peak_inds(used_msacs),backlag,forwardlag,0,all_trialvec);

[all_msac_grayback_avgs,lags,all_msac_grayback_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_msacs_grayback),backlag,forwardlag,nboot,all_trialvec);
[all_msac_imback_avgs,lags,all_msac_imback_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_msacs_imback),backlag,forwardlag,nboot,all_trialvec);
[all_gsac_grayback_avgs,lags,all_gsac_grayback_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_gsacs_grayback),backlag,forwardlag,nboot,all_trialvec);
[all_gsac_imback_avgs,lags,all_gsac_imback_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_gsacs_imback),backlag,forwardlag,nboot,all_trialvec);

[all_msac_ver_avgs,lags,all_msac_ver_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_msacs_ver),backlag,forwardlag,nboot,all_trialvec);
[all_gsac_ver_avgs,lags,all_gsac_ver_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_gsacs_ver),backlag,forwardlag,nboot,all_trialvec);
[all_simsac_ver_avgs,lags,all_simsac_ver_std] = get_event_trig_avg(all_binned_mua_sm,simsacs_ver,backlag,forwardlag,nboot,all_trialvec);
[all_msac_hor_avgs,lags,all_msac_hor_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_msacs_hori),backlag,forwardlag,nboot,all_trialvec);
[all_gsac_hor_avgs,lags,all_gsac_hor_std] = get_event_trig_avg(all_binned_mua_sm,interp_sac_start_inds(used_gsacs_hori),backlag,forwardlag,nboot,all_trialvec);
[all_simsac_hor_avgs,lags,all_simsac_hor_std] = get_event_trig_avg(all_binned_mua_sm,simsacs_hori,backlag,forwardlag,nboot,all_trialvec);

%% SINGLE-UNIT ANALYSIS
backlag = round(0.25/dt);
forwardlag = round(0.55/dt);
nboot = 50;

for cc = 1:length(single_units)
    
    fprintf('Estimating trig avgs for SU %d of %d\n',cc,length(single_units));
    
    use_expts = find(CellList(cur_expt_set,single_units(cc),1) > 0);
    use_inds = find(ismember(all_exptvec,use_expts));
    su_avg_rates(cc) = mean(all_binned_sua(use_inds,cc));
    
    cur_sim_sacs = find(ismember(use_inds,all_sim_sacs));
    cur_gsacs = find(ismember(use_inds,interp_sac_start_inds(used_gsacs)));
    cur_msacs = find(ismember(use_inds,interp_sac_start_inds(used_msacs)));
    
    [su_simsac_navgs(:,cc),lags,su_simsac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_sim_sacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_gsac_navgs(:,cc),lags,su_gsac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_msac_navgs(:,cc),lags,su_msac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));    
    
    [su_simsac_avgs(:,cc),lags,su_simsac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_sim_sacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_gsac_avgs(:,cc),lags,su_gsac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_msac_avgs(:,cc),lags,su_msac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));    
 
    %for horizontal expts
    use_hexpts = use_expts(expt_bar_ori(cur_expt_set(use_expts)) == 0);
    use_inds = find(ismember(all_exptvec,use_hexpts));
    cur_sim_sacs = find(ismember(use_inds,all_sim_sacs));
    cur_gsacs = find(ismember(use_inds,interp_sac_start_inds(used_gsacs)));
    cur_msacs = find(ismember(use_inds,interp_sac_start_inds(used_msacs)));
    
    [su_hor_simsac_avgs(:,cc),lags,su_hor_simsac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_sim_sacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_hor_gsac_avgs(:,cc),lags,su_hor_gsac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_hor_msac_avgs(:,cc),lags,su_hor_msac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));    
    [su_hor_simsac_navgs(:,cc),lags,su_hor_simsac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_sim_sacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_hor_gsac_navgs(:,cc),lags,su_hor_gsac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_hor_msac_navgs(:,cc),lags,su_hor_msac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));    
    
    %for vertical expts
    use_vexpts = use_expts(expt_bar_ori(cur_expt_set(use_expts)) == 90);
    use_inds = find(ismember(all_exptvec,use_vexpts));
    cur_sim_sacs = find(ismember(use_inds,all_sim_sacs));
    cur_gsacs = find(ismember(use_inds,interp_sac_start_inds(used_gsacs)));
    cur_msacs = find(ismember(use_inds,interp_sac_start_inds(used_msacs)));
    
    [su_ver_simsac_avgs(:,cc),lags,su_ver_simsac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_sim_sacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_ver_gsac_avgs(:,cc),lags,su_ver_gsac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_ver_msac_avgs(:,cc),lags,su_ver_msac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_ver_simsac_navgs(:,cc),lags,su_ver_simsac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_sim_sacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_ver_gsac_navgs(:,cc),lags,su_ver_gsac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_ver_msac_navgs(:,cc),lags,su_ver_msac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    
    %for gray-back expts
    use_gbexpts = use_expts(ismember(use_expts,grayback_gs_expts));
    use_inds = find(ismember(all_exptvec,use_gbexpts));
    cur_gsacs = find(ismember(use_inds,interp_sac_start_inds(used_gsacs)));
    cur_msacs = find(ismember(use_inds,interp_sac_start_inds(used_msacs)));
    
    [su_gb_gsac_avgs(:,cc),lags,su_gb_gsac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_gb_msac_avgs(:,cc),lags,su_gb_msac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_gb_gsac_navgs(:,cc),lags,su_gb_gsac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_gb_msac_navgs(:,cc),lags,su_gb_msac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));

    %for image-back expts
    use_ibexpts = use_expts(ismember(use_expts,imback_gs_expts));
    use_inds = find(ismember(all_exptvec,use_ibexpts));
    cur_gsacs = find(ismember(use_inds,interp_sac_start_inds(used_gsacs)));
    cur_msacs = find(ismember(use_inds,interp_sac_start_inds(used_msacs)));
    
    [su_ib_gsac_avgs(:,cc),lags,su_ib_gsac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_ib_msac_avgs(:,cc),lags,su_ib_msac_std(:,cc)] = get_event_trig_avg(all_binned_sua_sm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_ib_gsac_navgs(:,cc),lags,su_ib_gsac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_gsacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    [su_ib_msac_navgs(:,cc),lags,su_ib_msac_nstd(:,cc)] = get_event_trig_avg(all_binned_sua_norm(use_inds,cc),...
        cur_msacs,backlag,forwardlag,nboot,all_trialvec(use_inds));
    
end

if ~isempty(single_units)
    sname = [anal_dir '/su_sacmod_data.mat'];
    save(sname,'su_*','lags','dt');
end

%% EXTRACT SAC-MOD STATS
poss_gsac_enh = find(lags*dt >= 0.05 & lags*dt <= 0.15);
poss_msac_enh = find(lags*dt >= 0.05 & lags*dt <= 0.25);
poss_simsac_enh = find(lags*dt >= 0.05 & lags*dt <= 0.25);
poss_gsac_sup = find(lags*dt >= 0 & lags*dt <= 0.1);
poss_msac_sup = find(lags*dt >= 0 & lags*dt <= 0.1);
poss_simsac_sup = find(lags*dt >= 0 & lags*dt <= 0.2);

avg_rates = mean(all_binned_mua);

[simsac_sup_rate,simsac_sup_time] = min(bsxfun(@rdivide,all_simsac_avgs(poss_simsac_sup,:),avg_rates)); 
simsac_sup_time = lags(simsac_sup_time + poss_simsac_sup(1) - 1)*dt;
[gsac_sup_rate,gsac_sup_time] = min(bsxfun(@rdivide,all_gsac_avgs(poss_gsac_sup,:),avg_rates));
gsac_sup_time = lags(gsac_sup_time + poss_gsac_sup(1) - 1)*dt;
[msac_sup_rate,msac_sup_time] = min(bsxfun(@rdivide,all_msac_avgs(poss_msac_sup,:),avg_rates));
msac_sup_time = lags(msac_sup_time + poss_msac_sup(1) - 1)*dt;

[simsac_enh_rate,simsac_enh_time] = max(bsxfun(@rdivide,all_simsac_avgs(poss_simsac_enh,:),avg_rates));
simsac_enh_time = lags(simsac_enh_time + poss_simsac_enh(1) - 1)*dt;
[gsac_enh_rate,gsac_enh_time] = max(bsxfun(@rdivide,all_gsac_avgs(poss_gsac_enh,:),avg_rates));
gsac_enh_time = lags(gsac_enh_time + poss_gsac_enh(1) - 1)*dt;
[msac_enh_rate,msac_enh_time] = max(bsxfun(@rdivide,all_msac_avgs(poss_msac_enh,:),avg_rates));
msac_enh_time = lags(msac_enh_time + poss_msac_enh(1) - 1)*dt;

%now for SUs
if ~isempty(single_units)
%     [su_simsac_sup_rate,su_simsac_sup_time] = min(bsxfun(@rdivide,su_simsac_avgs(poss_simsac_sup,:),su_avg_rates));
    [su_simsac_sup_rate,su_simsac_sup_time] = min(su_simsac_navgs(poss_simsac_sup,:));
    su_simsac_sup_time = lags(su_simsac_sup_time + poss_simsac_sup(1) - 1)*dt;
%     [su_gsac_sup_rate,su_gsac_sup_time] = min(bsxfun(@rdivide,su_gsac_avgs(poss_gsac_sup,:),su_avg_rates));
    [su_gsac_sup_rate,su_gsac_sup_time] = min(su_gsac_navgs(poss_gsac_sup,:));
    su_gsac_sup_time = lags(su_gsac_sup_time + poss_gsac_sup(1) - 1)*dt;
%     [su_msac_sup_rate,su_msac_sup_time] = min(bsxfun(@rdivide,su_msac_avgs(poss_msac_sup,:),su_avg_rates));
    [su_msac_sup_rate,su_msac_sup_time] = min(su_msac_navgs(poss_msac_sup,:));
    su_msac_sup_time = lags(su_msac_sup_time + poss_msac_sup(1) - 1)*dt;
    
%     [su_simsac_enh_rate,su_simsac_enh_time] = max(bsxfun(@rdivide,su_simsac_avgs(poss_simsac_enh,:),su_avg_rates));
    [su_simsac_enh_rate,su_simsac_enh_time] = max(su_simsac_navgs(poss_simsac_enh,:));
    su_simsac_enh_time = lags(su_simsac_enh_time + poss_simsac_enh(1) - 1)*dt;
%     [su_gsac_enh_rate,su_gsac_enh_time] = max(bsxfun(@rdivide,su_gsac_avgs(poss_gsac_enh,:),su_avg_rates));
    [su_gsac_enh_rate,su_gsac_enh_time] = max(su_gsac_navgs(poss_gsac_enh,:));
    su_gsac_enh_time = lags(su_gsac_enh_time + poss_gsac_enh(1) - 1)*dt;
%     [su_msac_enh_rate,su_msac_enh_time] = max(bsxfun(@rdivide,su_msac_avgs(poss_msac_enh,:),su_avg_rates));
    [su_msac_enh_rate,su_msac_enh_time] = max(su_msac_navgs(poss_msac_enh,:));
    su_msac_enh_time = lags(su_msac_enh_time + poss_msac_enh(1) - 1)*dt;
    
    sname = [anal_dir '/su_sacmod_data.mat'];
    save(sname,'su_*','lags','dt');
end
%%
sname = [anal_dir '/sacmod_data.mat'];
save(sname,'all_*sac*','lags','dt','*_avgs','*_avg_rates');

%% PRINT FIGURES FOR EACH INDIVIDUAL UNIT (APPENDING TO SINGLE FILE)
close all
for cc = 1:96
    cur_mean_rate = mean(all_binned_mua(tr_inds,cc))/dt;
    fprintf('Cell %d of %d\n',cc,96);
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    hold on
    shadedErrorBar(lags*dt,all_gsac_avgs(:,cc)/dt,all_gsac_std(:,cc)/dt,{'color','r'});
    shadedErrorBar(lags*dt,all_msac_avgs(:,cc)/dt,all_msac_std(:,cc)/dt,{'color','b'});
    shadedErrorBar(lags*dt,all_simsac_avgs(:,cc)/dt,all_simsac_std(:,cc)/dt);
    plot(lags*dt,all_gsac_stop_avgs(:,cc)/dt,'r--');
    plot(lags*dt,all_msac_stop_avgs(:,cc)/dt,'b--');
%     plot(lags*dt,all_gsac_peak_avgs(:,cc)/dt,'r--');
%     plot(lags*dt,all_msac_peak_avgs(:,cc)/dt,'b--');
    xlim([-backlag*dt forwardlag*dt])
    xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Firing rate (Hz)','fontsize',16);
    title(sprintf('Unit%d',cc));
    if to_print == 1
        print(fname1,'-append','-dpsc');
        close all
    else
        pause
        close all
    end
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    subplot(3,1,1)
    shadedErrorBar(lags*dt,all_gsac_ver_avgs(:,cc)/dt,all_gsac_ver_std(:,cc)/dt);
    hold on
    shadedErrorBar(lags*dt,all_gsac_hor_avgs(:,cc)/dt,all_gsac_hor_std(:,cc)/dt,{'color','r'});
xlim([-backlag*dt forwardlag*dt])
    xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Relative firing rate','fontsize',16);
    title(sprintf('Unit%d Guided sacs',cc));
    subplot(3,1,2)
    shadedErrorBar(lags*dt,all_msac_ver_avgs(:,cc)/dt,all_msac_ver_std(:,cc)/dt);
    hold on
    shadedErrorBar(lags*dt,all_msac_hor_avgs(:,cc)/dt,all_msac_hor_std(:,cc)/dt,{'color','r'});
 xlim([-backlag*dt forwardlag*dt])
   xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Relative firing rate','fontsize',16);
    title(sprintf('Unit%d Micro sacs',cc));
    subplot(3,1,3)
    shadedErrorBar(lags*dt,all_simsac_hor_avgs(:,cc)/dt,all_simsac_hor_std(:,cc)/dt);
    hold on
    shadedErrorBar(lags*dt,all_simsac_ver_avgs(:,cc)/dt,all_simsac_hor_std(:,cc)/dt,{'color','r'});
xlim([-backlag*dt forwardlag*dt])
    xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Relative firing rate','fontsize',16);
    title(sprintf('Unit%d Sim sacs',cc));
    if to_print == 1
        fillPage(gcf,'papersize',[8 15]);
        print(fname2,'-append','-dpsc');
        close all
    end
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    subplot(2,1,1)
    shadedErrorBar(lags*dt,all_gsac_grayback_avgs(:,cc)/dt,all_gsac_grayback_std(:,cc)/dt);
    hold on
    shadedErrorBar(lags*dt,all_gsac_imback_avgs(:,cc)/dt,all_gsac_imback_std(:,cc)/dt,{'color','r'});
xlim([-backlag*dt forwardlag*dt])
    xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Relative firing rate','fontsize',16);
    title(sprintf('Unit%d Guided sacs',cc));
    subplot(2,1,2)
    shadedErrorBar(lags*dt,all_msac_grayback_avgs(:,cc)/dt,all_msac_grayback_std(:,cc)/dt);
    hold on
    shadedErrorBar(lags*dt,all_msac_imback_avgs(:,cc)/dt,all_msac_imback_std(:,cc)/dt,{'color','r'});
xlim([-backlag*dt forwardlag*dt])
    xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Relative firing rate','fontsize',16);
    title(sprintf('Unit%d Micro sacs',cc));
    if to_print == 1
        fillPage(gcf,'papersize',[8 12]);
        print(fname3,'-append','-dpsc');
        close all
    end
    
end

%% PRINT FIGURES FOR SUs (APPENDING TO SINGLE FILE)
close all
to_print = 0;
for cc = 1:length(single_units)
    cur_mean_rate = su_avg_rates(cc)/dt;
    fprintf('SU %d of %d\n',cc,length(single_units));
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    subplot(2,1,1)
    hold on
    shadedErrorBar(lags*dt,su_gsac_avgs(:,cc)/dt,su_gsac_std(:,cc)/dt,{'color','r'});
    shadedErrorBar(lags*dt,su_msac_avgs(:,cc)/dt,su_msac_std(:,cc)/dt,{'color','b'});
    shadedErrorBar(lags*dt,su_simsac_avgs(:,cc)/dt,su_simsac_std(:,cc)/dt);
    xlim([-backlag*dt forwardlag*dt])
    xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Firing rate (Hz)','fontsize',16);
    title(sprintf('SU %d',single_units(cc)));
    subplot(2,1,2)
    hold on
    shadedErrorBar(lags*dt,su_gsac_navgs(:,cc),su_gsac_nstd(:,cc),{'color','r'});
    shadedErrorBar(lags*dt,su_msac_navgs(:,cc),su_msac_nstd(:,cc),{'color','b'});
    shadedErrorBar(lags*dt,su_simsac_navgs(:,cc),su_simsac_nstd(:,cc));
    xlim([-backlag*dt forwardlag*dt])
    xl = xlim();yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Firing rate (Hz)','fontsize',16);
    title(sprintf('SU %d',single_units(cc)));
    if to_print == 1
        print(fname1,'-append','-dpsc');
        close all
    end
    
%     if to_print == 1
%         figure('visible','off');
%     else
%         figure
%     end
%     subplot(3,1,1)
%     shadedErrorBar(lags*dt,su_ver_gsac_avgs(:,cc)/dt,su_ver_gsac_std(:,cc)/dt);
%     hold on
%     shadedErrorBar(lags*dt,su_hor_gsac_avgs(:,cc)/dt,su_hor_gsac_std(:,cc)/dt,{'color','r'});
%     xlim([-backlag*dt forwardlag*dt])
%     xl = xlim();yl = ylim();
%     line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
%     xlabel('Time since saccade onset (s)','fontsize',16);
%     ylabel('Relative firing rate','fontsize',16);
%     title(sprintf('SU %d Guided sacs',single_units(cc)));
%     subplot(3,1,2)
%     shadedErrorBar(lags*dt,su_ver_msac_avgs(:,cc)/dt,su_ver_msac_std(:,cc)/dt);
%     hold on
%     shadedErrorBar(lags*dt,su_hor_msac_avgs(:,cc)/dt,su_hor_msac_std(:,cc)/dt,{'color','r'});
%     xlim([-backlag*dt forwardlag*dt])
%     xl = xlim();yl = ylim();
%     line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
%     xlabel('Time since saccade onset (s)','fontsize',16);
%     ylabel('Relative firing rate','fontsize',16);
%     title(sprintf('SU %d Micro sacs',single_units(cc)));
%     subplot(3,1,3)
%     shadedErrorBar(lags*dt,su_ver_simsac_avgs(:,cc)/dt,su_ver_simsac_std(:,cc)/dt);
%     hold on
%     shadedErrorBar(lags*dt,su_hor_simsac_avgs(:,cc)/dt,su_hor_simsac_std(:,cc)/dt,{'color','r'});
%     xlim([-backlag*dt forwardlag*dt])
%     xl = xlim();yl = ylim();
%     line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
%     xlabel('Time since saccade onset (s)','fontsize',16);
%     ylabel('Relative firing rate','fontsize',16);
%     title(sprintf('SU %d Sim sacs',single_units(cc)));
%     if to_print == 1
%         fillPage(gcf,'papersize',[8 15]);
%         print(fname2,'-append','-dpsc');
%         close all
%     end
%     
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    subplot(2,1,1)
    shadedErrorBar(lags*dt,su_gb_gsac_avgs(:,cc)/dt,su_gb_gsac_std(:,cc)/dt);
    hold on
    shadedErrorBar(lags*dt,su_ib_gsac_avgs(:,cc)/dt,su_ib_gsac_std(:,cc)/dt,{'color','r'});
    xlim([-backlag*dt forwardlag*dt])
    xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Relative firing rate','fontsize',16);
    title(sprintf('SU %d Guided sacs',single_units(cc)));
    subplot(2,1,2)
    shadedErrorBar(lags*dt,su_gb_msac_avgs(:,cc)/dt,su_gb_msac_std(:,cc)/dt);
    hold on
    shadedErrorBar(lags*dt,su_ib_msac_avgs(:,cc)/dt,su_ib_msac_std(:,cc)/dt,{'color','r'});
    xlim([-backlag*dt forwardlag*dt])
    xl = xlim();yl = ylim();
    line(xl,cur_mean_rate([1 1]),'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',16);
    ylabel('Relative firing rate','fontsize',16);
    title(sprintf('SU %d Micro sacs',single_units(cc)));
    if to_print == 1
        fillPage(gcf,'papersize',[8 12]);
        print(fname3,'-append','-dpsc');
        close all
    end
    
    if to_print == 0
        pause
        close all
    end
end

%%
sac_backlag = 0.25;
sac_forlag = 0.55;
sac_bin_width = 1;
sac_binspace = sac_bin_width*dt;
sac_bin_edges = -(sac_backlag-sac_binspace/2):sac_binspace:(sac_forlag+sac_binspace/2);
sac_bin_cents = 0.5*sac_bin_edges(1:end-1) + 0.5*sac_bin_edges(2:end);
n_sac_bins = length(sac_bin_cents);

gsac_inds = interp_sac_start_inds(used_gsacs);
msac_inds = interp_sac_start_inds(used_msacs);
burst_msac_inds = interp_sac_start_inds(burst_msacs);
noburst_msac_inds = interp_sac_start_inds(non_burst_msacs);

trial_simsac_mat = zeros(length(all_t_axis),n_sac_bins);
trial_gsac_mat = zeros(length(all_t_axis),n_sac_bins);
trial_msac_mat = zeros(length(all_t_axis),n_sac_bins);
trial_burst_msac_mat = zeros(length(all_t_axis),n_sac_bins);
trial_noburst_msac_mat = zeros(length(all_t_axis),n_sac_bins);
for i = 1:n_sac_bins
    for cc = 1:sac_bin_width
        cur_inds = all_sim_sacs + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
        cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
        trial_simsac_mat(cur_inds,i) = 1;
        cur_inds = gsac_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
        cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
        trial_gsac_mat(cur_inds,i) = 1;
        cur_inds = msac_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
        cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
        trial_msac_mat(cur_inds,i) = 1;
%         cur_inds = burst_msac_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
%         cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
%         trial_burst_msac_mat(cur_inds,i) = 1;
%         cur_inds = noburst_msac_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
%         cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
%         trial_noburst_msac_mat(cur_inds,i) = 1;
    end
end

%% MODEL BASED SACCADE FILTER ESTIMATION
L2_params = create_L2_params([],[1 n_sac_bins],n_sac_bins);
L2_params = create_L2_params(L2_params,n_sac_bins + [1 n_sac_bins],n_sac_bins);
L2_params = create_L2_params(L2_params,2*n_sac_bins + [1 n_sac_bins],n_sac_bins);
cur_Xmat = [trial_gsac_mat(tr_inds,:) trial_msac_mat(tr_inds,:) trial_simsac_mat(tr_inds,:)];
for cc = 1:96
    fprintf('Estimating optimal kernel for unit %d of %d\n',cc,96);
    Robs = all_binned_mua(tr_inds,cc);
    [fitmod] = regGLM_fit(cur_Xmat,Robs,L2_params,ones(length(L2_params),1)*1000,[],[],1);
    gsac_kern(cc,:) = fitmod.K((1:n_sac_bins));
    msac_kern(cc,:) = fitmod.K((1:n_sac_bins) + n_sac_bins);
    simsac_kern(cc,:) = fitmod.K((1:n_sac_bins) + 2*n_sac_bins);
    %     [fitmod] = regGLM_fit(trial_burst_msac_mat(tr_inds,:),Robs,L2_params,2500,[],[],1);
    %     burst_msac_kern(cc,:) = fitmod.K;
    %     [fitmod] = regGLM_fit(trial_noburst_msac_mat(tr_inds,:),Robs,L2_params,2500,[],[],1);
    %     noburst_msac_kern(cc,:) = fitmod.K;
end

for cc = 1:96
    subplot(2,1,1)
    plot(sac_bin_cents,msac_kern(cc,:))
    hold on
    plot(sac_bin_cents,gsac_kern(cc,:),'r')
    plot(sac_bin_cents,simsac_kern(cc,:),'k')
    %     plot(sac_bin_cents,burst_msac_kern(cc,:),'k')
    %     plot(sac_bin_cents,noburst_msac_kern(cc,:),'r')
    xlim([-sac_backlag sac_forlag])
    xl = xlim(); line(xl,[0 0],'color','k','linestyle','--');
    subplot(2,1,2)
    plot(lags*dt,msac_avgs(:,cc))
    hold on
    plot(lags*dt,gsac_avgs(:,cc),'r')
    plot(lags*dt,simsac_avgs(:,cc),'k')
    xlim([-sac_backlag sac_forlag])
    xl = xlim(); line(xl,[1 1],'color','k','linestyle','--');
    pause
    clf
end
%%
% figure;hold on
% shadedErrorBar(lags*dt,mean(msac_gray_avgs,2),std(msac_gray_avgs,[],2)/sqrt(96),{'color','b'});
% shadedErrorBar(lags*dt,mean(msac_im_avgs,2),std(msac_im_avgs,[],2)/sqrt(96),{'color','r'});
%     xl = xlim();yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
% 
% figure;hold on
% shadedErrorBar(lags*dt,mean(gsac_gray_avgs,2),std(gsac_gray_avgs,[],2)/sqrt(96),{'color','b'});
% shadedErrorBar(lags*dt,mean(gsac_im_avgs,2),std(gsac_im_avgs,[],2)/sqrt(96),{'color','r'});
%     xl = xlim();yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');


%%
% figure
% shadedErrorBar(trial_lags*dt,mean(simsac_trial_avgs,2),std(simsac_trial_avgs,[],2)/sqrt(96));
% hold on
% shadedErrorBar(trial_lags*dt,mean(grayback_trial_avgs,2),std(grayback_trial_avgs,[],2)/sqrt(96),{'color','r'});
% shadedErrorBar(trial_lags*dt,mean(imback_trial_avgs,2),std(imback_trial_avgs,[],2)/sqrt(96),{'color','b'});

%%
% figure
% shadedErrorBar(lags*dt,mean(all_simsac_avgs,2),std(all_simsac_avgs,[],2)/sqrt(96));
% hold on
% cmap = jet(length(sim_sac_times));
% for ii = 1:length(sim_sac_times)
% shadedErrorBar(lags*dt,mean(simsac_avg{ii},2),std(simsac_avg{ii},[],2)/sqrt(96),{'color',cmap(ii,:)});
% end

use_units = find(ov_avg_rates > 2);

% figure
hold on
shadedErrorBar(lags*dt,mean(gsac_avgs(:,use_units),2),std(gsac_avgs(:,use_units),[],2)/sqrt(length(use_units)),{'color','k'});
shadedErrorBar(lags*dt,mean(msac_avgs(:,use_units),2),std(msac_avgs(:,use_units),[],2)/sqrt(length(use_units)),{'color','g'});

xlabel('Time since saccade onset (s)','fontsize',16)
ylabel('Relative firing rate','fontsize',16)
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
shg
