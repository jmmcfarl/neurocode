clear all
close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';

Expt_num = [270];

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz

dt = 0.005;

min_trial_dur = 0.5;
beg_buffer = 0.3;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

mua_sm_sig = (0.0075/dt);

n_probes = 24;

clear expt* included_type

Expt_name = sprintf('M%d',Expt_num);
data_dir = [dir_prefix '/Data/bruce/' Expt_name];

cd(data_dir);

load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct

save_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end
save_name = 'trig_avg_data';

%%
if Expt_num > 260
    include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs'};
    for ii = 1:length(Expts)
        if ~isempty(Expts{ii})
            expt_names{ii} = Expts{ii}.Header.expname;
            expt_dds(ii) = Expts{ii}.Stimvals.dd;
            expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
            expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
            expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
            expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
            included_type(ii) = any(strcmp(expt_names{ii},include_expts));
        else
            included_type(ii) = 0;
        end
    end
    expt_binoc = zeros(size(expt_bar_ori));
    %         expt_has_ds(isnan(expt_has_ds)) = 0;
    %         expt_has_ds(expt_has_ds == -1) = 0;
    expt_binoc(isnan(expt_binoc)) = 0;
    cur_block_set = find(included_type & ~expt_binoc & expt_Fr == 1);
    
    if strcmp(Expt_name,'M270')
        cur_block_set(cur_block_set == 5) = [];
    end
    
else
    load('./random_bar_eyedata_ftime.mat','bar_expts');
    load(['lemM' num2str(Expt_num) 'Expts.mat']);
    cur_block_set = bar_expts;
end

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Fs = [];
all_trial_wi = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
    all_spk_times = cell(96,1);
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
for ee = 1:length(cur_block_set);
    cur_block = cur_block_set(ee);
    fprintf('Expt %d of %d\n',ee,length(cur_block_set));
    fname = sprintf('Expt%dClusterTimes.mat',cur_block);
    load(fname);
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times' + cur_toffset);
    end
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    
    [un_ids,id_inds] = unique(trial_ids);
    if length(un_ids) < length(trial_ids)
        fprintf('Warning, repeat trial inds detected!\n');
    end
    
    use_trials = find(trial_durs >= min_trial_dur);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    
    if isfield(Expts{cur_block}.Trials(1),'Fs')
        trial_Fs = [Expts{cur_block}.Trials(:).Fs];
    else
        trial_Fs = nan(1,length(trial_durs));
    end
    %         trial_fs = trial_fs(id_inds);
    all_trial_Fs = cat(1,all_trial_Fs,trial_Fs(use_trials)');
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        if length(cur_stim_times) == 1
            cur_stim_times = trial_start_times(use_trials(tt)):1/stim_fs:trial_end_times(use_trials(tt));
        end
        cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        all_stim_times = [all_stim_times; cur_stim_times' + cur_toffset];
        all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
        all_tsince_start = [all_tsince_start; cur_tsince_start];
        all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    end
    trial_cnt = trial_cnt + n_trials;
    trial_toffset(ee) = all_t_bin_edges(end);
    cur_toffset = trial_toffset(ee);
end

%% SPIKE DATA
    used_inds = find(all_tsince_start >= beg_buffer);

    all_binned_mua = nan(length(all_t_axis),n_probes);
%for only-MU probes
for cc = 1:n_probes
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc},all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end
all_mua_rate = nan(size(all_binned_mua));
for cc = 1:24
    all_mua_rate(:,cc) = jmm_smooth_1d_cor(all_binned_mua(:,cc),mua_sm_sig);
end
all_mua_rate_norm = nan(size(all_binned_mua));

block_mean_muarates = nan(length(cur_block_set),n_probes);
block_n_muaspikes = nan(length(cur_block_set),n_probes);
for ee = 1:length(cur_block_set)
    cur_block_inds = used_inds(all_blockvec(used_inds) == ee);
    block_mean_muarates(ee,:) = mean(all_mua_rate(cur_block_inds,:));
    block_n_muaspikes(ee,:) = sum(all_binned_mua(cur_block_inds,:));
end

mua_avg_rates = nanmean(block_mean_muarates);
mua_tot_nspikes = nansum(block_n_muaspikes);

%normalized firing rates (smoothed)
for ee = 1:length(cur_block_set)
    cur_block_inds = used_inds(all_blockvec(used_inds) == ee);
%     cur_block_inds = (all_blockvec == ee);
    all_mua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_mua_rate(cur_block_inds,:),...
        block_mean_muarates(ee,:));
end

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
all_eye_blockvec = [];
eye_smooth = 3;
for ee = 1:length(cur_block_set);
    
    emfile = sprintf('lem%s.%d.em.mat',Expt_name,cur_block_set(ee));
    load(emfile);
    
    fprintf('Loading ET data for expt %d, block %d of %d\n',Expt_num,ee,length(cur_block_set));
    cur_set = find(all_blockvec==ee);
    if ee > 1
        cur_toffset = trial_toffset(ee-1);
    else
        cur_toffset = 0;
    end
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])) - cur_toffset);
    
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
    all_eye_ts = [all_eye_ts; eye_ts_interp' + cur_toffset];
    all_eye_blockvec = [all_eye_blockvec; ee*ones(size(eye_speed))];
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
sac_ampx = [saccades(:).post_Lx]-[saccades(:).pre_Lx];
sac_ampy = [saccades(:).post_Ly]-[saccades(:).pre_Ly];
sac_dir = Expts{cur_block_set(1)}.Stimvals.or*pi/180;
sac_par_amps = abs(sac_ampx*cos(sac_dir)) + abs(sac_ampy*sin(sac_dir));
if Expt_num > 260
    sac_thresh = 1;
else
    sac_thresh = 0.5;
end
is_gsac = sac_par_amps > sac_thresh;

is_micro = sac_amps < 1;

sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
sac_start_expts = all_blockvec(interp_sac_start_inds);

if Expt_num < 260
    sim_sac_expts = [];
    expt_imback = zeros(1,max(cur_block_set));
else
    for ii = 1:length(cur_block_set)
        cur_eyevar(ii) = var(interp_eye_vals(all_blockvec == ii,2));
    end
    sim_sac_expts = find(cur_eyevar < 1);
end

all_sim_sacs = [];
sim_expt_inds = find(ismember(all_blockvec,sim_sac_expts));
sim_sacs = cell(length(sim_sac_times),1);
for ii = 1:length(sim_sac_times)
    sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
        all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
    all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
end
all_sim_msacs = all_sim_sacs(all_trial_Fs(all_trialvec(all_sim_sacs)) == 0.2);
all_sim_gsacs = all_sim_sacs(all_trial_Fs(all_trialvec(all_sim_sacs)) == 3);

imback_gs_expts = find(~ismember(1:length(cur_block_set),sim_sac_expts) & expt_imback(cur_block_set));
grayback_gs_expts = find(~ismember(1:length(cur_block_set),sim_sac_expts) & ~expt_imback(cur_block_set));


%define which saccades to use
used_msacs = find(is_micro & ~is_blink);
used_gsacs = find(is_gsac & ~is_blink);

trial_start_inds = find(all_tsince_start < dt);
simsac_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),sim_sac_expts));

used_msacs_grayback = used_msacs(ismember(all_blockvec(interp_sac_start_inds(used_msacs)),grayback_gs_expts));
used_msacs_imback = used_msacs(ismember(all_blockvec(interp_sac_start_inds(used_msacs)),imback_gs_expts));
used_gsacs_grayback = used_gsacs(ismember(all_blockvec(interp_sac_start_inds(used_gsacs)),grayback_gs_expts));
used_gsacs_imback = used_gsacs(ismember(all_blockvec(interp_sac_start_inds(used_gsacs)),imback_gs_expts));
%microsacs excluding bursts
non_burst_msacs = used_msacs(~is_sacburst(used_msacs));
burst_msacs = used_msacs(is_sacburst(used_msacs));

msac_start_inds = interp_sac_start_inds(used_msacs);
gsac_start_inds = interp_sac_start_inds(used_gsacs);

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
    
    cur_uset = find(all_blockvec == ee);
    uinds = find(expt_lfp_t_axis >= all_t_axis(cur_uset(1)) & expt_lfp_t_axis <= all_t_axis(cur_uset(end)));
    full_lfps = cat(1,full_lfps,expt_lfps(uinds,:));
    full_lfp_taxis = cat(1,full_lfp_taxis,expt_lfp_t_axis(uinds));
    
    cur_toffset = trial_toffset(ee);
end
%%
lfp_trialvec = round(interp1(all_t_axis,all_trialvec,full_lfp_taxis));
lfp_trial_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_trial_start_times));
lfp_trial_start_inds(isnan(lfp_trial_start_inds)) = [];
lfp_trial_stop_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_trial_end_times));
lfp_trial_stop_inds(isnan(lfp_trial_stop_inds)) = [];
lfp_sac_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),sac_start_times));
lfp_simsac_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_t_axis(all_sim_sacs)));
lfp_sim_msac_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_t_axis(all_sim_msacs)));
lfp_sim_gsac_start_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_t_axis(all_sim_gsacs)));

trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
trial_start_inds(isnan(trial_start_inds)) = [];
trial_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));
trial_stop_inds(isnan(trial_stop_inds)) = [];
sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
simsac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_t_axis(all_sim_sacs)));
sim_msac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_t_axis(all_sim_msacs)));
sim_gsac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_t_axis(all_sim_gsacs)));

%%
inc_prev_trial = 0;
nboot = [];
backlag = round(0.3*Fsd);
forwardlag = round(0.6*Fsd);

[trial_onset_avg,lags] = get_event_trig_avg(full_lfps,lfp_trial_start_inds,backlag,forwardlag,nboot,lfp_trialvec,1);
[trial_offset_avg,lags] = get_event_trig_avg(full_lfps,lfp_trial_stop_inds,backlag,forwardlag,nboot,lfp_trialvec,1);

[gsac_avg,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(used_gsacs),backlag,forwardlag,nboot,lfp_trialvec,inc_prev_trial);
[gsac_grayback,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(used_gsacs_grayback),backlag,forwardlag,nboot,lfp_trialvec,inc_prev_trial);
[gsac_imback,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(used_gsacs_imback),backlag,forwardlag,nboot,lfp_trialvec,inc_prev_trial);

[msac_avg,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(used_msacs),backlag,forwardlag,nboot,lfp_trialvec,inc_prev_trial);
[msac_grayback,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(used_msacs_grayback),backlag,forwardlag,nboot,lfp_trialvec,inc_prev_trial);
[msac_imback,lags] = get_event_trig_avg(full_lfps,lfp_sac_start_inds(used_msacs_imback),backlag,forwardlag,nboot,lfp_trialvec,inc_prev_trial);

[simsac_avg,lags] = get_event_trig_avg(full_lfps,lfp_simsac_start_inds,backlag,forwardlag,nboot,lfp_trialvec,inc_prev_trial);
[sim_msac_avg,lags] = get_event_trig_avg(full_lfps,lfp_sim_msac_start_inds,backlag,forwardlag,nboot,lfp_trialvec,inc_prev_trial);
[sim_gsac_avg,lags] = get_event_trig_avg(full_lfps,lfp_sim_gsac_start_inds,backlag,forwardlag,nboot,lfp_trialvec,inc_prev_trial);

[gsac_rate,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(used_gsacs),backlag,forwardlag,nboot,all_trialvec,inc_prev_trial);
[gsac_grayback_rate,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(used_gsacs_grayback),backlag,forwardlag,nboot,all_trialvec,inc_prev_trial);
[gsac_imback_rate,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(used_gsacs_imback),backlag,forwardlag,nboot,all_trialvec,inc_prev_trial);

[msac_rate,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(used_msacs),backlag,forwardlag,nboot,all_trialvec,inc_prev_trial);
[msac_grayback_rate,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(used_msacs_grayback),backlag,forwardlag,nboot,all_trialvec,inc_prev_trial);
[msac_imback_rate,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(used_msacs_imback),backlag,forwardlag,nboot,all_trialvec,inc_prev_trial);

[simsac_rate,lags] = get_event_trig_avg(all_mua_rate_norm,simsac_start_inds,backlag,forwardlag,nboot,all_trialvec,inc_prev_trial);
[sim_msac_rate,lags] = get_event_trig_avg(all_mua_rate_norm,sim_msac_start_inds,backlag,forwardlag,nboot,all_trialvec,inc_prev_trial);
[sim_gsac_rate,lags] = get_event_trig_avg(all_mua_rate_norm,sim_gsac_start_inds,backlag,forwardlag,nboot,all_trialvec,inc_prev_trial);


%%
csd_params.Fs = Fsd; %sample freq
csd_params.BrainBound = 1; %first channel that is in the brain
csd_params.ChanSep = 0.05; %channel sep in mm
csd_params.diam = 2; %current disc diameter (in mm)

% [trial_onset_csd,lags] = get_event_trig_csd(full_lfps,trial_start_inds,backlag,forwardlag,csd_params,lfp_trialvec,1);
% [trial_offset_csd,lags] = get_event_trig_csd(full_lfps,trial_stop_inds,backlag,forwardlag,csd_params,lfp_trialvec,1);
[trial_onset_csd,lags] = get_event_trig_csd(full_lfps,trial_start_inds,backlag,forwardlag,csd_params);
[trial_offset_csd,lags] = get_event_trig_csd(full_lfps,trial_stop_inds,backlag,forwardlag,csd_params);

[gsac_csd,lags] = get_event_trig_csd(full_lfps,lfp_sac_start_inds(used_gsacs),backlag,forwardlag,csd_params,lfp_trialvec,inc_prev_trial);
[msac_csd,lags] = get_event_trig_csd(full_lfps,lfp_sac_start_inds(used_msacs),backlag,forwardlag,csd_params,lfp_trialvec,inc_prev_trial);
[simsac_csd,lags] = get_event_trig_csd(full_lfps,lfp_simsac_start_inds,backlag,forwardlag,csd_params,lfp_trialvec,inc_prev_trial);

%%
cd(save_dir);
save(save_name,'lags','Fsd','csd_params','*_csd','*_rate','*_avg');

%%
figure('Visible','off');
subplot(2,1,1)
imagesc(lags/Fsd,1:24,trial_onset_avg');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
xlim([-0.1 0.3]);
line([0 0],yl,'color','k');
xlabel('Time since trial onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
subplot(2,1,2)
imagesc(lags/Fsd,1:24,trial_offset_avg');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
xlim([-0.1 0.3]);
line([0 0],yl,'color','k');
xlabel('Time since trial offset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
fillPage(gcf,'papersize',[6 10]);
fname = 'trial_trig_lfps';
print(fname,'-dpdf','-painters');
close

figure('Visible','off');
subplot(3,1,1);
imagesc(lags/Fsd,1:24,gsac_avg');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
line([0 0],yl,'color','k');
xlim([-0.2 0.5])
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Guided saccades','fontsize',14);
subplot(3,1,2);
imagesc(lags/Fsd,1:24,msac_avg');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
xlim([-0.2 0.5])
line([0 0],yl,'color','k');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Micro saccades','fontsize',14);
subplot(3,1,3);
imagesc(lags/Fsd,1:24,simsac_avg');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar;
yl = ylim();
xlim([-0.2 0.5])
line([0 0],yl,'color','k');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Simulated saccades','fontsize',14);
fillPage(gcf,'papersize',[6 15]);
fname = 'sac_trig_lfps';
print(fname,'-dpdf','-painters');
close

figure('Visible','off');
subplot(3,1,1);
imagesc(lags/Fsd,1:24,gsac_avg');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
line([0 0],yl,'color','k');
xlim([-0.2 0.5])
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Guided saccades','fontsize',14);
subplot(3,1,2);
imagesc(lags/Fsd,1:24,msac_avg');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
xlim([-0.2 0.5])
line([0 0],yl,'color','k');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Micro saccades','fontsize',14);
subplot(3,1,3);
imagesc(lags/Fsd,1:24,simsac_avg');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar;
yl = ylim();
xlim([-0.2 0.5])
line([0 0],yl,'color','k');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Simulated saccades','fontsize',14);
fillPage(gcf,'papersize',[6 15]);
fname = 'sac_trig_lfps';
print(fname,'-dpdf','-painters');
close

figure('Visible','off');
subplot(2,1,1);
imagesc(lags/Fsd,1:24,gsac_grayback');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
line([0 0],yl,'color','k');
xlim([-0.2 0.4])
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Guided saccades','fontsize',14);
subplot(2,1,2);
imagesc(lags/Fsd,1:24,msac_grayback');
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
xlim([-0.2 0.4])
line([0 0],yl,'color','k');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Micro saccades','fontsize',14);
fillPage(gcf,'papersize',[6 12]);
fname = 'sac_trig_lfps_grayback';
print(fname,'-dpdf','-painters');
close

%%
figure('Visible','off');
subplot(2,1,1)
imagesc(lags/Fsd,1:24,trial_onset_csd);
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.75);
colorbar
yl = ylim();
line([0 0],yl,'color','k');
xlabel('Time since trial onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Trial onset','fontsize',14);
xlim([-0.1 0.3]);
subplot(2,1,2)
imagesc(lags/Fsd,1:24,trial_offset_csd);
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.75);
colorbar
yl = ylim();
line([0 0],yl,'color','k');
xlabel('Time since trial offset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
xlim([-0.1 0.3]);
title('Trial offset','fontsize',14);
fillPage(gcf,'papersize',[6 10]);
fname = 'trial_trig_csd';
print(fname,'-dpdf','-painters');
close

figure('Visible','off');
subplot(3,1,1);
imagesc(lags/Fsd,1:24,gsac_csd);
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
xlim([-0.2 0.5])
line([0 0],yl,'color','k');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Guided saccades','fontsize',14);
subplot(3,1,2);
imagesc(lags/Fsd,1:24,msac_csd);
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
xlim([-0.2 0.5])
line([0 0],yl,'color','k');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Micro saccades','fontsize',14);
subplot(3,1,3);
imagesc(lags/Fsd,1:24,simsac_csd);
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.85);
colorbar
yl = ylim();
xlim([-0.2 0.5])
line([0 0],yl,'color','k');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
set(gca,'fontsize',10,'fontname','arial');
title('Simulated saccades','fontsize',14);
fillPage(gcf,'papersize',[6 15]);
fname = 'sac_trig_csd';
print(fname,'-dpdf','-painters');
close