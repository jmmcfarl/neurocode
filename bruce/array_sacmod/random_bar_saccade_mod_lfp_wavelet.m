clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G086';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/G081/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

recompute_prepdata = 1;
prep_data_name = [anal_dir '/sac_prepdata.mat'];
recompute_sacdata = 1;
sac_data_name = [anal_dir '/sac_sacdata.mat'];

%% PARSE TRIAL DATA STRUCTURES

stim_fs = 100; %in Hz

dsf = 2; %lfps originally sampled at 400Hz
use_lfps = 1:16:96;
Fs_app =400;
nwfreqs = 35;
min_freq = 2; max_freq = 100;
min_scale = 1/max_freq*Fs_app;
max_scale = 1/min_freq*Fs_app;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fs_app);

%%
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
%% COMPUTE TRIAL DATA
cd(data_dir);

%%
fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_tsince_start = [];
all_exptvec = [];
all_trialvec = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_V = [];
all_ampgrams = [];
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
    
    trial_start_times = trial_start_times(use_trials);
    trial_end_times = trial_end_times(use_trials);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times');
    n_trials = length(use_trials);
    
    %load lfps
    lfp_fname = sprintf('Expt%d_LFP.mat',cur_expt);
    load(lfp_fname);
    Fs = lfp_params.Fsd;
    cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
    cur_ampgrams = nan(size(cur_lfps,1),length(wfreqs),length(use_lfps));
    for ll = 1:length(use_lfps)
        temp = cwt(cur_lfps(:,ll),scales,'cmor1-1');
        cur_ampgrams(:,:,ll) = abs(temp)';
    end
    
    if dsf > 1
        [filt_b,filt_a] = butter(4,0.8/dsf,'low');
        for ii = 1:length(use_lfps)
            cur_lfps(:,ii) = filtfilt(filt_b,filt_a,cur_lfps(:,ii));
        end
        cur_lfps = downsample(cur_lfps,dsf);
        cur_lfp_t = downsample(lfp_t_ax,dsf);
        cur_ampgrams = downsample(cur_ampgrams,dsf);
    end
    all_ampgrams = cat(1,all_ampgrams,cur_ampgrams);
    all_V = cat(1,all_V,cur_lfps);
    all_t_axis = cat(1,all_t_axis,cur_lfp_t');
    
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

%% NORMALIZE BY POWER LAW FIT

norm_used_inds = find(all_tsince_start > 0.5);
avg_pow_spectra = squeeze(mean(all_ampgrams(norm_used_inds,:,:)));
std_pow_spectra = squeeze(std(all_ampgrams(norm_used_inds,:,:)));

% all_ampgrams_norm = all_ampgrams;
% 
% un_wfreqs = linspace(wfreqs(1),wfreqs(end),100)';
% avg_pow_int = interp1(wfreqs,avg_pow_spectra,un_wfreqs);
% % powfun = @(a,x)(a(1)*x.^a(2));
% powfun = @(a,x)(a(1)*exp(-x*a(2)));
% for ll = 1:length(use_lfps)
% %     BETA0 = [(avg_pow_int(end,ll))/(un_wfreqs(end)^-2) -2];
%     BETA0 = [(avg_pow_int(end,ll))/(exp(-un_wfreqs(end)*2)) 2];
%     BETA = nlinfit(un_wfreqs,avg_pow_int(:,ll),powfun,BETA0);
% %     powfit = un_wfreqs.^BETA(2)*BETA(1);
%     powfit = BETA(1)*exp(-un_wfreqs*BETA(2));
%     interp_powfit = interp1(un_wfreqs,powfit,wfreqs);
%     all_ampgrams_norm(:,:,ll) = bsxfun(@rdivide,all_ampgrams(:,:,ll),interp_powfit);
% end
% avg_pow_spectra_norm = squeeze(mean(all_ampgrams_norm(norm_used_inds,:,:)));

all_ampgrams_norm = bsxfun(@minus,all_ampgrams,reshape(avg_pow_spectra,[1 length(wfreqs) length(use_lfps)]));
all_ampgrams_norm = bsxfun(@rdivide,all_ampgrams_norm,reshape(std_pow_spectra,[1 length(wfreqs) length(use_lfps)]));

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)

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

%% PICK OUT SACCADES FOR ANALYSIS
sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > 1/Fsd | isnan(all_trialvec(interp_sac_start_inds))');
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
is_sacburst = (isis(1:end-1) < 0.15 | isis(2:end) < 0.15);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps < 1;

sac_post_Lx = [saccades(:).post_Lx];
sac_post_Ly = [saccades(:).post_Ly];
sac_pre_Lx = [saccades(:).pre_Lx];
sac_pre_Ly = [saccades(:).pre_Ly];
sac_deltaX = sac_post_Lx - sac_pre_Lx;
sac_deltaY = sac_post_Ly - sac_pre_Ly;

sac_dirs = [0 90];
sac_start_expts = all_exptvec(interp_sac_start_inds);
delta_sacpar = nan(size(saccades));
for bb = 1:2
    cur_bar_expts = find(expt_bar_ori(cur_expt_set) == sac_dirs(bb));
    cur_sac_set = find(ismember(sac_start_expts,cur_bar_expts));
    
    cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));
    delta_sacpar(cur_sac_set) = cur_delta_sac_par;
end
is_gsac = delta_sacpar' > 1;

sim_sac_expts = find(~expt_has_ds(cur_expt_set));

sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
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

% all_sim_sacs(all_tsince_start(all_sim_sacs) <= 1) = [];
% used_msacs(all_tsince_start(interp_sac_start_inds(used_msacs)) <= 1) = [];
% used_gsacs(all_tsince_start(interp_sac_start_inds(used_gsacs)) <= 1) = [];

imback_gs_expts = find(expt_has_ds(cur_expt_set) & expt_imback(cur_expt_set)');
grayback_gs_expts = find(expt_has_ds(cur_expt_set) & ~expt_imback(cur_expt_set)' & ~expt_binoc(cur_expt_set));
hori_expts = find(expt_bar_ori(cur_expt_set) == 0);
ver_expts = find(expt_bar_ori(cur_expt_set) == 90);

trial_start_inds = find(all_tsince_start < 1/Fsd);
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
nboot = [];

%trial averages
[simsac_trial_avgs,trial_lags] = get_event_trig_avg(all_ampgrams_norm,simsac_trial_starts,0,round(4*Fsd));
[simsac_trial_avgs_raw,trial_lags] = get_event_trig_avg(all_V,simsac_trial_starts,0,round(4*Fsd));
% [grayback_trial_avgs,trial_lags] = get_event_trig_avg(all_ampgrams_norm,grayback_trial_starts,0,round(4*Fsd));
% [imback_trial_avgs,trial_lags] = get_event_trig_avg(all_ampgrams_norm,imback_trial_starts,0,round(4*Fsd));

backlag = round(0.4*Fsd);
forwardlag = round(0.6*Fsd);

%general averages
[msac_avgs,lags] = get_event_trig_avg(all_ampgrams_norm,interp_sac_start_inds(used_msacs),backlag,forwardlag,nboot,all_trialvec);
[gsac_avgs,lags] = get_event_trig_avg(all_ampgrams_norm,interp_sac_start_inds(used_gsacs),backlag,forwardlag,nboot,all_trialvec);
[simsac_avgs,lags] = get_event_trig_avg(all_ampgrams_norm,all_sim_sacs,backlag,forwardlag,nboot,all_trialvec);

% %stop- and peak-triggered averages
% [msac_stop_avgs,lags] = get_event_trig_avg(all_V,interp_sac_stop_inds(used_msacs),backlag,forwardlag,nboot,all_trialvec);
% [gsac_stop_avgs,lags] = get_event_trig_avg(all_V,interp_sac_stop_inds(used_gsacs),backlag,forwardlag,nboot,all_trialvec);
% [msac_peak_avgs,lags] = get_event_trig_avg(all_V,interp_sac_peak_inds(used_msacs),backlag,forwardlag,nboot,all_trialvec);
% [gsac_peak_avgs,lags] = get_event_trig_avg(all_V,interp_sac_peak_inds(used_gsacs),backlag,forwardlag,nboot,all_trialvec);
% 
% %background dependent
% [msac_gray_avgs,lags] = get_event_trig_avg(all_V,interp_sac_start_inds(used_msacs_grayback),backlag,forwardlag,nboot,all_trialvec);
% [gsac_gray_avgs,lags] = get_event_trig_avg(all_V,interp_sac_start_inds(used_gsacs_grayback),backlag,forwardlag,nboot,all_trialvec);
% [msac_im_avgs,lags] = get_event_trig_avg(all_V,interp_sac_start_inds(used_msacs_imback),backlag,forwardlag,nboot,all_trialvec);
% [gsac_im_avgs,lags] = get_event_trig_avg(all_V,interp_sac_start_inds(used_gsacs_imback),backlag,forwardlag,nboot,all_trialvec);
% 
% %saccade direction dependent
% [msac_ver_avgs,lags] = get_event_trig_avg(all_V,interp_sac_start_inds(used_msacs_ver),backlag,forwardlag,nboot,all_trialvec);
% [gsac_ver_avgs,lags] = get_event_trig_avg(all_V,interp_sac_start_inds(used_gsacs_ver),backlag,forwardlag,nboot,all_trialvec);
% [simsac_ver_avgs,lags] = get_event_trig_avg(all_V,simsacs_ver,backlag,forwardlag,nboot,all_trialvec);
% [msac_hor_avgs,lags] = get_event_trig_avg(all_V,interp_sac_start_inds(used_msacs_hori),backlag,forwardlag,nboot,all_trialvec);
% [gsac_hor_avgs,lags] = get_event_trig_avg(all_V,interp_sac_start_inds(used_gsacs_hori),backlag,forwardlag,nboot,all_trialvec);
% [simsac_hor_avgs,lags] = get_event_trig_avg(all_V,simsacs_hori,backlag,forwardlag,nboot,all_trialvec);

%%
% sname = [anal_dir '/sacmod_LFP_data.mat'];
% save(sname,'lags','Fsd','*_avgs');

%% PRINT ENSEMBLE AVERAGES
to_print = 1;

shadedErrorBar(trial_lags/Fsd,mean(simsac_trial_avgs_raw,2),std(simsac_trial_avgs_raw,[],2)/sqrt(96),{'color','r'});
ylim([-4 4]*1e-5)
xl = xlim();yl = ylim();
line(xl,[0 0],'color','k','linestyle','--');
for ii = 1:length(sim_sac_times)
    line(sim_sac_times([ii ii]),yl,'color','k','linestyle','--')
end


fname0 = [anal_dir '/simsac_fulltrial_all_LFPs'];
if to_print == 1
    figure('visible','off');
else
    figure
end
plot(trial_lags/Fsd,simsac_trial_avgs,'k')
hold on
shadedErrorBar(trial_lags/Fsd,mean(simsac_trial_avgs,2),std(simsac_trial_avgs,[],2)/sqrt(96),{'color','r'});
ylim([-4 4]*1e-5)
xl = xlim();yl = ylim();
line(xl,[0 0],'color','k','linestyle','--');
for ii = 1:length(sim_sac_times)
    line(sim_sac_times([ii ii]),yl,'color','k','linestyle','--')
end
if to_print == 1
    print(fname0,'-dpsc');
    close all
end

% hold on
% shadedErrorBar(trial_lags/Fsd,mean(grayback_trial_avgs,2),std(grayback_trial_avgs,[],2)/sqrt(96),{'color','r'});
% shadedErrorBar(trial_lags/Fsd,mean(imback_trial_avgs,2),std(imback_trial_avgs,[],2)/sqrt(96),{'color','b'});


fname1 = [anal_dir '/sacmod_all_LFPs'];
if to_print == 1
    figure('visible','off');
else
    figure
end
shadedErrorBar(lags/Fsd,mean(simsac_avgs,2),std(simsac_avgs,[],2)/sqrt(96));
hold on
shadedErrorBar(lags/Fsd,mean(gsac_avgs,2),std(gsac_avgs,[],2)/sqrt(96),{'color','r'});
shadedErrorBar(lags/Fsd,mean(msac_avgs,2),std(msac_avgs,[],2)/sqrt(96),{'color','b'});
xlim([-backlag/Fsd forwardlag/Fsd])
xl = xlim();yl = ylim();
line(xl,[0 0],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
if to_print == 1
    print(fname1,'-dpsc');
    close all
end

fname2 = [anal_dir '/sacmod_all_LFPs_dircompare'];
if to_print == 1
    figure('visible','off');
else
    figure
end
subplot(3,1,1)
shadedErrorBar(lags/Fsd,mean(gsac_ver_avgs,2),std(gsac_ver_avgs,[],2)/sqrt(96));
hold on
shadedErrorBar(lags/Fsd,mean(gsac_hor_avgs,2),std(gsac_hor_avgs,[],2)/sqrt(96),{'color','r'});
xlim([-backlag*dt forwardlag*dt])
xl = xlim();yl = ylim();
line(xl,[0 0],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
subplot(3,1,2)
shadedErrorBar(lags/Fsd,mean(msac_ver_avgs,2),std(msac_ver_avgs,[],2)/sqrt(96));
hold on
shadedErrorBar(lags/Fsd,mean(msac_hor_avgs,2),std(msac_hor_avgs,[],2)/sqrt(96),{'color','r'});
xlim([-backlag/Fsd forwardlag/Fsd])
xl = xlim();yl = ylim();
line(xl,[0 0],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
subplot(3,1,3)
shadedErrorBar(lags/Fsd,mean(simsac_hor_avgs,2),std(simsac_hor_avgs,[],2)/sqrt(96));
hold on
shadedErrorBar(lags/Fsd,mean(simsac_ver_avgs,2),std(simsac_ver_avgs,[],2)/sqrt(96),{'color','r'});
xlim([-backlag/Fsd forwardlag/Fsd])
xl = xlim();yl = ylim();
line(xl,[0 0],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
if to_print == 1
    fillPage(gcf,'papersize',[8 15]);
    print(fname2,'-dpsc');
    close all
end

fname3 = [anal_dir '/sacmod_all_LFPs_backcomp'];
if to_print == 1
    figure('visible','off');
else
    figure
end
subplot(2,1,1)
shadedErrorBar(lags/Fsd,mean(gsac_gray_avgs,2),std(gsac_gray_avgs,[],2)/sqrt(96));
hold on
shadedErrorBar(lags/Fsd,mean(gsac_im_avgs,2),std(gsac_im_avgs,[],2)/sqrt(96),{'color','r'});
xlim([-backlag*dt forwardlag*dt])
xl = xlim();yl = ylim();
line(xl,[0 0],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
subplot(2,1,2)
shadedErrorBar(lags/Fsd,mean(msac_gray_avgs,2),std(msac_gray_avgs,[],2)/sqrt(96));
hold on
shadedErrorBar(lags/Fsd,mean(msac_im_avgs,2),std(msac_im_avgs,[],2)/sqrt(96),{'color','r'});
xlim([-backlag*dt forwardlag*dt])
xl = xlim();yl = ylim();
line(xl,[0 0],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',16);
ylabel('Relative firing rate','fontsize',16);
title('Overall average');
if to_print == 1
    fillPage(gcf,'papersize',[8 12]);
    print(fname3,'-dpsc');
    close all
end
