% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_num = 92;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

%%
dsf = 1; %lfps originally sampled at 400Hz
use_lfps = 1:2:96;
lcf = 0.5;

cur_block_set = [28:43];
Fsd = 400;
start_buffer = round(Fsd*0.5);
min_trial_dur = 0.5;
%%
fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_tsince_start = [];
all_ttill_end = [];
all_blockvec = [];
all_trialvec = [];
all_trial_durs = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_V = [];
all_trial_exptvec = [];
all_trial_nph = [];
all_trial_jv = [];
all_trial_st = [];
all_trial_sl = [];
all_trial_opt = [];
all_trial_tf = [];
for ee = 1:length(cur_block_set);
    fprintf('Expt %d of %d\n',ee,length(cur_block_set));
    cur_expt = cur_block_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    trial_jv = [Expts{cur_expt}.Trials(:).jv];
    trial_opt = [Expts{cur_expt}.Trials(:).optionb];
    trial_st = [Expts{cur_expt}.Trials(:).st];
    trial_sl = [Expts{cur_expt}.Trials(:).sl];
    trial_nph = [Expts{cur_expt}.Trials(:).nph];
    trial_tf = [Expts{cur_expt}.Trials(:).tf];
    
    [un_ids,id_inds] = unique(trial_ids);
    use_trials = id_inds(trial_durs(id_inds) >= min_trial_dur);
    
    trial_start_times = trial_start_times(use_trials);
    trial_end_times = trial_end_times(use_trials);
    trial_durs = trial_durs(use_trials);
    trial_jv = trial_jv(use_trials);
    trial_tf = trial_tf(use_trials);
    trial_opt = trial_opt(use_trials);
    trial_nph = trial_nph(use_trials);
    trial_st = trial_st(use_trials);
    trial_sl = trial_sl(use_trials);
    
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times');
    all_trial_durs = cat(1,all_trial_durs,trial_durs');
    all_trial_nph = cat(1,all_trial_nph,trial_nph');
    all_trial_tf = cat(1,all_trial_tf,trial_tf');
    all_trial_jv = cat(1,all_trial_jv,trial_jv');
    all_trial_st = cat(1,all_trial_st,trial_st');
    all_trial_sl = cat(1,all_trial_sl,trial_sl');
    all_trial_opt = cat(1,all_trial_opt,trial_opt');
    
    all_trial_exptvec = cat(1,all_trial_exptvec,ee*ones(length(use_trials),1));
    
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
    expt_ttill_end = nan(length(cur_lfp_t),1);
    expt_exptvec = nan(length(cur_lfp_t),1);
    expt_trialvec = nan(length(cur_lfp_t),1);
    for tt = 1:n_trials
        cur_samples = find(cur_lfp_t >= trial_start_times(tt) & cur_lfp_t <= trial_end_times(tt));
        
        expt_tsince_start(cur_samples) = cur_lfp_t(cur_samples) - trial_start_times(tt);
        expt_ttill_end(cur_samples) = trial_end_times(tt)-cur_lfp_t(cur_samples);
        expt_exptvec(cur_samples) = ee;
        expt_trialvec(cur_samples) = tt + trial_cnt;
    end
    trial_cnt = trial_cnt + n_trials;
    
    all_blockvec = cat(1,all_blockvec,expt_exptvec);
    all_tsince_start = cat(1,all_tsince_start,expt_tsince_start);
    all_ttill_end = cat(1,all_ttill_end,expt_ttill_end);
    all_trialvec = cat(1,all_trialvec,expt_trialvec);
end
Fsd = Fs/dsf;

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)

emfile = ['jbe' Expt_name '.em.mat'];
load(emfile);

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
all_eye_blockvec = [];
eye_smooth = 3;
for ee = 1:length(cur_block_set);
    fprintf('Loading ET data for expt %d, block %d of %d\n',Expt_num,ee,length(cur_block_set));
    cur_set = find(all_blockvec==ee);
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

sac_dbuff = round(0.05/eye_dt);
sac_pre_pos = nan(length(saccade_inds),4);
sac_post_pos = nan(length(saccade_inds),4);
for ii = 1:length(saccade_inds)
    pre_inds = (sac_start_inds(ii) - sac_dbuff):sac_start_inds(ii);
    pre_inds(pre_inds < 1) = 1;
    sac_pre_pos(ii,:) = median(all_eye_vals(pre_inds,:),1);
    post_inds = sac_stop_inds(ii):(sac_stop_inds(ii) + sac_dbuff);
    post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
    sac_post_pos(ii,:) = median(all_eye_vals(post_inds,:),1);
end

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
dt = 1/Fsd;
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

poss_use_sacs = find(~is_sacburst' & ~is_blink);

small_sacs = poss_use_sacs(sac_amps(poss_use_sacs) < 0.2);
bigger_sacs = poss_use_sacs(sac_amps(poss_use_sacs) > 1);
medium_sacs = poss_use_sacs(sac_amps(poss_use_sacs) < 0.7 & sac_amps(poss_use_sacs) > 0.3);

sac_dir = atan2(sac_delta_pos(:,2),sac_delta_pos(:,1));
hor_sacs = poss_use_sacs(sac_dir(poss_use_sacs) > -pi/4 & sac_dir(poss_use_sacs) < pi/4);
ver_sacs = poss_use_sacs(sac_dir(poss_use_sacs) < -3*pi/4 | sac_dir(poss_use_sacs) > 3*pi/4);

%% CLASSIFY TRIALS
used_inds = find(all_tsince_start > start_buffer/Fsd);
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));


all_trial_jv(all_trial_st == 3) = 0;
conditions = [all_trial_opt(:) all_trial_st(:) all_trial_sl(:) all_trial_tf(:) all_trial_nph(:) all_trial_jv(:)];
[C,IA,IC] = unique(conditions,'rows');
cond_spec{1} = 'opt';
cond_spec{2} = 'st';
cond_spec{3} = 'sl';
cond_spec{4} = 'tf';
cond_spec{5} = 'nph';
cond_spec{6} = 'jv';
fprintf('%d unique conditions\n',length(C));

%% TRIG AVG WAVELET SPECGRAMS

%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 40;
min_freq = 2; max_freq = 125;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

forlag = round(0.6*Fsd);
backlag = round(0.3*Fsd);
lags = -backlag:forlag;

inc_prev_trial = 1;
nboot = 0;
n_conds = length(C);

% % use_conds = [2 3 12 13 14 15 16];
% % used_trials = find(ismember(IC,use_conds) & all_trial_durs > start_buffer/Fsd);
% % used_inds = find(ismember(all_trialvec,used_trials));
% % used_sacs = find(ismember(interp_sac_start_inds,used_inds));
% % used_sacs(is_blink(used_sacs)) = [];
% % use_sac_inds = find(ismember(used_inds,interp_sac_start_inds(used_sacs)));
% % fprintf('Using %d sacs\n',length(use_sac_inds));
% cond_trig_spec = nan(length(use_lfps),n_conds,length(lags),length(wfreqs));
% cond_trig_lfp = nan(length(use_lfps),n_conds,length(lags));
% for cc = 1:n_conds
%     
%     used_trials = find(IC == cc & all_trial_durs > start_buffer/Fsd);
%     used_inds = find(ismember(all_trialvec,used_trials));
%     used_sacs = find(ismember(interp_sac_start_inds,used_inds));
%     used_sacs(is_blink(used_sacs)) = [];
%     use_sac_inds = find(ismember(used_inds,interp_sac_start_inds(used_sacs)));
%     fprintf('Using %d sacs\n',length(use_sac_inds));
%     
%     for ll = 1:length(use_lfps);
%         fprintf('LFP %d of %d\n',ll,length(use_lfps));
%         temp = cwt(all_V(used_inds,ll),scales,'cmor1-1');
%         cur_ampgram = abs(temp)';
%         
%         cond_trig_spec(ll,cc,:,:) = get_event_trig_avg(cur_ampgram,use_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
%         cond_trig_lfp(ll,cc,:) = get_event_trig_avg(all_V(used_inds,ll),use_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
%     end
% end
% %%
% cd ~/Analysis/bruce/ 
% load gamma_powfits
% znorm_cond_spec = bsxfun(@minus,cond_trig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,1,length(wfreqs)]));
% znorm_cond_spec = bsxfun(@rdivide,znorm_cond_spec,reshape(all_std_spectrum,[length(use_lfps),1,1,length(wfreqs)]));
% 
% avg_znorm_cond_spec = squeeze(mean(znorm_cond_spec));
% 
% for cc = 1:n_conds
%     fprintf('Condition %d of %d\n',cc,n_conds);
%     disp(C(cc,:))
%     subplot(3,1,[1 2]);
%     pcolor(lags/Fsd,wfreqs,squeeze(avg_znorm_cond_spec(cc,:,:))');shading flat
%     caxis([-0.75 0.75]);
%     xlim([-backlag forlag]/Fsd)
%     subplot(3,1,3);
%     shadedErrorBar(lags/Fsd,squeeze(mean(cond_trig_lfp(:,cc,:))),squeeze(std(cond_trig_lfp(:,cc,:))));
%     xlim([-backlag forlag]/Fsd)
%     pause
%     clf
% end
% 
% %%
% use_conds = [2 3 12 13 14 15 16];
% % use_conds = [2 3 12 13 15 16];
% % use_conds = [14];
% used_trials = find(ismember(IC,use_conds) & all_trial_durs > start_buffer/Fsd);
% used_inds = find(ismember(all_trialvec,used_trials));
% % used_sacs = find(ismember(interp_sac_start_inds,used_inds));
% % used_sacs(is_blink(used_sacs)) = [];
% % use_sac_inds = find(ismember(used_inds,interp_sac_start_inds(used_sacs)));
% used_sacs = find(ismember(interp_sac_stop_inds,used_inds));
% used_sacs(is_blink(used_sacs)) = [];
% use_sac_inds = find(ismember(used_inds,interp_sac_stop_inds(used_sacs)));
% 
% big_sac_inds = use_sac_inds(ismember(used_sacs,bigger_sacs));
% med_sac_inds = use_sac_inds(ismember(used_sacs,medium_sacs));
% small_sac_inds = use_sac_inds(ismember(used_sacs,small_sacs));
% hor_sac_inds = use_sac_inds(ismember(used_sacs,hor_sacs));
% ver_sac_inds = use_sac_inds(ismember(used_sacs,ver_sacs));
% fprintf('Using %d big sacs  %d med sacs  %d small sacs  %d hor sacs  %d ver sacs\n',...
%     length(big_sac_inds),length(med_sac_inds),length(small_sac_inds),length(hor_sac_inds),length(ver_sac_inds));
% 
% small_trig_spec = nan(length(use_lfps),length(lags),length(wfreqs));
% small_trig_lfp = nan(length(use_lfps),length(lags));
% med_trig_spec = nan(length(use_lfps),length(lags),length(wfreqs));
% med_trig_lfp = nan(length(use_lfps),length(lags));
% big_trig_spec = nan(length(use_lfps),length(lags),length(wfreqs));
% big_trig_lfp = nan(length(use_lfps),length(lags));
% % hor_trig_spec = nan(length(use_lfps),length(lags),length(wfreqs));
% % hor_trig_lfp = nan(length(use_lfps),length(lags));
% % ver_trig_spec = nan(length(use_lfps),length(lags),length(wfreqs));
% % ver_trig_lfp = nan(length(use_lfps),length(lags));
% 
% for ll = 1:length(use_lfps);
%     fprintf('LFP %d of %d\n',ll,length(use_lfps));
%     temp = cwt(all_V(used_inds,ll),scales,'cmor1-1');
%     cur_ampgram = abs(temp)';
%     
%     small_trig_spec(ll,:,:) = get_event_trig_avg(cur_ampgram,small_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
%     small_trig_lfp(ll,:) = get_event_trig_avg(all_V(used_inds,ll),small_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
%     med_trig_spec(ll,:,:) = get_event_trig_avg(cur_ampgram,med_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
%     med_trig_lfp(ll,:) = get_event_trig_avg(all_V(used_inds,ll),med_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
%     big_trig_spec(ll,:,:) = get_event_trig_avg(cur_ampgram,big_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
%     big_trig_lfp(ll,:) = get_event_trig_avg(all_V(used_inds,ll),big_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
% %     hor_trig_spec(ll,:,:) = get_event_trig_avg(cur_ampgram,hor_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
% %     hor_trig_lfp(ll,:) = get_event_trig_avg(all_V(used_inds,ll),hor_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
% %     ver_trig_spec(ll,:,:) = get_event_trig_avg(cur_ampgram,ver_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
% %     ver_trig_lfp(ll,:) = get_event_trig_avg(all_V(used_inds,ll),ver_sac_inds,backlag,forlag,nboot,all_trialvec,inc_prev_trial);
% end
% 
% znorm_small_cond_spec = bsxfun(@minus,small_trig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,length(wfreqs)]));
% znorm_small_cond_spec = bsxfun(@rdivide,znorm_small_cond_spec,reshape(all_std_spectrum,[length(use_lfps),1,length(wfreqs)]));
% znorm_med_cond_spec = bsxfun(@minus,med_trig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,length(wfreqs)]));
% znorm_med_cond_spec = bsxfun(@rdivide,znorm_med_cond_spec,reshape(all_std_spectrum,[length(use_lfps),1,length(wfreqs)]));
% znorm_big_cond_spec = bsxfun(@minus,big_trig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,length(wfreqs)]));
% znorm_big_cond_spec = bsxfun(@rdivide,znorm_big_cond_spec,reshape(all_std_spectrum,[length(use_lfps),1,length(wfreqs)]));
% % znorm_hor_cond_spec = bsxfun(@minus,hor_trig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,length(wfreqs)]));
% % znorm_hor_cond_spec = bsxfun(@rdivide,znorm_hor_cond_spec,reshape(all_std_spectrum,[length(use_lfps),1,length(wfreqs)]));
% % znorm_ver_cond_spec = bsxfun(@minus,ver_trig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,length(wfreqs)]));
% % znorm_ver_cond_spec = bsxfun(@rdivide,znorm_ver_cond_spec,reshape(all_std_spectrum,[length(use_lfps),1,length(wfreqs)]));
% 
% avg_znorm_small_spec = squeeze(mean(znorm_small_cond_spec));
% avg_znorm_med_spec = squeeze(mean(znorm_med_cond_spec));
% avg_znorm_big_spec = squeeze(mean(znorm_big_cond_spec));
% % avg_znorm_hor_spec = squeeze(mean(znorm_hor_cond_spec));
% % avg_znorm_ver_spec = squeeze(mean(znorm_ver_cond_spec));
% 
% ca = [-0.25 0.75];
% xl = [-0.2 0.4];
% subplot(3,2,1)
% pcolor(lags/Fsd,wfreqs,avg_znorm_small_spec');shading flat
% caxis(ca);
% xlim(xl)
% subplot(3,2,2)
% shadedErrorBar(lags/Fsd,mean(small_trig_lfp),std(small_trig_lfp));
% xlim(xl)
% subplot(3,2,3)
% pcolor(lags/Fsd,wfreqs,avg_znorm_med_spec');shading flat
% caxis(ca);
% xlim(xl)
% subplot(3,2,4)
% shadedErrorBar(lags/Fsd,mean(med_trig_lfp),std(small_trig_lfp));
% xlim(xl)
% subplot(3,2,5)
% pcolor(lags/Fsd,wfreqs,avg_znorm_big_spec');shading flat
% caxis(ca);
% xlim(xl)
% subplot(3,2,6)
% shadedErrorBar(lags/Fsd,mean(big_trig_lfp),std(small_trig_lfp));
% xlim(xl)
% 
% %%
% for ll = 1:length(use_lfps)
%     ca = [-0.25 0.75];
%     xl = [-0.2 0.4];
%     subplot(3,2,1)
%     pcolor(lags/Fsd,wfreqs,squeeze(znorm_small_cond_spec(ll,:,:))');shading flat
%     caxis(ca);
%     xlim(xl)
%     subplot(3,2,2)
%     plot(lags/Fsd,small_trig_lfp(ll,:));
%     xlim(xl)
%     subplot(3,2,3)
%     pcolor(lags/Fsd,wfreqs,squeeze(znorm_med_cond_spec(ll,:,:))');shading flat
%     caxis(ca);
%     xlim(xl)
%     subplot(3,2,4)
%     plot(lags/Fsd,med_trig_lfp(ll,:));
%     xlim(xl)
%     subplot(3,2,5)
%     pcolor(lags/Fsd,wfreqs,squeeze(znorm_big_cond_spec(ll,:,:))');shading flat
%     caxis(ca);
%     xlim(xl)
%     subplot(3,2,6)
%     plot(lags/Fsd,big_trig_lfp(ll,:));
%     xlim(xl)
%     pause
%     clf
% end
% 
% % subplot(2,1,1)
% % pcolor(lags/Fsd,wfreqs,avg_znorm_hor_spec');shading flat
% % caxis([-0.25 0.7]);
% % xlim([-0.1 0.3])
% % subplot(2,1,2)
% % pcolor(lags/Fsd,wfreqs,avg_znorm_ver_spec');shading flat
% % caxis([-0.25 0.7]);
% % xlim([-0.1 0.3])

%%
cd ~/Analysis/bruce/ 
load gamma_powfits

tbacklag = 0;
tforlag = round(Fsd*4);

% use_conds = [2 3 12 13 14 15 16];
% use_conds = [2 3 12 13 15 16];
use_conds = [14];
used_trials = find(ismember(IC,use_conds) & all_trial_durs > start_buffer/Fsd);
used_inds = find(ismember(all_trialvec,used_trials));
usable_sacs = find(ismember(interp_sac_stop_inds,used_inds) & ismember(interp_sac_start_inds,used_inds));
usable_sacs(is_blink(usable_sacs) | is_sacburst(usable_sacs)') = [];
usable_sacs(all_ttill_end(interp_sac_stop_inds(usable_sacs)) < 0.2) = [];
usable_sac_stop_inds = find(ismember(used_inds,interp_sac_stop_inds(usable_sacs)));
usable_sac_start_inds = find(ismember(used_inds,interp_sac_start_inds(usable_sacs)));

used_trial_start_inds = find((ismember(used_inds,all_trial_start_inds(used_trials))));

n_prctile_bins = 2;
sac_amp_prctiles = prctile(sac_amps(usable_sacs),[linspace(0,100,n_prctile_bins)]);

sacamp_btrig_spec = nan(length(use_lfps),n_prctile_bins-1,length(lags),length(wfreqs));
sacamp_btrig_lfp = nan(length(use_lfps),n_prctile_bins-1,length(lags));
sacamp_etrig_spec = nan(length(use_lfps),n_prctile_bins-1,length(lags),length(wfreqs));
sacamp_etrig_lfp = nan(length(use_lfps),n_prctile_bins-1,length(lags));
for ll = 1:length(use_lfps);
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
%     temp = cwt(all_V(used_inds,ll),scales,'cmor1-1');
%     cur_ampgram = abs(temp)';
    
    for aa = 1:n_prctile_bins-1
        cur_sacs = find(sac_amps(usable_sacs) >= sac_amp_prctiles(aa) & sac_amps(usable_sacs) < sac_amp_prctiles(aa+1));
        fprintf('bin %d, %d sacs\n',aa,length(cur_sacs));
        
%         sacamp_btrig_spec(ll,aa,:,:) = get_event_trig_avg(cur_ampgram,usable_sac_start_inds(cur_sacs),backlag,forlag,nboot,all_trialvec,inc_prev_trial);
        sacamp_btrig_lfp(ll,aa,:) = get_event_trig_avg(all_V(used_inds,ll),usable_sac_start_inds(cur_sacs),backlag,forlag,nboot,all_trialvec,inc_prev_trial);
%         sacamp_etrig_spec(ll,aa,:,:) = get_event_trig_avg(cur_ampgram,usable_sac_stop_inds(cur_sacs),backlag,forlag,nboot,all_trialvec,inc_prev_trial);
        sacamp_etrig_lfp(ll,aa,:) = get_event_trig_avg(all_V(used_inds,ll),usable_sac_stop_inds(cur_sacs),backlag,forlag,nboot,all_trialvec,inc_prev_trial);
    end
    
%     [cond_trig_spec(ll,:,:),tlags] = get_event_trig_avg(cur_ampgram,used_trial_start_inds,tbacklag,tforlag,nboot,all_trialvec,inc_prev_trial);
    [cond_trig_lfp(ll,:),tlags] = get_event_trig_avg(all_V(used_inds,ll),used_trial_start_inds,tbacklag,tforlag,nboot,all_trialvec,inc_prev_trial);

end

znorm_sacamp_btrig_spec = bsxfun(@minus,sacamp_btrig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,1,length(wfreqs)]));
znorm_sacamp_btrig_spec = bsxfun(@rdivide,znorm_sacamp_btrig_spec,reshape(all_std_spectrum,[length(use_lfps),1,1,length(wfreqs)]));
znorm_sacamp_etrig_spec = bsxfun(@minus,sacamp_etrig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,1,length(wfreqs)]));
znorm_sacamp_etrig_spec = bsxfun(@rdivide,znorm_sacamp_etrig_spec,reshape(all_std_spectrum,[length(use_lfps),1,1,length(wfreqs)]));

avg_znorm_sacamp_btrig_spec = squeeze(nanmean(znorm_sacamp_btrig_spec));
avg_znorm_sacamp_etrig_spec = squeeze(nanmean(znorm_sacamp_etrig_spec));

avg_sacamp_btrig_lfp = squeeze(nanmean(sacamp_btrig_lfp));
avg_sacamp_btrig_std = squeeze(nanstd(sacamp_btrig_lfp));
avg_sacamp_etrig_lfp = squeeze(nanmean(sacamp_etrig_lfp));
avg_sacamp_etrig_std = squeeze(nanstd(sacamp_etrig_lfp));

znorm_cond_trig_spec = bsxfun(@minus,cond_trig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,length(wfreqs)]));
znorm_cond_trig_spec = bsxfun(@rdivide,znorm_cond_trig_spec,reshape(all_std_spectrum,[length(use_lfps),1,length(wfreqs)]));
avg_cond_trig_spec = squeeze(mean(znorm_cond_trig_spec));
avg_cond_trig_lfp = squeeze(mean(cond_trig_lfp));
std_cond_trig_lfp = squeeze(std(cond_trig_lfp));

stim_start_times = 0.24*(1:16);
% stim_start_inds

%%
% close all
ca = [-0.25 0.75];
xl = [-0.3 0.5];
% for aa = 1:n_prctile_bins-1
%     fprintf('%.2f to %.2f\n',sac_amp_prctiles(aa),sac_amp_prctiles(aa+1));
%     subplot(2,1,1)
%     pcolor(lags/Fsd,wfreqs,squeeze(avg_znorm_sacamp_btrig_spec(aa,:,:))');shading flat
%     caxis(ca);
%     xlim(xl)
%     subplot(2,1,2)
%     shadedErrorBar(lags/Fsd,avg_sacamp_btrig_lfp(aa,:),avg_sacamp_btrig_std(aa,:));
%     xlim(xl)
%     ylim([-10 5]*1e-5);
%     pause
%     clf
% end

close all
    figure
    pcolor(lags/Fsd,wfreqs,squeeze((avg_znorm_sacamp_btrig_spec))');shading flat
    caxis(ca);
    xlim(xl)
    
    figure
    shadedErrorBar(lags/Fsd,avg_sacamp_btrig_lfp,avg_sacamp_btrig_std);
    
    figure
    pcolor(tlags/Fsd,wfreqs,avg_cond_trig_spec');shading flat
    caxis(ca);
%     xlim([0 2])

