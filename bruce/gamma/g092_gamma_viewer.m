% close all
clear all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G092';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

%%
dsf = 1; %lfps originally sampled at 400Hz
use_lfps = 1:4:96;
lcf = 0.5;

cur_expt_set = [28:44];
min_trial_dur = 1;
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
all_trial_exptvec = [];
all_trial_nph = [];
all_trial_jv = [];
all_trial_st = [];
all_trial_sl = [];
all_trial_opt = [];
all_trial_tf = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
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
start_buffer = round(Fsd*0.5);

%     %% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
%     cd(data_dir)
%     emfile = ['jbe' Expt_name '.em.mat'];
%     load(emfile);
%     
%     all_eye_vals = [];
%     all_eye_speed = [];
%     all_eye_ts = [];
%     all_eye_blockvec = [];
%     eye_smooth = 3;
%     for ee = 1:length(cur_expt_set);
%         fprintf('Loading ET data for expt block %d of %d\n',ee,length(cur_expt_set));
%         cur_set = find(all_exptvec==ee);
%         [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));
%         
%         eye_dt = Expt.Header.CRrates(1);
%         eye_fs = 1/eye_dt;
%         lEyeXY = eye_vals_interp(:,1:2);
%         rEyeXY = eye_vals_interp(:,3:4);
%         
%         %slight smoothing before computing speed
%         sm_avg_eyepos = lEyeXY; eye_vel = lEyeXY; %initialization
%         sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),eye_smooth);
%         sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),eye_smooth);
%         eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
%         eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
%         
%         eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
%         
%         all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
%         all_eye_speed = [all_eye_speed; eye_speed];
%         all_eye_ts = [all_eye_ts; eye_ts_interp'];
%         all_eye_blockvec = [all_eye_blockvec; ee*ones(size(eye_speed))];
%     end
%     
%     back_pts = 1 + find(diff(all_eye_ts) <= 0);
%     double_samples = [];
%     for i = 1:length(back_pts)
%         next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
%         double_samples = [double_samples back_pts(i):next_forward];
%     end
%     all_eye_ts(double_samples) = [];
%     all_eye_speed(double_samples) = [];
%     all_eye_vals(double_samples,:) = [];
%     
%     interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis);
%     interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);
%     
%     %%
%     sac_thresh = 10;
%     peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
%     saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);
%     
%     peri_thresh = 3; %threshold eye speed for defining saccade boundary inds
%     thresh_cross_up = 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh);
%     thresh_cross_down = 1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh);
%     sac_start_inds = nan(size(saccade_inds));
%     sac_stop_inds = nan(size(saccade_inds));
%     for ii = 1:length(saccade_inds)
%         next_tc = find(thresh_cross_down > saccade_inds(ii),1,'first');
%         if ~isempty(next_tc)
%             sac_stop_inds(ii) = thresh_cross_down(next_tc);
%         end
%         prev_tc = find(thresh_cross_up < saccade_inds(ii),1,'last');
%         if ~isempty(prev_tc)
%             sac_start_inds(ii) = thresh_cross_up(prev_tc);
%         end
%         
%     end
%     
%     %get rid of double-peaks
%     min_isi = 0.05; max_isi = Inf;
%     isis = [Inf; diff(sac_start_inds)]/eye_fs;
%     bad_isis = (isis < min_isi | isis > max_isi);
%     bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
%     saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];
%     
%     saccade_times = all_eye_ts(saccade_inds);
%     sac_start_times = all_eye_ts(sac_start_inds);
%     sac_stop_times = all_eye_ts(sac_stop_inds);
%     sac_durs = sac_stop_times - sac_start_times;
%     
%     sac_dbuff = round(0.005/eye_dt);
%     pre_inds = saccade_inds - sac_dbuff;
%     pre_inds(pre_inds < 1) = 1;
%     sac_pre_pos = all_eye_vals(pre_inds,:);
%     post_inds = saccade_inds + sac_dbuff;
%     post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
%     sac_post_pos = all_eye_vals(post_inds,:);
%     
%     %use only left-eye signal here
%     sac_delta_pos = sac_post_pos(:,1:2) - sac_pre_pos(:,1:2);
%     sac_amps = sqrt(sum(sac_delta_pos.^2,2));
%     sac_dirs = atan2(sac_delta_pos(:,2),sac_delta_pos(:,1));
%     
%     temp = ones(length(saccade_times),1);
%     saccades = struct('peak_time',mat2cell(saccade_times,temp),'start_time',mat2cell(sac_start_times,temp),...
%         'stop_time',mat2cell(sac_stop_times,temp),'isi',mat2cell(isis,temp),...
%         'duration',mat2cell(sac_durs,temp),'amplitude',mat2cell(sac_amps,temp),'direction',mat2cell(sac_dirs,temp),...
%         'pre_Lx',mat2cell(sac_pre_pos(:,1),temp),'post_Lx',mat2cell(sac_post_pos(:,1),temp),...
%         'pre_Ly',mat2cell(sac_pre_pos(:,2),temp),'post_Ly',mat2cell(sac_post_pos(:,2),temp),...
%         'pre_Rx',mat2cell(sac_pre_pos(:,3),temp),'post_Rx',mat2cell(sac_post_pos(:,3),temp),...
%         'pre_Ry',mat2cell(sac_pre_pos(:,4),temp),'post_Ry',mat2cell(sac_post_pos(:,4),temp));

%% CLASSIFY TRIALS
used_inds = find(all_tsince_start > start_buffer/Fsd);
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));
all_trial_durs = (all_trial_end_inds - all_trial_start_inds)/Fsd;

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

%%
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));
start_buffer = round(0.2*Fsd);
params.Fs = Fsd;
params.tapers = [3 5];
movingwin = [1.6 1.6];
sMarkers = [all_trial_start_inds(:)+start_buffer all_trial_end_inds(:)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];
clear S_cond
for cc = 1:size(C,1)
    fprintf('Cond %d of %d\n',cc,size(C,1));
    cur_trials = find(IC == cc & all_trial_durs >= 2);
    for ll = 1:length(use_lfps)
        
        [S_cond(cc,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
        
    end
end

%%
close all
lS_cond = log10(S_cond);
h1=plot(f,squeeze(mean(lS_cond(2,:,:),2)),'k','linewidth',2);
hold on
% h2=plot(f,squeeze(mean(lS_cond(3,:,:),2)),'k','linewidth',2);
% h3=plot(f,squeeze(mean(lS_cond(4,:,:),2)),'k','linewidth',1);
h4=plot(f,squeeze(mean(lS_cond(18,:,:),2)),'r','linewidth',2);
h5=plot(f,squeeze(mean(lS_cond(19,:,:),2)),'g','linewidth',2);
h6=plot(f,squeeze(mean(lS_cond(20,:,:),2)),'b','linewidth',2);
xlim([0 120]);
legend([h1 h4 h5 h6],{'Drifting grating','RLS jv=0.5','RLS jv=1.0','RLS jv=2.1'});
xlabel('Frequency (Hz)','fontsize',14);
ylabel('Log power','fontsize',14);

% h1=shadedErrorBar(f,squeeze(mean(lS_cond(2,:,:),2)),squeeze(std(lS_cond(2,:,:),[],2))/sqrt(length(use_lfps)),{'color','k'});
% hold on
% h1=shadedErrorBar(f,squeeze(mean(lS_cond(3,:,:),2)),squeeze(std(lS_cond(3,:,:),[],2))/sqrt(length(use_lfps)),{'color','k'});
% h1=shadedErrorBar(f,squeeze(mean(lS_cond(4,:,:),2)),squeeze(std(lS_cond(4,:,:),[],2))/sqrt(length(use_lfps)),{'color','k'});
% h2=shadedErrorBar(f,squeeze(mean(lS_cond(18,:,:),2)),squeeze(std(lS_cond(18,:,:),[],2))/sqrt(length(use_lfps)),{'color','r'});
% h3=shadedErrorBar(f,squeeze(mean(lS_cond(19,:,:),2)),squeeze(std(lS_cond(19,:,:),[],2))/sqrt(length(use_lfps)),{'color','g'});
% h4=shadedErrorBar(f,squeeze(mean(lS_cond(20,:,:),2)),squeeze(std(lS_cond(20,:,:),[],2))/sqrt(length(use_lfps)),{'color','b'});
% h3=shadedErrorBar(f,squeeze(mean(lS_cond(21,:,:),2)),squeeze(std(lS_cond(21,:,:),[],2))/sqrt(length(use_lfps)),{'color','c'});
% h3=shadedErrorBar(f,squeeze(mean(lS_cond(22,:,:),2)),squeeze(std(lS_cond(22,:,:),[],2))/sqrt(length(use_lfps)),{'color','m'});
% legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Drifting grating','RLS jv=0.5','RLS jv=1.0','RLS jv=2.1'})
% xlim([0 120]);
%%
% close all
% rp_grate = find(C(:,1) == 36 & C(:,2) == 3 & C(:,3) == 0 & C(:,4) == 0);
% drift_grate = find(C(:,1) == 36 & C(:,2) == 3 & C(:,3) == 0 & C(:,4) == 4.17);
% randjump_grate = find(C(:,1) == 36 & C(:,2) == 3 & C(:,3) == 24 & C(:,4) == 0 & C(:,5) == 360);
% 
% figure;hold on
% plot(f,squeeze(mean(squeeze(log10(S_cond(rp_grate,:,:))))),'k','linewidth',2);
% plot(f,squeeze(mean(squeeze(log10(S_cond(drift_grate,:,:))))),'r','linewidth',2);
% plot(f,squeeze(mean(squeeze(log10(S_cond(randjump_grate,:,:))))),'b','linewidth',2);
% xlim([0 150]);
% ylim([-14 -9.5]);
% xlabel('Frequency (Hz)','fontsize',14);
% ylabel('Power (dB)','fontsize',14);

%%
%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 40;
min_freq = 2; max_freq = 120;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

inc_prev_trial = 1;
nboot = 0;

forlag = round(Fsd*0.75);
backlag = round(0.2*Fsd);
inc_prev_trial = 1;
ll = 1;
% for ll = 1:length(use_lfps);
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    temp = cwt(all_V(:,ll),scales,'cmor1-1');
    cur_ampgram = abs(temp)';
    ov_avg_ampgram(ll,:) = mean(cur_ampgram(used_inds,:));
    ov_std_ampgram(ll,:) = std(cur_ampgram(used_inds,:));
    
    for cc = 1:length(C)
        fprintf('Condition %d of %d\n',cc,length(C));
        
        cur_trials = find(IC == cc & all_trial_durs > start_buffer/Fsd);
        [cond_trig_spec(ll,cc,:,:),lags] = get_event_trig_avg(cur_ampgram,all_trial_start_inds(cur_trials),backlag,forlag,nboot,all_trialvec,inc_prev_trial);
        [cond_trig_lfp(ll,cc,:),lags] = get_event_trig_avg(all_V(:,ll),all_trial_start_inds(cur_trials),backlag,forlag,nboot,all_trialvec,inc_prev_trial);        
        
    end
% end

%%
cd ~/Analysis/bruce/ 
load('gamma_powfits');
% znorm_cond_spec = bsxfun(@minus,cond_trig_spec,reshape(all_avg_spectrum,[length(use_lfps),1,1,length(wfreqs)]));
% znorm_cond_spec = bsxfun(@rdivide,znorm_cond_spec,reshape(all_std_spectrum,[length(use_lfps),1,1,length(wfreqs)]));
zamp_gram = bsxfun(@minus,cur_ampgram,all_avg_spectrum(ll,:));
zamp_gram = bsxfun(@rdivide,zamp_gram,all_std_spectrum(ll,:));

%% DRIFT GRATING EXAMPLE
close all
cond_num = 3;
cur_t_set = find(IC == cond_num);
comp_t_set = cur_t_set(all_trial_durs(cur_t_set) >= 4);
cur_lfp = 1;
use_t = comp_t_set(5);
cur_t_inds = (all_trial_start_inds(use_t)-backlag):(all_trial_end_inds(use_t)+forlag);

cur_sacs = saccade_times(saccade_times >= all_t_axis(cur_t_inds(1)) & saccade_times < all_t_axis(cur_t_inds(end)));
cur_sacs = cur_sacs - all_t_axis(all_trial_start_inds(use_t));

cur_eye_samps = find(all_eye_ts >= all_t_axis(cur_t_inds(1)) & all_eye_ts <= all_t_axis(cur_t_inds(end)));
subplot(4,1,1)
plot(all_eye_ts(cur_eye_samps)- all_t_axis(all_trial_start_inds(use_t)),all_eye_speed(cur_eye_samps));
axis tight
xlabel('Time since trial onset (s)','fontsize',20);
ylabel('Eye speed (deg/sec)','fontsize',20);

subplot(4,1,2)
plot(all_t_axis(cur_t_inds) - all_t_axis(all_trial_start_inds(use_t)),all_V(cur_t_inds,cur_lfp),'k'); axis tight
ylim([-2 2]*1e-4);
yl = ylim();
for ii = 1:length(cur_sacs)
    line(cur_sacs([ii ii]),yl,'color','r','linestyle','--');
end
xlabel('Time since trial onset (s)','fontsize',20);
ylabel('LFP amplitude (V)','fontsize',20);

subplot(4,1,[3 4])
pcolor(all_t_axis(cur_t_inds) - all_t_axis(all_trial_start_inds(use_t)),wfreqs,zamp_gram(cur_t_inds,:)');shading flat
caxis([-1 3])
yl = ylim();
for ii = 1:length(cur_sacs)
    line(cur_sacs([ii ii]),yl,'color','w','linestyle','--');
end
xlabel('Time since trial onset (s)','fontsize',20);
ylabel('Frequency (Hz)','fontsize',20);

%% JUMP GRATING EXAMPLE
close all
cond_num = 12;
cur_t_set = find(IC == cond_num);
comp_t_set = cur_t_set(all_trial_durs(cur_t_set) >= 4);
cur_lfp = 1;
use_t = comp_t_set(2);
stim_times = [0:0.24:4];
cur_t_inds = (all_trial_start_inds(use_t)-backlag):(all_trial_end_inds(use_t)+forlag);

xl = [1 3];

subplot(3,1,1)
plot(all_t_axis(cur_t_inds) - all_t_axis(all_trial_start_inds(use_t)),all_V(cur_t_inds,cur_lfp),'k'); axis tight
ylim([-2 2]*1e-4);
yl = ylim();
for ii = 1:length(stim_times)
    line(stim_times([ii ii]),yl,'color','k');
end
for ii = 1:length(cur_sacs)
    line(cur_sacs([ii ii]),yl,'color','r','linestyle','--');
end
xlim(xl);
xlabel('Time since trial onset (s)','fontsize',14);
ylabel('LFP amplitude (V)','fontsize',14);

subplot(3,1,[2 3])
pcolor(all_t_axis(cur_t_inds) - all_t_axis(all_trial_start_inds(use_t)),wfreqs,zamp_gram(cur_t_inds,:)');shading flat
caxis([-1 4])
yl = ylim();
for ii = 1:length(stim_times)
    line(stim_times([ii ii]),yl,'color','w','linestyle','-');
end
xlim(xl);
xlabel('Time since trial onset (s)','fontsize',14);
ylabel('Frequency (Hz)','fontsize',14);

% for ii = 1:length(cur_sacs)
%     line(cur_sacs([ii ii]),yl,'color','w','linestyle','--');
% end

%%
close all
cond_num = 8;
cur_t_set = find(IC == cond_num);
comp_t_set = cur_t_set(all_trial_durs(cur_t_set) >= 4);
cur_lfps = 1:8:length(use_lfps);
space = 2e-4;
cmap = jet(length(cur_lfps));
use_t = comp_t_set(2);
backlag = round(0.5*Fsd);

cnt = 0;
cur_t_inds = all_trial_start_inds(use_t):all_trial_end_inds(use_t);
for ii = 1:length(cur_lfps)
    plot(all_t_axis(cur_t_inds) - all_t_axis(cur_t_inds(1)) - backlag/Fsd,all_V(cur_t_inds,cur_lfps(ii))+cnt,'color',cmap(ii,:));
    hold on
    cnt = cnt + space;
end







