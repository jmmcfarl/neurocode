% close all
clear all


dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G095';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

%%
dsf = 1; %lfps originally sampled at 400Hz
use_lfps = 1:4:96;
lcf = 0.5;

min_trial_dur = 1;
cur_expt_set = [2:7];

stim_dt = 0.5;
%%
fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_tsince_start = [];
all_exptvec = [];
all_trial_result = [];
all_trialvec = [];
all_trial_durs = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_V = [];
all_trial_exptvec = [];
all_stim_times = [];
all_stim_sl = [];
all_stim_jv = [];
all_stim_jcomp = [];
all_stim_tonset = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    trial_result = [Expts{cur_expt}.Trials(:).Result];
    
%     [un_ids,id_inds] = unique(trial_ids);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times');
    all_trial_durs = cat(1,all_trial_durs,trial_durs');
    all_trial_result = cat(1,all_trial_result,trial_result');
    all_trial_exptvec = cat(1,all_trial_exptvec,ee*ones(length(trial_ids),1));
    
    n_trials = length(trial_ids);
    
    for tt = 1:n_trials
        cur_tstart = all_trial_start_times(tt);
        cur_tend = all_trial_end_times(tt);
        cur_stim_times = (cur_tstart:stim_dt:cur_tend)';
        cur_nstims = length(cur_stim_times);
        cur_stim_sl = [Expts{cur_expt}.Trials(tt).sl];
        cur_stim_jv = [Expts{cur_expt}.Trials(tt).jv];
        cur_stim_jcomp = [Expts{cur_expt}.Trials(tt).jcomp];
        cur_tonset = zeros(size(cur_stim_times));
        cur_tonset(1) = 1;
            
        all_stim_times = [all_stim_times;  cur_stim_times(:)];
        all_stim_sl = [all_stim_sl;  repmat(cur_stim_sl,cur_nstims,1)];
        all_stim_jv = [all_stim_jv;  repmat(cur_stim_jv,cur_nstims,1)];
        all_stim_jcomp = [all_stim_jcomp;  repmat(cur_stim_jcomp,cur_nstims,1)];
        all_stim_tonset = [all_stim_tonset; cur_tonset];
    end
    
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


%%
event_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_stim_times));
event_types = zeros(size(event_inds));
bad_events = find(isnan(event_inds));
event_inds(isnan(event_inds)) = 1; %just to avoid errors

start_buffer = round(Fsd*0.3);
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));


used_inds = find(all_tsince_start > 0.2);
%%
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));

params.Fs = Fsd;
params.tapers = [3 5];
movingwin = [2 2];
sMarkers = [all_trial_start_inds(:)+start_buffer all_trial_end_inds(:)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];

t_onsets = find(all_stim_tonset == 1);
% drift_trials = find(all_stim_jv(t_onsets) == 5);
% tf_trials = find(all_stim_jv(t_onsets) == 0);
drift_trials = find(all_stim_jv(t_onsets) == 5 & all_stim_sl(t_onsets) == 0);
tf_trials = find(all_stim_jv(t_onsets) == 0 & all_stim_sl(t_onsets) == 0);

drift_trials(used_trial_durs(drift_trials) < movingwin(1)) = [];
tf_trials(used_trial_durs(tf_trials) < movingwin(1)) = [];
clear S_*
for ll = 1:length(use_lfps)
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    [S_drift(ll,:), f]= mtspectrumc_unequal_length_trials(all_V(:,ll), movingwin, params, sMarkers(drift_trials,:));
    [S_tf(ll,:), f]= mtspectrumc_unequal_length_trials(all_V(:,ll), movingwin, params, sMarkers(tf_trials,:));    
end

%%
lS_drift = log10(S_drift);
lS_tf = log10(S_tf);
h1=plot(f,squeeze(mean(lS_drift)),'k','linewidth',2);
hold on
h2=plot(f,squeeze(mean(lS_tf)),'r','linewidth',2);
xlim([0 120])
ylim([-13 -9.5])
legend([h1 h2],{'Uniform drift','Uniform TF'});
xlabel('Frequency (Hz)','fontsize',14);
ylabel('Log power','fontsize',14);

%%
%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 35;
min_freq = 2; max_freq = 100;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

forlag = round(Fsd*0.5);
backlag = round(0.2*Fsd);
inc_prev_trial = 1;
% for ll = 1:length(use_lfps);
for ll = 2
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    temp = cwt(all_V(:,ll),scales,'cmor1-1');
    cur_ampgram = abs(temp)';
    ov_avg_ampgram(ll,:) = mean(cur_ampgram(used_inds,:));
    ov_std_ampgram(ll,:) = std(cur_ampgram(used_inds,:));
    
    
%     un_wfreqs = linspace(wfreqs(1),wfreqs(end),100);
%     avg_pow_int = interp1(wfreqs,ov_avg_ampgram(ll,:),un_wfreqs);
%     
%     powfun = @(a,x)(a(1)*x.^a(2));
%     BETA = nlinfit(un_wfreqs,avg_pow_int,powfun,[max(ov_avg_ampgram(ll,:)) -2]);
%     powfit = un_wfreqs.^BETA(2)*BETA(1);
%     interp_powfit = interp1(un_wfreqs,powfit,wfreqs);
    
%     cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram,interp_powfit);
%     cur_ampgram_norm = bsxfun(@minus,cur_ampgram,ov_avg_ampgram(ll,:));
    cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram,ov_std_ampgram(ll,:));
    
    cur_events = find(all_stim_tonset == 0 & all_stim_sl == 0);
    
    [cond_trig_spec(ll,:,:),lags] = get_event_trig_avg(cur_ampgram_norm,event_inds(cur_events),backlag,forlag,0,all_trialvec,inc_prev_trial);
    [cond_trig_lfp(ll,:),lags] = get_event_trig_avg(all_V(:,ll),event_inds(cur_events),backlag,forlag,0,all_trialvec,inc_prev_trial);
    
    %%
%     for cc = 1:n_events
%         fprintf('Event type %d of %d\n',cc,n_events);
%         use_events = event_types{cc};
%         use_events(ismember(use_events,bad_events)) = [];
%         [cond_trig_spec(ll,cc,:,:),lags] = get_event_trig_avg(cur_ampgram,event_inds(use_events),backlag,forlag,0,all_trialvec,inc_prev_trial);
%         [cond_trig_lfp(ll,cc,:),lags] = get_event_trig_avg(all_V(:,ll),event_inds(use_events),backlag,forlag,0,all_trialvec,inc_prev_trial);
%     end
end

%%
