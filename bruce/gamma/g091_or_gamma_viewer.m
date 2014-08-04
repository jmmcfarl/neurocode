% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G091';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

%%
dsf = 1; %lfps originally sampled at 400Hz
use_lfps = 1:2:96;
lcf = 0.5;

% cur_expt_set = [39:41];
cur_expt_set = [44 45 46 47 48 49];

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
all_trial_nph = [];
all_trial_or = [];
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
%     trial_or = [Expts{cur_expt}.Trials(:).or];
    trial_nph = [Expts{cur_expt}.Trials(:).nph];
    trial_tf = [Expts{cur_expt}.Trials(:).tf];
    
    [un_ids,id_inds] = unique(trial_ids);
    use_trials = id_inds(trial_durs(id_inds) >= 0.5);
    
    trial_start_times = trial_start_times(use_trials);
    trial_end_times = trial_end_times(use_trials);
    trial_durs = trial_durs(use_trials);
%     trial_or = trial_or(use_trials);
    trial_tf = trial_tf(use_trials);
    trial_nph = trial_nph(use_trials);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times');
    all_trial_durs = cat(1,all_trial_durs,trial_durs');
    all_trial_nph = cat(1,all_trial_nph,trial_nph');
    all_trial_tf = cat(1,all_trial_tf,trial_tf');
%     all_trial_or = cat(1,all_trial_or,trial_or');
    
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

%% CLASSIFY TRIALS
used_inds = find(all_tsince_start > start_buffer/Fsd);
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));
all_trial_durs = (all_trial_end_inds - all_trial_start_inds)/Fsd;

% conditions = [all_trial_or(:) all_trial_nph(:)];
conditions = [all_trial_tf(:) all_trial_nph(:)];
cond_spec{1} = 'or'; 
cond_spec{2} = 'nph';
[C,IA,IC] = unique(conditions,'rows');
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
    cur_trials = find(IC == cc & all_trial_durs >= 2);
    for ll = 1:length(use_lfps)
        
        [S_cond(cc,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
        
    end
end

%%
close all
drift_set = find(C(:,2) == 0);
rp_set = find(C(:,2) == 20);
% drift_set = circshift(drift_set,4);
% rp_set = circshift(rp_set,4);
cmap = jet(6);
subplot(1,2,1); hold on
for cc = 1:6
%     shadedErrorBar(f,squeeze(mean(log10(S_cond(drift_set(cc),:,:)),2)),squeeze(std(log10(S_cond(drift_set(cc),:,:)),[],2))/sqrt(length(use_lfps)),{'color',cmap(cc,:)});
    plot(f,squeeze(mean(log10(S_cond(drift_set(cc),:,:)),2)),'color',cmap(cc,:),'linewidth',2);
end
legend('0 Hz','1 Hz','2 Hz','4 Hz','8 Hz','16 Hz');
hold on
%     plot(f,squeeze(mean(log10(S_cond(drift_set(1),:,:)),2)),'k','linewidth',2);
plot(f,squeeze(mean(squeeze(mean(log10(S_cond(rp_set,:,:)))))),'k','linewidth',3);
xlim([0 150]);
ylim([-14 -9.5]);
xlabel('Frequency (Hz)','fontsize',14);
ylabel('Power (dB)','fontsize',14);
subplot(1,2,2); hold on
for cc = 1:6
    plot(f,squeeze(mean(log10(S_cond(rp_set(cc),:,:)),2)),'color',cmap(cc,:),'linewidth',2);
end
xlim([0 150]);
ylim([-14 -9.5]);
xlabel('Frequency (Hz)','fontsize',14);
ylabel('Power (dB)','fontsize',14);
legend('0 Hz','1 Hz','2 Hz','4 Hz','8 Hz','16 Hz');


%%
%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 40;
min_freq = 2; max_freq = 120;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

forlag = round(Fsd*0.75);
backlag = round(0.2*Fsd);
inc_prev_trial = 1;
for ll = 1:length(use_lfps);
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    temp = cwt(all_V(:,ll),scales,'cmor1-1');
    cur_ampgram = abs(temp)';
    ov_avg_ampgram(ll,:) = mean(cur_ampgram);
    ov_std_ampgram(ll,:) = std(cur_ampgram);    
    
    %%
    for cc = 1:n_events
        fprintf('Event type %d of %d\n',cc,n_events);
        use_events = event_types{cc};
        use_events(ismember(use_events,bad_events)) = [];
        [cond_trig_spec(ll,cc,:,:),lags] = get_event_trig_avg(cur_ampgram,event_inds(use_events),backlag,forlag,0,all_trialvec,inc_prev_trial);
        [cond_trig_lfp(ll,cc,:),lags] = get_event_trig_avg(all_V(:,ll),event_inds(use_events),backlag,forlag,0,all_trialvec,inc_prev_trial);
    end
end

%%
close all
cond_num = 7;
cur_t_set = find(IC == cond_num);
comp_t_set = cur_t_set(all_trial_durs(cur_t_set) >= 4);
cur_lfp = 1;
use_t = comp_t_set(1);

temp = cwt(all_V(:,cur_lfp),scales,'cmor1-1');
cur_ampgram = abs(temp)';

cur_t_inds = (all_trial_start_inds(use_t)-backlag):(all_trial_end_inds(use_t)+forlag);
subplot(3,1,1)
plot(all_t_axis(cur_t_inds) - all_t_axis(all_trial_start_inds(use_t)),all_V(cur_t_inds,cur_lfp),'k'); axis tight
subplot(3,1,[2 3])
pcolor(all_t_axis(cur_t_inds) - all_t_axis(all_trial_start_inds(use_t)),wfreqs,cur_ampgram(cur_t_inds,:)');shading flat


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







