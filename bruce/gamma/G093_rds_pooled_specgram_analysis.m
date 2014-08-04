% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G093';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

%%
dsf = 1; %lfps originally sampled at 400Hz
use_lfps = 1:2:96;
lcf = 0.5;

cur_expt_set = [42];

start_buffer = round(Fsd*0.5);
min_trial_dur = 0.5;
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
all_trial_sl = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    trial_sl = [Expts{cur_expt}.Trials(:).sl];
    
    [un_ids,id_inds] = unique(trial_ids);
    use_trials = id_inds(trial_durs(id_inds) >= min_trial_dur);
    
    trial_start_times = trial_start_times(use_trials);
    trial_end_times = trial_end_times(use_trials);
    trial_durs = trial_durs(use_trials);
    trial_sl = trial_sl(use_trials);
    
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times');
    all_trial_durs = cat(1,all_trial_durs,trial_durs');
    all_trial_sl = cat(1,all_trial_sl,trial_sl');
    
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

%% CLASSIFY TRIALS
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));

conditions = [all_trial_sl(:)];
[C,IA,IC] = unique(conditions,'rows');
cond_spec{1} = 'sl'; 
fprintf('%d unique conditions\n',length(C));

%% TRIG AVG WAVELET SPECGRAMS

%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 40;
min_freq = 2; max_freq = 125;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

forlag = round(Fsd*4);
backlag = round(0.2*Fsd);

inc_prev_trial = 1;
nboot = 0;

for ll = 1:length(use_lfps);
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    temp = cwt(all_V(:,ll),scales,'cmor1-1');
    cur_ampgram = abs(temp)';
    ov_avg_ampgram(ll,:) = mean(cur_ampgram);
    ov_std_ampgram(ll,:) = std(cur_ampgram);
    
    for cc = 1:length(C)
        fprintf('Condition %d of %d\n',cc,length(C));
        
        cur_trials = find(IC == cc & all_trial_durs > start_buffer/Fsd);
        [cond_trig_spec(ll,cc,:,:),lags] = get_event_trig_avg(cur_ampgram,all_trial_start_inds(cur_trials),backlag,forlag,nboot,all_trialvec,inc_prev_trial);
        [cond_trig_lfp(ll,cc,:),lags] = get_event_trig_avg(all_V(:,ll),all_trial_start_inds(cur_trials),backlag,forlag,nboot,all_trialvec,inc_prev_trial);        
        
    end
end


%%
save_dir = ['~/Analysis/bruce/' Expt_name];
cd(save_dir);
if ~exist('gamma','dir')
    system('mkdir gamma');
end
cd('gamma');
save rdscond_trig_avg_gamma use_lfps wfreqs lags Fsd cond_* C ov_*

