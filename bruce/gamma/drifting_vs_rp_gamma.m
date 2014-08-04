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
% all_trial_or = [];
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

%% CLASSIFY TRIALS
start_buffer = round(Fsd*0.5);
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));

all_drifting_grate_trials = find(all_trial_nph == 0);
all_rp_grate_trials = find(all_trial_nph > 0);

%% TRIAL-AVERAGED POWER SPECTRA
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));

params.Fs = Fsd;
params.tapers = [3 5];
movingwin = [1.5 1.5];
sMarkers = [all_trial_start_inds(:)+start_buffer all_trial_end_inds(:)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];


clear S_rp_or S_dg_or
un_tfs = unique(all_trial_tf);
for ii = 1:length(un_tfs)
    fprintf('TF %d of %d\n',ii,length(un_tfs));
    cur_tf_trials = find(all_trial_tf == un_tfs(ii) & used_trial_durs > movingwin(1));
    
    cur_dg_trials = all_drifting_grate_trials(ismember(all_drifting_grate_trials,cur_tf_trials));
    cur_rp_trials = all_rp_grate_trials(ismember(all_rp_grate_trials,cur_tf_trials));
    
    for ll = 1:length(use_lfps)
        
        [S_dg_or(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_dg_trials,:));
        [S_rp_or(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_rp_trials,:));
        
    end
end

%%
cmap = jet(length(un_tfs));

for ii = 1:length(un_tfs)
    subplot(1,2,1)
%     shadedErrorBar(f,squeeze(mean(log(S_dg_or(ii,:,:)),2)),squeeze(std(log(S_dg_or(ii,:,:)),[],2))/sqrt(length(use_lfps)),{'color',cmap(ii,:)});
    plot(f,squeeze(mean(log(S_dg_or(ii,:,:)),2)),'color',cmap(ii,:),'linewidth',2);
    hold on
    
    subplot(1,2,2)
%     shadedErrorBar(f,squeeze(mean(log(S_rp_or(ii,:,:)),2)),squeeze(std(log(S_rp_or(ii,:,:)),[],2))/sqrt(length(use_lfps)),{'color',cmap(ii,:)});
    plot(f,squeeze(mean(log(S_rp_or(ii,:,:)),2)),'color',cmap(ii,:),'linewidth',2);
    hold on
    
end

subplot(1,2,1)
xlim([1 150]); %ylim([-30 -22])
% set(gca,'xscale','log');
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Log power','fontsize',16)
legend('0','1','2','4','8.3','16.7');
title('Drifting grating','fontsize',16)

subplot(1,2,2)
xlim([1 150]); %ylim([-30 -22])
% set(gca,'xscale','log');
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Log power','fontsize',16)
legend('0','1','2','4','8.3','16.7');
title('Random-phase grating','fontsize',16)

%%  FOR SPECTROGRAMS
% params.Fs = Fsd;
% params.tapers = [2 3];
% movingwin = [0.1 0.1];
% sMarkers = [all_trial_start_inds(:) all_trial_end_inds(:)];
% used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
% params.err = [0];
% params.trialave = 1;
% trig_win = [0 3.75];
% 
% clear S_rp_or S_dg_or
% un_tfs = unique(all_trial_tf);
% for ii = 1:length(un_tfs)
%     fprintf('TF %d of %d\n',ii,length(un_tfs));
%     cur_tf_trials = find(all_trial_tf == un_tfs(ii) & used_trial_durs > trig_win(2));
%     
%     cur_dg_trials = all_drifting_grate_trials(ismember(all_drifting_grate_trials,cur_tf_trials));
%     cur_rp_trials = all_rp_grate_trials(ismember(all_rp_grate_trials,cur_tf_trials));
%     
%     cur_dg_trig_inds = all_trial_start_inds(cur_dg_trials);
%     cur_rp_trig_inds = all_trial_start_inds(cur_rp_trials);
%     
%     for ll = 1:length(use_lfps)
%         
%         [S,t,f]=mtspecgramtrigc( all_V(:,ll),cur_rp_trig_inds/Fsd,trig_win,movingwin,params);
% 
% 
% [S_dg_or(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_dg_trials,:));
%         [S_rp_or(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_rp_trials,:));
%         
%     end
% end

%% COMPUTE JUMP INDS FOR RP TRIALS
for ii = 1:length(un_tfs)
    cur_reltimeset = 1/un_tfs(ii):1/un_tfs(ii):(4-1/un_tfs(ii));
    cur_relindset = round(cur_reltimeset*Fsd);
    
    cur_start_ind_set{ii} = [];
    cur_tf_trials = find(all_trial_tf == un_tfs(ii) & used_trial_durs > movingwin(1));
    cur_rp_trials = all_rp_grate_trials(ismember(all_rp_grate_trials,cur_tf_trials));
    for tt = 1:length(cur_rp_trials)
        cur_jump_inds = all_trial_start_inds(cur_rp_trials(tt)) + cur_relindset;
%         cur_start_times = all_trial_start_times(cur_rp_trials(tt)) + cur_reltimeset;
%         bad = find(cur_start_times(2:end) > all_trial_start_times(cur_rp_trials(tt)));
%         cur_start_times(bad) = [];
%         cur_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_start_times));
        cur_start_ind_set{ii} = cat(1,cur_start_ind_set{ii},cur_jump_inds');
    end
end

%% TRIG AVG WAVELET SPECGRAMS

%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 35;
min_freq = 2; max_freq = 100;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

forlag = round(Fsd*4);
backlag = round(0.2*Fsd);
jforlag = round(Fsd*1);
jbacklag = round(0.2*Fsd);


un_tfs = unique(all_trial_tf);
un_tfs(1) = [];

for ll = 1:length(use_lfps);
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    temp = cwt(all_V(:,ll),scales,'cmor1-1');
    cur_ampgram = abs(temp)';
    
    base_trials = find(all_trial_tf == 0 & all_trial_nph == 0);
    base_inds = find(ismember(all_trialvec,base_trials));
    
    avg_pow_base = mean(cur_ampgram(base_inds,:,:));
    std_pow_base = std(cur_ampgram(base_inds,:,:));
    
    cur_ampgram_norm = bsxfun(@minus,cur_ampgram,avg_pow);
    cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram_norm,std_pow);

%     un_wfreqs = linspace(wfreqs(1),wfreqs(end),100);
%     avg_pow_int = interp1(wfreqs,avg_pow,un_wfreqs);
%     powfun = @(a,x)(a(1)*x.^a(2));
%     BETA = nlinfit(un_wfreqs,avg_pow_int,powfun,[max(avg_pow) -2]);
%     powfit = un_wfreqs.^BETA(2)*BETA(1);
%     interp_powfit = interp1(un_wfreqs,powfit,wfreqs);
    
%     cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram,interp_powfit);
    % [trig_avg,lags] = get_event_trig_avg(cur_ampgram_norm,dg_tstart_inds,backlag,forlag,[]);
    % [trig_avg2,lags] = get_event_trig_avg(cur_ampgram_norm,rp_tstart_inds,backlag,forlag,[]);
    
    %%
    % clear or_trig*
    for ii = 1:length(un_tfs)
        fprintf('Orientation %d of %d\n',ii,length(un_tfs));
        cur_tf_trials = find(all_trial_tf == un_tfs(ii) & used_trial_durs > movingwin(1));
        cur_dg_trials = all_drifting_grate_trials(ismember(all_drifting_grate_trials,cur_tf_trials));
        cur_rp_trials = all_rp_grate_trials(ismember(all_rp_grate_trials,cur_tf_trials));
        
        [dg_tf_trig_avg(ll,ii,:,:),lags] = get_event_trig_avg(cur_ampgram_norm,all_trial_start_inds(cur_dg_trials),backlag,forlag,[]);
        [rp_tf_trig_avg(ll,ii,:,:),lags] = get_event_trig_avg(cur_ampgram_norm,all_trial_start_inds(cur_rp_trials),backlag,forlag,[]);
        
        [rp_jumptrig_avg(ll,ii,:,:),jlags] = get_event_trig_avg(cur_ampgram_norm,cur_start_ind_set{ii},jbacklag,jforlag,[]);
        
    end
    
end





