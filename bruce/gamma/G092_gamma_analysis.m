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

min_trial_dur = 0.5;
cur_expt_set = [28:44];

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

%% CLASSIFY TRIALS
static_grate_trials = find(all_trial_opt == 36 & all_trial_st == 3 & all_trial_sl == 0 ...
    & all_trial_nph == 360 & all_trial_tf == 0);

grate_drift_trials = find(all_trial_opt == 36 & all_trial_st == 3 & all_trial_sl == 0 ...
    & all_trial_nph == 0 & all_trial_tf > 0);

grate_cp_trials = find(all_trial_opt == 52 & all_trial_st == 3 & all_trial_sl == 0 ...
    & all_trial_nph == 0 & all_trial_tf > 0);

grate_drift_detjump_trials = find(all_trial_opt == 36 & all_trial_st == 3 & all_trial_sl > 0 & ...
    all_trial_nph == 0 & all_trial_tf > 0);

grate_drift_randjump_trials = find(all_trial_opt == 36 & all_trial_st == 3 & all_trial_sl > 0 & ...
    all_trial_nph == 360 & all_trial_tf > 0);

grate_stat_randjump_trials = find(all_trial_opt == 36 & all_trial_st == 3 & all_trial_sl > 0 & ...
    all_trial_nph == 360 & all_trial_tf == 0);

rls_stat_randjump_trials = find(all_trial_opt == 36 & all_trial_st == 15 & all_trial_sl > 1 & ...
    all_trial_nph == 0 & all_trial_tf > 0 & all_trial_jv == 0);

rls_drift_trials = find(all_trial_opt == 36 & all_trial_st == 15 & all_trial_sl == 1 & ...
    all_trial_nph == 0 & all_trial_tf > 0 & all_trial_jv > 0);

rls_cp_trials = find(all_trial_opt == 52 & all_trial_st == 15 & all_trial_sl == 1 & ...
    all_trial_nph == 0 & all_trial_tf > 0 & all_trial_jv == 0);

all_dsdrift_trials = find(all_trial_opt == 36 & all_trial_st == 3 & all_trial_sl > 0 & all_trial_nph == 0);

all_trial_jv(all_trial_st == 3) = 0;
cond_mat = [all_trial_opt(:) all_trial_st(:) all_trial_sl(:) all_trial_tf(:) all_trial_nph(:) all_trial_jv(:)];
[C,IA,IC] = unique(cond_mat,'rows');
tabulate(IC)


%% TRIAL-AVERAGED POWER SPECTRA
start_buffer = round(Fsd*0.5);
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));

addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));

params.Fs = Fsd;
params.tapers = [3 5];
movingwin = [3.25 3.25];
sMarkers = [all_trial_start_inds(:)+start_buffer all_trial_end_inds(:)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];

clear S_cond
for cc = 1:length(C)
    fprintf('Condition %d of %d\n',cc,length(C));
    cur_trials = find(IC == cc & used_trial_durs > movingwin(2));
    fprintf('Using %d trials\n',length(cur_trials));
    for ll = 1:length(use_lfps)
        [S_cond(ll,cc,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
    end
end

%%
close all
for cc = 2:length(C)
    fprintf('opt: %d  st: %d  sl: %.2f  tf: %.2f  nph: %d  jv: %.2f\n',C(cc,:));
    shadedErrorBar(f,squeeze(nanmean(log10(S_cond(:,cc,:)))),squeeze(nanstd(log10(S_cond(:,cc,:))))/sqrt(length(use_lfps)),{'color','r'});
    hold on
    plot(f,squeeze(nanmean(log10(S_cond(:,1,:)))),'color','k','linewidth',2)
    xlabel('Frequency (Hz)','fontsize',14);
    ylabel('Log power (dB)','fontsize',14);
    xlim([0 150]);
    pause
    clf
end

%%
%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 35;
min_freq = 2; max_freq = 100;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

forlag = round(Fsd*4);
backlag = round(0*Fsd);
jforlag = round(Fsd*1);
jbacklag = round(0.2*Fsd);

for ll = 1:length(use_lfps);
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    temp = cwt(all_V(:,ll),scales,'cmor1-1');
    cur_ampgram = abs(temp)';
    avg_pow = mean(cur_ampgram);
    std_pow = std(cur_ampgram);
    
%     base_trials = static_grate_trials;
%     base_inds = find(ismember(all_trialvec,base_trials));
%     
%     avg_pow_base = mean(cur_ampgram(base_inds,:,:));
%     std_pow_base = std(cur_ampgram(base_inds,:,:));
%     
%     cur_ampgram_norm = bsxfun(@minus,cur_ampgram,avg_pow_base);
%     cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram_norm,std_pow_base);

    
    un_wfreqs = linspace(wfreqs(1),wfreqs(end),100);
    avg_pow_int = interp1(wfreqs,avg_pow,un_wfreqs);
    
    powfun = @(a,x)(a(1)*x.^a(2));
    BETA = nlinfit(un_wfreqs,avg_pow_int,powfun,[max(avg_pow) -2]);
    powfit = un_wfreqs.^BETA(2)*BETA(1);
    interp_powfit = interp1(un_wfreqs,powfit,wfreqs);
    
    cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram,interp_powfit);

    %%
    for cc = 1:length(C)
        fprintf('Condition %d of %d\n',cc,length(C));
        cur_trials = find(IC == cc);
        fprintf('Using %d trials\n',length(cur_trials));
        [trig_spec_cond(ll,cc,:,:),lags] = get_event_trig_avg(cur_ampgram_norm,all_trial_start_inds(cur_trials),backlag,forlag,0,all_trialvec);
        [trig_sig_cond(ll,cc,:),lags] = get_event_trig_avg(all_V(:,ll),all_trial_start_inds(cur_trials),backlag,forlag,0,all_trialvec);
    end
    
end

avg_trig_spec_cond = squeeze(mean(trig_spec_cond));
stpt = round(0.96*Fsd);
ept = stpt + round(Fsd*0.96*3);
windur = round(Fsd*0.96);
avg_trig_spec_reshape = permute(avg_trig_spec_cond,[2 1 3]);
avg_trig_spec_reshape = reshape(avg_trig_spec_reshape(stpt+1:ept,:,:),[windur 3 length(C) length(wfreqs)]);
avg_trig_spec_win = squeeze(mean(avg_trig_spec_reshape,2));

avg_trig_sig_cond = squeeze(mean(trig_sig_cond));
avg_trig_sig_reshape = permute(avg_trig_sig_cond,[2 1]);
avg_trig_sig_reshape = reshape(avg_trig_sig_reshape(stpt+1:ept,:),[windur 3 length(C)]);
avg_trig_sig_win = squeeze(mean(avg_trig_sig_reshape,2));

sem_trig_sig_reshape = permute(trig_sig_cond,[3 1 2]);
sem_trig_sig_reshape = reshape(sem_trig_sig_reshape(stpt+1:ept,:),[windur 3 length(use_lfps) length(C)]);
sem_trig_sig_reshape = squeeze(mean(sem_trig_sig_reshape,2));
sem_trig_sig_win = squeeze(std(sem_trig_sig_reshape,[],2))/sqrt(length(use_lfps));

%%
close all
clc
for cc = 1:length(C)
    fprintf('opt: %d  st: %d  sl: %.2f  tf: %.2f  nph: %d  jv: %.2f\n',C(cc,:));
%     plot(f,squeeze(mean(log(S_cond(:,1,:)))),'k','linewidth',2);hold on
%     plot(f,squeeze(mean(log(S_cond(:,cc,:)))),'r')

subplot(2,1,1)
    pcolor(lags(1:windur)/Fsd,wfreqs,squeeze(avg_trig_spec_win(:,cc,:))');shading flat
%     caxis([-1 1])
caxis([0.5 1.7])
xlim([0 windur/Fsd]);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
    subplot(2,1,2)
shadedErrorBar(lags(1:windur)/Fsd,avg_trig_sig_win(:,cc),sem_trig_sig_win(:,cc));
xlim([0 windur/Fsd]);
xlabel('Time (s)'); ylabel('Amplitude (V)');
xl = xlim();
line(xl,[0 0],'color','k','linestyle','--');

% %     pcolor(lags/Fsd,wfreqs,squeeze(avg_trig_spec_cond(cc,:,:))');shading flat
% %     caxis([-1 1])

pause
    clf
end

%%
%COMPARE DRIFTING AND d-JUMPING GRATINGS

fprintf('TF %d of %d\n',ii,length(un_tfs));
cur_tf_trials = find(all_trial_tf == 4.17 & used_trial_durs > movingwin(1));
% cur_tf_trials = find(all_trial_tf == 8.33 & used_trial_durs > movingwin(1));

cur_dg_trials = grate_drift_trials(ismember(grate_drift_trials,cur_tf_trials));
cur_dg_cp_trials = grate_cp_trials(ismember(grate_cp_trials,cur_tf_trials));
for ll = 1:length(use_lfps)
    [S_dg(ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_dg_trials,:));
    [S_cpgrate(ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_dg_cp_trials,:));
end

cur_sl_trials = find(all_trial_sl == 24 & used_trial_durs > movingwin(1));
cur_grate_randjump_trials = grate_stat_randjump_trials(ismember(grate_stat_randjump_trials,cur_sl_trials));
for ll = 1:length(use_lfps)
    [S_gr_randjump(ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_grate_randjump_trials,:));
end

cur_jvtrials = find(all_trial_jv == 1 & used_trial_durs > movingwin(1));
cur_drls_trials = rls_drift_trials(ismember(rls_drift_trials,cur_jvtrials));
cur_tftrials = find(all_trial_tf == 4.17 & used_trial_durs > movingwin(1));
cur_rls_cp_trials = rls_cp_trials(ismember(rls_cp_trials,cur_tftrials));
for ll = 1:length(use_lfps)
    [S_drls(ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_drls_trials,:));
    [S_cprls(ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_rls_cp_trials,:));
end

cur_dg_djump_trials = grate_drift_detjump_trials(ismember(grate_drift_detjump_trials,cur_tf_trials));
un_sl_vals = unique(all_trial_sl(cur_dg_djump_trials));
for sl = 1:length(un_sl_vals)
    sl
    cur_trials = cur_dg_djump_trials(all_trial_sl(cur_dg_djump_trials) == un_sl_vals(sl) & used_trial_durs(cur_dg_djump_trials) > movingwin(1));
    for ll = 1:length(use_lfps)
        [S_dg_djump(ll,sl,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
    end
end
%% TRIG AVG WAVELET SPECGRAMS






