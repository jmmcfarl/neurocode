% close all
clear all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G096';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

for ii = 1:27
    expt_name{ii} = Expts{ii}.Header.expname;
end
%%
dsf = 1; %lfps originally sampled at 400Hz
use_lfps = 1:4:96;
lcf = 0.5;

min_trial_dur = 1;
cur_expt_set = [6:19];

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
all_trial_sf = [];
all_trial_tf = [];
all_trial_or = [];
all_trial_sz = [];
all_trial_jv = [];
all_trial_nsf = [];

for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    trial_result = [Expts{cur_expt}.Trials(:).Result];
    if isfield(Expts{cur_expt}.Trials,'sf')
        trial_sf = [Expts{cur_expt}.Trials(:).sf];
    else
        trial_sf = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.sf;
    end
    if isfield(Expts{cur_expt}.Trials,'tf')
        trial_tf = [Expts{cur_expt}.Trials(:).tf];
    else
        trial_tf = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.tf;
    end
    if isfield(Expts{cur_expt}.Trials,'or')
        trial_or = [Expts{cur_expt}.Trials(:).or];
    else
        trial_or = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.or;
    end
    trial_sz = [Expts{cur_expt}.Trials(:).sz];
    if isfield(Expts{cur_expt}.Trials,'jv')
        trial_jv = [Expts{cur_expt}.Trials(:).jv];
    else
        trial_jv = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.jv;
    end
    if isfield(Expts{cur_expt}.Trials,'nsf')
        trial_nsf = reshape([Expts{cur_expt}.Trials(:).nsf],10,length(trial_start_times));
    else
        trial_nsf = [];
    end
    
%     [un_ids,id_inds] = unique(trial_ids);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times');
    all_trial_durs = cat(1,all_trial_durs,trial_durs');
    all_trial_result = cat(1,all_trial_result,trial_result');
    all_trial_exptvec = cat(1,all_trial_exptvec,ee*ones(length(trial_ids),1));
    
    all_trial_sf = cat(1,all_trial_sf,trial_sf');
    all_trial_tf = cat(1,all_trial_tf,trial_tf');
    all_trial_or = cat(1,all_trial_or,trial_or');
    all_trial_sz = cat(1,all_trial_sz,trial_sz');
    all_trial_jv = cat(1,all_trial_jv,trial_jv');
    all_trial_nsf = cat(1,all_trial_nsf,trial_nsf');
    
    n_trials = length(trial_ids);

    
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
start_buffer = round(Fsd*0.5);
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));

% stims = [all_trial_sf(:) all_trial_tf(:)];
stims = [all_trial_nsf];
[C,IA,IC] = unique(stims,'rows');
n_types = size(C,1);
%%
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));
params.Fs = Fsd;
params.tapers = [3 5];
params.fpass = [0 125];
movingwin = [1.5 1.5];
sMarkers = [all_trial_start_inds(:)+start_buffer all_trial_end_inds(:)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];

clear S
for ii = 1:n_types
    fprintf('Condition %d of %d\n',ii,n_types);
    cur_trials = find(IC == ii & used_trial_durs > movingwin(1));
    
    for ll = 1:length(use_lfps)
        [S(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
    end
        
end

%%
base_C = 5;
logS = log10(S);
avg_pow = squeeze(mean(logS(base_C,:,:),2))';

beta = polyfit(f,avg_pow,1);
pred_pow = polyval(beta,f);

% norm_S = bsxfun(@minus,logS,reshape(pred_pow,[1 1 length(f)]));
norm_S = bsxfun(@minus,logS,reshape(avg_pow,[1 1 length(f)]));
%%
close all
for ii = 1:n_types
    plot(f,squeeze(mean(log10(S(ii,:,:)),2)));
    pause
    clf
end

%%
cmap =[1 0 0; 0 1 0; 0 0 1];
for ii = 2:4
    
    plot(f,squeeze(mean(log10(S(ii,:,:)),2)),'color',cmap(ii-1,:),'linewidth',2);
hold on
end
%%
close all
cmap =[1 0 0; 0 1 0; 0 0 1];
unique_sfs = unique(all_trial_sf);
unique_tfs = unique(all_trial_tf);
for jj = 1:length(unique_tfs);
    use_tf = unique_tfs(jj);
    subplot(2,3,jj)
    for ii = 1:length(unique_sfs)
    cur_C = find(C(:,1) == unique_sfs(ii) & C(:,2) == use_tf);
%     plot(f,squeeze(mean(norm_S(cur_C,:,:),2)),'color',cmap(ii,:),'linewidth',2);
    plot(f,squeeze(mean(logS(cur_C,:,:),2)),'color',cmap(ii,:),'linewidth',2);
    hold on
    end
title(sprintf('TF = %d',use_tf));
legend('SF = 1','SF = 2','SF = 4');
xlabel('Frequency (Hz)','fontsize',12);
ylabel('Log power','fontsize',12);
ylim([-13.5 -9.5]);
xlim([0 120])
end