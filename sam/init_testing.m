clear all

data_dir = '/Users/james/Data/sam/12_16_14/';
cd(data_dir);

%%
load spikeSorting.mat
load driftingGrating_visStimInfo.mat

%%
n_units = size(ic,2);
for cc = 1:n_units
   cur_ind_range = ic(3,cc):ic(4,cc);
   spike_times{cc} = t(cur_ind_range);
end

%%
trial_dur = 5;
dt = 0.025;
trial_tbin_edges = (0:dt:trial_dur)*1e3;
trial_tbin_cents = 0.5*trial_tbin_edges(1:end-1) + 0.5*trial_tbin_edges(2:end);

n_trials = length(trialStartT);
tbt_spikes = cell(n_trials,n_units);
psths = zeros(n_units,length(trial_tbin_cents));
for ii = 1:n_trials
    for cc = 1:n_units
       cur_trial_spikes = spike_times{cc}(spike_times{cc} >= trialStartTMS(ii) & ...
           spike_times{cc} <= trialEndTMS(ii));
       tbt_spikes{ii,cc} = cur_trial_spikes - trialStartTMS(ii);
       
       cur_hist = histc(tbt_spikes{ii,cc},trial_tbin_edges);
       psths(cc,:) = psths(cc,:) + cur_hist(1:end-1);
    end
end
psths = psths/n_trials;

trial_spike_counts = cellfun(@(x) length(x),tbt_spikes);
ov_trial_avgs = mean(trial_spike_counts);

%%
sm_sig = round(0.1/dt);
norm_psths = bsxfun(@rdivide,psths/dt,ov_trial_avgs');
for ii = 1:n_units
   norm_psths(ii,:) = jmm_smooth_1d_cor(norm_psths(ii,:),sm_sig); 
end

uset = find(ov_trial_avgs >= 0.01);
f1 = figure();
shadedErrorBar(trial_tbin_cents/1e3,mean(norm_psths(uset,:)),std(norm_psths(uset,:))/sqrt(length(uset)));
xlabel('Time (s)');
ylabel('Rate (Hz)');

%%
un_oris = unique(angleOrder);
ori_rates = nan(length(un_oris),n_units);
for ii = 1:length(un_oris)
   cur_trials = find(angleOrder == un_oris(ii));
   ori_rates(ii,:) = mean(trial_spike_counts(cur_trials,:));
end
ori_nrates = bsxfun(@rdivide,ori_rates,ov_trial_avgs);

%%
un_sfs = unique(spatialFreq);
sf_rates = nan(length(un_sfs),n_units);
for ii = 1:length(un_sfs)
   cur_trials = find(spatialFreq == un_sfs(ii));
   sf_rates(ii,:) = mean(trial_spike_counts(cur_trials,:));
end
sf_nrates = bsxfun(@rdivide,sf_rates,ov_trial_avgs);

%%
un_tfs = unique(tempFreq);
tf_rates = nan(length(un_tfs),n_units);
for ii = 1:length(un_tfs)
   cur_trials = find(tempFreq == un_tfs(ii));
   tf_rates(ii,:) = mean(trial_spike_counts(cur_trials,:));
end
tf_nrates = bsxfun(@rdivide,tf_rates,ov_trial_avgs);

%% LOAD LFP DATA 
use_chs = 1:4:120;
lfp_dsf = 6;
load('headerInfo')
LFP_Fs = hdr.filteredFS;
LFP_Fsd = LFP_Fs/lfp_dsf;
[bb,aa] = butter(2,[0.1 0.8*LFP_Fsd]/LFP_Fs);

lfp_files = dir('chunk*.mat');
n_lfp_files = length(lfp_files);
n_use_files = 200;
all_lfps = [];
for ii = 1:n_use_files
    cur_fname = sprintf('chunk%d.mat',ii);
    fprintf('Loading %s\n',cur_fname);
    load(cur_fname);
    cur_lfps = double(lfp(use_chs,:)');
    cur_lfps = filtfilt(bb,aa,cur_lfps);
    all_lfps = cat(1,all_lfps,downsample(cur_lfps,lfp_dsf));
end
all_lfps = zscore(all_lfps);

lfp_taxis = (1:size(all_lfps,1))/LFP_Fsd;
%%

clear sig_var
tstart_inds = round(interp1(lfp_taxis,1:length(lfp_taxis),trialStartTMS/1e3));
tstart_inds(tstart_inds > length(lfp_taxis) | isnan(tstart_inds)) = [];

backlag = round(0*LFP_Fsd);
forlag = round(3*LFP_Fsd);
[temp,templags,~,~,tempmat] = get_event_trig_avg_v3(all_lfps,tstart_inds,backlag,forlag,2);
tempmat = tempmat';
tempmat = reshape(tempmat,[forlag+backlag+1, size(all_lfps,2), size(tempmat,2)]);

%%
tempFreq_bin_edges = prctile(tempFreq,[0:8:100]);
[~,tfreq_bin] = histc(tempFreq,tempFreq_bin_edges);
%%
params.Fs = LFP_Fsd;
params.tapers = [5 9];
params.trialave = 0;
for cc = 1:size(tempmat,2)
    cc
cur_data = squeeze(tempmat(:,cc,:));
[S,f] = mtspectrumc(cur_data,params);

%%
tot_trials = size(S,2);
un_oris = unique(angleOrder);
ori_Savg = nan(length(un_oris),length(f));
for oo = 1:length(un_oris)
    tset = find(angleOrder == un_oris(oo));
    tset(tset > tot_trials) = [];
    ori_Savg(oo,:) = nanmean(log10(S(:,tset)),2);
end
ori_Savg = bsxfun(@minus,ori_Savg,nanmean(log10(S),2)');
ori_Savg = bsxfun(@rdivide,ori_Savg,nanstd(log10(S),[],2)');
% f1 = figure();
% imagesc(f,un_oris,ori_Savg);
%%
tot_trials = size(S,2);
un_sfs = unique(spatialFreq);
sf_Savg = nan(length(un_sfs),length(f));
for oo = 1:length(un_sfs)
    tset = find(spatialFreq == un_sfs(oo));
    tset(tset > tot_trials) = [];
    sf_Savg(oo,:) = nanmean(log10(S(:,tset)),2);
end
sf_Savg = bsxfun(@minus,sf_Savg,nanmean(log10(S),2)');
sf_Savg = bsxfun(@rdivide,sf_Savg,nanstd(log10(S),[],2)');
% f1 = figure();
% imagesc(f,un_sfs,ori_Savg);

%%
tot_trials = size(S,2);
un_tfs = unique(tfreq_bin);
tf_Savg = nan(length(un_tfs),length(f));
for oo = 1:length(un_tfs)
    tset = find(tfreq_bin == un_tfs(oo));
    tset(tset > tot_trials) = [];
    tf_Savg(oo,:) = nanmean(log10(S(:,tset)),2);
end
tf_Savg = bsxfun(@minus,tf_Savg,nanmean(log10(S),2)');
tf_Savg = bsxfun(@rdivide,tf_Savg,nanstd(log10(S),[],2)');

%%
all_ori_Savg(cc,:,:) = ori_Savg;
all_sf_Savg(cc,:,:) = sf_Savg;
all_tf_Savg(cc,:,:) = tf_Savg;
end

%%
for cc = 1:size(tempmat,2)
%    imagesc(f,un_oris,squeeze(all_ori_Savg(cc,:,:)))
   imagesc(f,un_sfs,squeeze(all_sf_Savg(cc,:,:)))
%    imagesc(f,un_tfs,squeeze(all_tf_Savg(cc,:,:)))
   pause
   clf
   
    
end

%%
pow_range = find(f > 60);
ov_ori_pow = nanmean(all_ori_Savg(:,:,pow_range),3);
ov_sf_pow = nanmean(all_sf_Savg(:,:,pow_range),3);
ov_tf_pow = nanmean(all_tf_Savg(:,:,pow_range),3);
f1 = figure();
subplot(3,1,1);
shadedErrorBar(un_oris,nanmean(ov_ori_pow),nanstd(ov_ori_pow)/sqrt(size(ov_ori_pow,1)));
xlabel('Orientation');
ylabel('Avg pow (z)');
axis tight
subplot(3,1,2);
shadedErrorBar(un_sfs,nanmean(ov_sf_pow),nanstd(ov_sf_pow)/sqrt(size(ov_sf_pow,1)));
xlabel('spatial freq');
ylabel('Avg pow (z)');
axis tight
subplot(3,1,3);
shadedErrorBar(un_tfs,nanmean(ov_tf_pow),nanstd(ov_tf_pow)/sqrt(size(ov_tf_pow,1)));
xlabel('temporal freq');
ylabel('Avg pow (z)');
axis tight
%%
MEA_nums = hdr.MEAlayoutNums;
[mea_vals,mea_ord] = sort(MEA_nums(:));
[XX,YY] = meshgrid(1:12);
xvals = XX(mea_ord);
yvals = YY(mea_ord);
xvals(isnan(mea_vals)) = [];
yvals(isnan(mea_vals)) = [];

[~,xsort] = sort(xvals);
[~,ysort] = sort(yvals);
[~,dsort] = sort(xvals.^2 + yvals.^2);

%%
% for ii = 1:length(templags)
ii = 168
    templags(ii)/LFP_Fsd
Vq = griddata(xvals,yvals,temp(ii,:)',XX,YY);
imagescnan(Vq); caxis([-0.35 0.35])
xlabel('X'); ylabel('Y'); 
% pause(0.2)
% clf
% end

%%
params.Fs = LFP_Fsd;
params.tapers = [2 3];
params.trialave = 1;
movingwin = [0.25 0.1];
win = [-0.5 5];
E = trialStartTMS/1e3;
E(E > lfp_taxis(end)) = [];
E(end) = [];
for cc = 1:size(all_lfps,2)
[S,t,f]=mtspecgramtrigc(all_lfps(:,cc),E,win,movingwin,params);
all_Sgram(cc,:,:) = log10(S);
end
