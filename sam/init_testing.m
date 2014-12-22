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