clear all
clc

%% Load Data
cd ~/Data/bruce/7_15_12/G034/

load ./Expt3_fixbased_data.mat trial* full_t
trial_start_times = full_t(trial_start_inds);
trial_stop_times = full_t(trial_stop_inds);
trial_durs = trial_stop_times-trial_start_times;
%%

load ./expt3_parsed_g034
mic_sacs = find(all_sac_amps < 1);
mac_sacs = find(all_sac_amps > 1);
all_sac_times = all_sac_times(mic_sacs);
% all_sac_times = all_sac_times(mac_sacs);
%%
eye_dt = median(diff(all_t));
full_dt = median(diff(full_t));
avg_sac_rate = length(all_sac_times)/length(all_t)/eye_dt;

n_bins = 30;
rel_t_bin_edges = linspace(0,0.5,n_bins+1);
rel_t_bin_cents = 0.5*rel_t_bin_edges(1:end-1) + 0.5*rel_t_bin_edges(2:end);
rel_bin_dt = median(diff(rel_t_bin_cents));

rel_sac_hist = zeros(n_bins,1);
rel_occ = zeros(n_bins,1);
for i = 1:length(trial_start_inds)
    cur_sac_set = find(all_sac_times >= trial_start_times(i) & all_sac_times < trial_stop_times(i));
    if ~isempty(cur_sac_set)
        rel_sac_times = all_sac_times(cur_sac_set) - trial_start_times(i);
        cur_hist = histc(rel_sac_times,rel_t_bin_edges);
        cur_hist = cur_hist(1:end-1);
        rel_sac_hist = rel_sac_hist + cur_hist(:);
    end
    cur_dur = trial_durs(i);
    last_occ_bin = find(rel_t_bin_edges < cur_dur,1,'last');
    last_occ_bin = min(last_occ_bin,n_bins);
    rel_occ(1:last_occ_bin) = rel_occ(1:last_occ_bin) + 1;
end
rel_sac_rate = rel_sac_hist./rel_occ/rel_bin_dt;