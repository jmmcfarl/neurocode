function [np_state_seq,t_axis,Fs] = get_np_state_sequence(raw_data,raw_Fs,Fs_desired)

dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;

niqf = raw_Fs/2;
lcf = 0.05;
hcf = 2;
[b,a] = butter(2,[lcf/niqf hcf/niqf]);

data = filtfilt(b,a,raw_data);
data = zscore(downsample(data,dsf));

mix=gmm(1,2,'full');
gmm_options(3) = 1e-15; %tolerance
gmm_options(5) = 1;%reset cov if singular values
gmm_options(14) = 200; %max iterations
[mix, netlaboptions] = gmmem(mix, data, gmm_options);

state_means = mix.centres;
[state_means,state_order] = sort(state_means);
state_vars = squeeze(mix.covars);
state_priors = mix.priors;
state_vars = state_vars(state_order);
state_priors = state_priors(state_order);

thresh_range = linspace(-3,3,1000);

state1_lik = state_priors(1)*normpdf(thresh_range,state_means(1),state_vars(1));
state2_lik = state_priors(2)*normpdf(thresh_range,state_means(2),state_vars(2));

[dummy,state1_peakloc] = findpeaks(state1_lik);

threshold_loc = state1_peakloc+find(state1_lik(state1_peakloc(1):end) < state2_lik(state1_peakloc(1):end),1,'first');
threshold = thresh_range(threshold_loc);

