function [hf_features,t_axis,Fs] = get_causal_hf_sig(raw_data,raw_Fs,Fs_desired,filter_range)


dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;

lcf = filter_range(1);
hcf = filter_range(2);
[b_hf,a_hf] = butter(2,[lcf/niqf hcf/niqf]);

%compute low-freq observations
obs_hf = filter(b_hf,a_hf,raw_data);
obs_hf = zscore(downsample(obs_hf,dsf));
hf_features = obs_hf(:);

t_axis = (1:length(hf_features))/Fs;