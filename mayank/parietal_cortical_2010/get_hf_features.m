function [hf_features,t_axis,Fs] = get_hf_features(raw_data,raw_Fs,Fs_desired,filter_range,hf_smooth)


dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;

lcf = filter_range(1);
hcf = filter_range(2);
[b_hf,a_hf] = butter(2,[lcf/niqf hcf/niqf]);

%compute low-freq observations
obs_hf = filtfilt(b_hf,a_hf,raw_data);
obs_hf = sqrt(jmm_smooth_1d_cor(obs_hf.^2,round(hf_smooth*raw_Fs))); %smoothed RMS power
% obs_hf = sqrt(boxcar_smooth(obs_hf.^2,round(hf_smooth*raw_Fs))); %smoothed RMS power
obs_hf = zscore(obs_hf);
obs_hf = log(obs_hf-min(obs_hf)+1); %log transform the rms power to normalize
obs_hf = downsample(obs_hf,dsf);
obs_hf = zscore(obs_hf);
hf_features = obs_hf(:);

t_axis = (1:length(hf_features))/Fs;