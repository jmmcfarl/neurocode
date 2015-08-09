function [hf_features,t_axis,Fs] = get_hf_lpower(raw_data,raw_Fs,Fs_desired,filter_range,hf_smooth)


dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;

lcf = filter_range(1);
hcf = filter_range(2);
[b_hf,a_hf] = butter(2,[lcf/niqf hcf/niqf]);

obs_hf = filtfilt(b_hf,a_hf,raw_data); %high-pass filter
eps = var(obs_hf)/100; %compute min amp was var(.)/10 
obs_hf = obs_hf.^2; %square signal
if hf_smooth > 0
    obs_hf = jmm_smooth_1d_cor(obs_hf,round(hf_smooth*raw_Fs)); %smooth with gaussian
end
obs_hf(obs_hf < eps) = eps; %impose minimum to avoid infinite values if there was any saturation-based flatness
obs_hf = log10(obs_hf); %log power
obs_hf = downsample(obs_hf,dsf);

hf_features = obs_hf(:);
t_axis = (1:length(hf_features))/Fs;