function [lf_features,t_axis,Fs] = get_lf_features_acaus(raw_data,raw_Fs,Fs_desired,cut_off_freqs)


% Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
lcf = cut_off_freqs(1); %low cut-off for the low-freq filter
hcf = cut_off_freqs(2); %high cut-off for the low-freq filter

[b,a] = butter(2,[lcf/niqf hcf/niqf]);

%compute low-freq observations
obs_lf = filter(b,a,(flipud(raw_data(:))));
obs_lf = flipud(obs_lf);
lf_features = zscore(downsample(obs_lf,dsf));
t_axis = (1:length(lf_features))/Fs;