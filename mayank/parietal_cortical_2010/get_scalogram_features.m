function [scal_features,wfreqs,t_axis,Fs] = get_scalogram_features(raw_data,raw_Fs,Fs_desired)


% Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
hf_lcf = 2;
hf_hcf = 100;

min_freq = 10; max_freq = 40; delta_j = 0.05;
k0 = 6; %wavenumber for morlet wavelet
fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));
min_scale = 1/max_freq/fourier_factor;
max_scale = 1/min_freq/fourier_factor;
n_scales = round(log2(max_scale/min_scale)/delta_j);

%compute filter coefficients
[b_hf,a_hf] = butter(2,[hf_lcf/niqf hf_hcf/niqf]);

obs_hf = filtfilt(b_hf,a_hf,raw_data);
dsf1 = 8;
obs_hf = downsample(obs_hf,dsf1);
Fsd = raw_Fs/dsf1;
dsf2 = dsf/dsf1;

[wav_trans,periods,scales,COI] = wavelet(obs_hf,1/Fsd,0,delta_j,min_scale,n_scales);
wfreqs = 1./periods;
inv_scales = 1./scales';
raw_scalogram = abs(wav_trans).^2;
freqs1 = find(wfreqs <= 20);
freqs2 = find(wfreqs > 20);
pow_series1 = mean(raw_scalogram(freqs1,:));
pow_series2 = mean(raw_scalogram(freqs2,:));
pow_series1 = zscore(log(pow_series1 + 1e-5));
pow_series2 = zscore(log(pow_series2 + 1e-5));
scal_features = [downsample(pow_series1',dsf2) downsample(pow_series2',dsf2)];

% n = size(raw_scalogram,2);
% sm_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*raw_scalogram,1,delta_j,scales,0.6);
% scal_features = downsample(sm_scalogram,round(0.5/delta_j));
% wfreqs = downsample(wfreqs,round(0.5/delta_j));
% scal_features = log(scal_features+1);
% scal_features = zscore(scal_features');
% scal_features = downsample(scal_features,dsf2);

t_axis = (1:size(scal_features,1))/Fs;