function [scal_features,wfreqs,t_axis,Fs] = get_full_scalogram_features(raw_data,raw_Fs,Fs_desired)


% Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
hf_lcf = 2;
hf_hcf = 100;

min_freq = 10; max_freq = 80; 
delta_j = 0.1;
%delta_j = 0.05;
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
n= size(raw_scalogram,2);
sm_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*raw_scalogram,1/Fsd,delta_j,scales,0.6);
w_dsf = round(0.2/delta_j);
sm_scalogram = downsample(sm_scalogram,w_dsf);
% sm_scalogram = downsample(raw_scalogram,w_dsf);
wfreqs = downsample(wfreqs,w_dsf);

if dsf2 > 1
    scal_features = [downsample(sm_scalogram',dsf2)];
else
    scal_features = sm_scalogram';
end

% n = size(raw_scalogram,2);
% sm_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*raw_scalogram,1,delta_j,scales,0.6);
% scal_features = downsample(sm_scalogram,round(0.5/delta_j));
% wfreqs = downsample(wfreqs,round(0.5/delta_j));
% scal_features = log(scal_features+1);
% scal_features = zscore(scal_features');
% scal_features = downsample(scal_features,dsf2);

t_axis = (1:size(scal_features,1))/Fs;