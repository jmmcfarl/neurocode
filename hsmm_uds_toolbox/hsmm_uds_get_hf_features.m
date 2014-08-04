function [hf_features,t_axis,Fs] = hsmm_uds_get_hf_features(sig,Fs_orig,Fs_desired,cut_off_freqs,smooth_sigma)

% compute 'High-frequency power' signal feature from input signal sig
%
% Input: 
%       sig: signal to be analyzed (normalized to 0 mean and unit variance)
%       Fs_orig:  sample frequency of sig (Hz)
%       dsf: desired down-sample factor (Default = 1)
%       cut_off_freqs: 2-element vector containing the desired low and
%         high-cutoff frequencies in Hz (Default: [0.05 2]).  Setting either
%         cutoff frequency to Nan will produce either high-pass or low-pass
%         filtering rather than band-pass filtering.  Setting both to nan
%         will result in no-filtering
%       smooth_sigma: Sigma for Gaussian smoothing to compute RMS power (in
%         seconds)
%              
% Output:
%        hf_features: Smoothed estimate of signal power in the defined
%           frequency band
%        t_axis: associated time-axis (in seconds)
%        Fs: sample-frequency of lf_features

if nargin < 2
    error('Need original sample frequency');
end
if nargin < 3
    dsf = 1;
end
if nargin < 4
    cut_off_freqs = [0.05 2];
end

%%
dsf = round(Fs_orig/Fs_desired);
Fs = Fs_orig/dsf;
niqf = Fs_orig/2;
lcf = cut_off_freqs(1); %low cut-off for the low-freq filter
hcf = cut_off_freqs(2); %high cut-off for the low-freq filter

%filter with a symmetric Butterworth filter (2-poles per cut-off)
if isnan(lcf) 
    error('Must input a low-cutoff filter for high-pass filtering.')
end
if ~isnan(hcf)    
    [b,a] = butter(2,[lcf hcf]/niqf);
    obs_hf = filtfilt(b,a,sig);
else
    [b,a] = butter(2,lcf/niqf,'high');
    obs_hf = filtfilt(b,a,sig);
end
    
obs_hf = hsmm_uds_gsmooth(obs_hf.^2,round(smooth_sigma*Fs_orig)); %smoothed RMS power (don't need to take sqrt because of log transformation and renormalization
obs_hf = zscore(obs_hf);
obs_hf = log(obs_hf-min(obs_hf)+1); %log transform the rms power to normalize
obs_hf = downsample(obs_hf,dsf); %normalize and down-sample
hf_features = zscore(obs_hf(:));

t_axis = (1:length(hf_features))/Fs;