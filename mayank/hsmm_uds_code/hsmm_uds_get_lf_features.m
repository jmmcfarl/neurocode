function [lf_features,t_axis,Fs] = hsmm_uds_get_lf_features(sig,Fs_orig,dsf,cut_off_freqs)

% compute 'low-frequency amplitude' signal feature from input signal sig
%
% Input: 
%       sig: signal to be analyzed (normalized to 0 mean and unit variance)
%       Fs_orig:  sample frequency of sig (Hz)
%       dsf: desired down-sample factor (Default = 1)
%       cut_off_freqs: 2-element vector containing the desired low and
%       high-cutoff frequencies in Hz (Default: [0.05 2]).  Setting either
%       cutoff frequency to Nan will produce either high-pass or low-pass
%       filtering rather than band-pass filtering.  Setting both to nan
%       will result in no-filtering
%              
% Output:
%        lf_features: filtered and down-sampled signal amplitude
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
Fs = Fs_orig/dsf;
niqf = Fs_orig/2;
lcf = cut_off_freqs(1); %low cut-off for the low-freq filter
hcf = cut_off_freqs(2); %high cut-off for the low-freq filter

%filter with a symmetric Butterworth filter (2-poles per cut-off)
if isnan(lcf) & ~isnan(hcf)
    [b,a] = butter(2,hcf/niqf,'low');
    obs_lf = filtfilt(b,a,sig);
elseif ~isnan(lcf) & isnan(hcf)
    [b,a] = butter(2,lcf/niqf,'high');
    obs_lf = filtfilt(b,a,sig); 
elseif ~isnan(lcf) & ~isnan(hcf)
    [b,a] = butter(2,[lcf hcf]/niqf);
    obs_lf = filtfilt(b,a,sig);
end
    
obs_lf = zscore(downsample(obs_lf,dsf)); %down-sample and normalize
lf_features = obs_lf(:);
t_axis = (1:length(lf_features))/Fs;

