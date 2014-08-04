% Compute the shift between the LFP and MP
clear all

if exist('all_eeg_data.mat')
    load all_eeg_data
elseif exist('part1_eeg_data.mat')
    load part1_eeg_data
elseif exist('part2_eeg_data.mat')
    load part2_eeg_data
else
    disp('ERROR DATA DOES NOT EXIST')
    return
end

clear CSC4* 
clear *NumberValidSamples *ChannelNumbers

set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);

load sync_times
if exist('spike_time.mat')
    load spike_time
elseif exist('spike_time_br.mat')
    load spike_time_br
else
    disp('ERROR NO SPIKE TIME DATA')
end

global synct

[len,widwcv] = size(CSC1_Samples);

% Get the data.

wcv = reshape(CSC1_Samples,len*widwcv,1);
wcv = -wcv; % this sign flip ensures that spikes go upwards.
wcv = detrend(wcv); % detrend wcv.
Fs = 32258;

save used_data_full_wcv wcv Fs