function [] = hip_wc_lfp_spk_shift_combined()

% Compute the shift between the LFP and MP
% clear all

if exist('./all_eeg_data2.mat','file')
    load ./all_eeg_data2
elseif exist('./all_eeg_data.mat','file')
    load ./all_eeg_data
else
    disp('ERROR DATA DOES NOT EXIST')
    return
end

clear *NumberValidSamples *ChannelNumbers

load ./sync_times
if exist('spike_time_jmm.mat')
    load spike_time_jmm
else
    disp('ERROR NO SPIKE TIME DATA')
end

% load ./used_data ep
% if exist('ep','var')
%     if length(synct1id) ~= ep
%         error('EP!')
%     end
% end

[len,widwcv] = size(CSC1_Samples);
[widlfp] = length(CSC8_Samples(1,:));
[widlf3] = length(CSC3_Samples(1,:));
[widlf2] = length(CSC2_Samples(1,:));

% Get the data.
% wcv = reshape(CSC1_Samples,len*widwcv,1)/gains(1)*gain_to_volts;
wcv = reshape(CSC1_Samples,len*widwcv,1);
wcv = wcv(synct1id);
wcv = -wcv; % this sign flip ensures that spikes go upwards.
wcv = detrend(wcv); % detrend wcv.
clear CSC1_Samples;

% lf2 = reshape(CSC2_Samples,len*widlf2,1)/gains(2)*gain_to_volts;
% lf3 = reshape(CSC3_Samples,len*widlf3,1)/gains(3)*gain_to_volts;
% lf4 = reshape(CSC4_Samples,len*widlfp,1)/gains(4)*gain_to_volts;
% lf5 = reshape(CSC5_Samples,len*widlfp,1)/gains(5)*gain_to_volts;
% lf6 = reshape(CSC6_Samples,len*widlfp,1)/gains(6)*gain_to_volts;
% lf7 = reshape(CSC7_Samples,len*widlfp,1)/gains(7)*gain_to_volts;
% lf8 = reshape(CSC8_Samples,len*widlfp,1)/gains(8)*gain_to_volts;
lf2 = reshape(CSC2_Samples,len*widlf2,1);
lf3 = reshape(CSC3_Samples,len*widlf3,1);
lf4 = reshape(CSC4_Samples,len*widlfp,1);
lf5 = reshape(CSC5_Samples,len*widlfp,1);
lf6 = reshape(CSC6_Samples,len*widlfp,1);
lf7 = reshape(CSC7_Samples,len*widlfp,1);
lf8 = reshape(CSC8_Samples,len*widlfp,1);

wcv_Fs = median(CSC1_SampleFrequencies);
lf3_Fs = median(CSC3_SampleFrequencies);
lf2_Fs = median(CSC3_SampleFrequencies);

if lf2_Fs == wcv_Fs
    lf2 = lf2(synct1id);
else
    lf2 = lf2(synct2id);
end
if lf3_Fs == wcv_Fs
    lf3 = lf3(synct1id);
else
    lf3 = lf3(synct2id);
end
lf4 = lf4(synct2id);
lf5 = lf5(synct2id);
lf6 = lf6(synct2id);
lf7 = lf7(synct2id);
lf8 = lf8(synct2id);

% clear *Samples *TimeStamps *SampleFrequencies
% clear CSC2_Samples CSC5_Samples CSC8_Samples;
clear *_Samples

%% This step should be unnecesary because of hardward filters, also because
%% of later high-pass filtering steps
% lf2 = detrend(lf2);
% lf3 = detrend(lf3);
% lf4 = detrend(lf4);
% lf5 = detrend(lf5);
% lf6 = detrend(lf6);
% lf7 = detrend(lf7);
% lf8 = detrend(lf8);

synct = synct*10^-6;
dt = median(diff(synct));


wcv_minus_spike = wcv;


%************* THIS IS MAYANK'S SPIKE SUBTRACTION ALGO
% for k = 1:5
%     wcv_minus_spike(spkid+k-1) = mean(wcv) + 2*std(wcv);
% end
%***************

%spike subtraction jmm
spkRmWidth = 6; %width of spike to remove
for i = 1:length(spkid)-1
    begPt = wcv(spkid(i)-1);
    endPt = wcv(spkid(i)+spkRmWidth+1);
    interpSlope = (endPt - begPt)/spkRmWidth;
    wcv_minus_spike(spkid(i):spkid(i)+spkRmWidth-1) = [1:spkRmWidth]*interpSlope+begPt;
end


save used_data lf8 lf3 lf4 lf5 lf2 lf6 lf7 wcv wcv_minus_spike

% dataViewer(lf8,wcv_minus_spike,synct,dt)


