clear all

% depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_1s_POS_CI';
depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_1s_POS_CI';
depol_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_1s_POS';
% depol_array{4} = 'C:\WC_Germany\april_09_data\2009-04-21_CWC_LFP\2009-4-21_12s_1s_POS_CI';
% depol_array{5} = 'C:\WC_Germany\april_09_data\2009-05-16_CWC_LFP\2009-5-16_1s_POS_CI';


ci_time = 1;
wait_time = 2;
first_pulse = 1;

Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[100/niqf 5e3/niqf]);
spkthresh = 5;

mp_min = -100;
mp_max = 10;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

for d = 1:2
d
    %% get depolarizing CI
    load(depol_array{d})

    time = (1:length(data))/Fs;
%     data = downsample(data,dsf);
%     time = downsample(time,dsf);
diff_data = [0;diff(data)];
alen = length(data);
filt_data = zscore(filtfilt(b,a,data));

spkid = find( (filt_data(1:alen-1) < spkthresh) & (filt_data(2:alen) >= spkthresh) );
% bad_spks = find(data(spkid) > data(spkid+25));
% spkid(bad_spks) = [];
bad_spks = find(diff(spkid) < 0.001*Fs);
spkid(bad_spks) = [];

    num_cis = floor((max(time)-2*first_pulse)/(ci_time+wait_time));
    stim_times = [];
    artifact_times = [];
    for i = 1:num_cis

        begpt = (i-1)*3*Fs+1+Fs;
        endpt = begpt+ci_time*Fs;
        stim_times = [stim_times begpt:endpt];
        artifact_times = [artifact_times (begpt-25):begpt endpt:(endpt+25)];
    end
    spkid(ismember(spkid,artifact_times)) = [];
    
            spkRmWidth = 300; %width of spike to remove
data_minus_spike{d} = data;
for i = 1:length(spkid)-1
    begPt = data(spkid(i)-20);
    endPt = data(spkid(i)+spkRmWidth-1);
    interpSlope = (endPt - begPt)/spkRmWidth;
    data_minus_spike{d}(spkid(i)-20:spkid(i)+spkRmWidth-21) = [1:spkRmWidth]*interpSlope+begPt;
end

% plot(data); hold on; plot(spkid,data(spkid),'r.')
% pause
% close all
    stim_spikes = ismember(spkid,stim_times);
    num_stim_spikes = sum(stim_spikes);
    num_non_stim_spikes = length(spkid)-num_stim_spikes;
    stim_time = num_cis*ci_time;
    non_stim_time = max(time) - stim_time;
    stim_rate(d) = num_stim_spikes/stim_time;
    non_stim_rate(d) = num_non_stim_spikes/non_stim_time;
    
end

stim_rate
non_stim_rate
% cd C:\WC_Germany\current_injection\persistent_ci
% save april_09_1s_ci_data tvec mean_* avg_amp* norm_mean* prior*


