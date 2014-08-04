clear all
cd C:\WC_Germany\current_injection\persistent_ci

mp_min = -1.0;
mp_max = 0.1;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

Fs = 500;
niqf = Fs/2;
[b,a] = butter(2,100/niqf,'high');
spkthresh = 1.5;


%% FOR OLD DATA
old_depol_array{1} = 'C:\WC_Germany\longdepol\2007-05-31_CWC_LFP\2007-05-31_CWC_LFP_Hekadata\sweep_data';
old_depol_array{2} = 'C:\WC_Germany\longdepol\2007-06-04_CWC_LFP_B\2007-06-04_CWC_LFP_B_Hekadata\sweep_data';
sweep_beg = Fs;
sweep_end = round(13*Fs);
for d = 1:2
       spkid = [];
    stim_times = [];
    artifact_times = [];
    load(old_depol_array{d})
   for s = 1:size(sweep_data,1)
        cur_data = sweep_data(s,:);
        diff_data = [0 diff(cur_data)];
        alen = length(cur_data);
        filt_cur_data = zscore(filtfilt(b,a,diff_data));
        cur_spkid = find((filt_cur_data(1:alen-1) < spkthresh) & (filt_cur_data(2:alen) >= spkthresh) );
%         bad_spks = find(cur_data(cur_spkid) > cur_data(cur_spkid+30));
%         for t = 1:length(cur_spkid)
%             if min(filt_cur_data(cur_spkid(t)-10:cur_spkid(t))) < -10
%                 bad_spks = [bad_spks t];
%             end
%         end
%         cur_spkid(unique(bad_spks)) = [];
        spkid = [spkid ((s-1)*size(sweep_data,2)+cur_spkid)];
        stim_times = [stim_times ((s-1)*size(sweep_data,2)+sweep_beg):((s-1)*size(sweep_data,2)+sweep_end)];
        artifact_times = [artifact_times ((s-1)*size(sweep_data,2)+sweep_beg-4):((s-1)*size(sweep_data,2)+sweep_beg+4) ...
            ((s-1)*size(sweep_data,2)+sweep_end-4):((s-1)*size(sweep_data,2)+sweep_end+4)];
   end
   spkid(ismember(spkid,artifact_times)) = [];
   
data = reshape(sweep_data',numel(sweep_data),1);

spkRmWidth = 5; %width of spike to remove
data_minus_spike{d} = data;
for i = 1:length(spkid)-1
    begPt = data(spkid(i)-1);
    endPt = data(spkid(i)+spkRmWidth-1);
    interpSlope = (endPt - begPt)/spkRmWidth;
    data_minus_spike{d}(spkid(i):spkid(i)+spkRmWidth-1) = [1:spkRmWidth]*interpSlope+begPt;
end

data_minus_spike{d} = reshape(data_minus_spike{d},size(sweep_data'));
data_minus_spike{d} = data_minus_spike{d}';

% plot(data), hold on, plot(spkid,data(spkid),'r.')
% pause
close all

    stim_spikes = ismember(spkid,stim_times);
    num_stim_spikes = sum(stim_spikes);
    num_non_stim_spikes = length(spkid)-num_stim_spikes;
    stim_time = size(sweep_data,1)*12;
    non_stim_time = size(sweep_data,1)*4;
    stim_rate(d) = num_stim_spikes/stim_time
    non_stim_rate(d) = num_non_stim_spikes/non_stim_time

end


%% FOR NEW DATA
clear all
Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[100/niqf 5e3/niqf]);
spkthresh = 3;

stim_dur = 12;
pause_dur = 18;
first_stim = 9;


% depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_12s_CI';
depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_12s_CI';
depol_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_12s_CI';
depol_array{3} = 'C:\WC_Germany\april_09_data\2009-04-13_B\2009-04-13_CWC_LFP_B_12s_CI_smaller';

for d = 1:3
    
load(depol_array{d})
    time = (1:length(data))/Fs;
diff_data = [0;diff(data)];
alen = length(data);
filt_data = zscore(filtfilt(b,a,data));

spkid = find( (filt_data(1:alen-1) < spkthresh) & (filt_data(2:alen) >= spkthresh) );
% bad_spks = find(data(spkid) > data(spkid+30));
% spkid(bad_spks) = [];
bad_spks = find(diff(spkid) < 0.001*Fs);
spkid(bad_spks) = [];

    num_cis = floor((max(time)-2*first_stim)/(stim_dur+pause_dur))+1;
stim_times = [];
artifact_times = [];
    for i = 1:num_cis

        begpt = (i-1)*30*Fs+1+round(9*Fs);
        endpt = begpt+12*Fs;
        if endpt > length(data)
            endpt = length(data);
            stim_times = [stim_times begpt:endpt];

        else
            stim_times = [stim_times begpt:endpt];
        end
        
        artifact_times = [artifact_times (begpt-40):(begpt+40) (endpt-40):(endpt+40)];
        
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

    
%     plot(data), hold on, plot(spkid,data(spkid),'r.')
%     pause, close all
    
       stim_spikes = ismember(spkid,stim_times);
    num_stim_spikes = sum(stim_spikes);
    num_non_stim_spikes = length(spkid)-num_stim_spikes;
    stim_time = num_cis*12;
    non_stim_time = num_cis*18;
    stim_rate(d+2) = num_stim_spikes/stim_time;
    non_stim_rate(d+2) = num_non_stim_spikes/non_stim_time;
 
end

