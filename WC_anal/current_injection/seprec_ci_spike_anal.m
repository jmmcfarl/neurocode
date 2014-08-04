clear all

depol_dir{1} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-05-24_CWC_LFP_B';
depol_dir{2} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-05-31_CWC_LFP';
depol_dir{3} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-06-04_CWC_LFP_B';
depol_dir{4} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-03_CWC_LFP_A';
depol_dir{5} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-03_CWC_LFP_B';
depol_dir{6} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-04_CWC_LFP_B';
depol_dir{7} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-04_CWC_LFP_D';
depol_dir{8} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-05_CWC_LFP';



mp_min = -1.0;
mp_max = 0.1;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

Fs = 2e4;
% dsf = 20;
% Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[100/niqf 5e3/niqf]);
spkthresh = 5;

wait_samps = round(0.5*Fs);
ci_dur = Fs;

for d = 1:8
    
    load(depol_dir{d})
    fdata = filtfilt(b,a,[0;diff(data)]);
    f_std = std(fdata);
    send_times = find(time==1.9999);
    sbeg_times = find(time==0);
if length(send_times) < length(sbeg_times)
    sbeg_times(end) = [];
end

    num_sweeps = length(sbeg_times);
    stim_times = [];
    artifact_times = [];
    spkid = [];
    if length(sbeg_times) >= 120
        end_sweep = 120;
    else
        end_sweep = 119;
    end
    for i = 61:119
        cur_data = data(sbeg_times(i):send_times(i));
        diff_data = [0;diff(cur_data)];
        alen = length(cur_data);
        filt_cur_data = (filtfilt(b,a,diff_data))/f_std;
        cur_spkid = find((filt_cur_data(1:alen-1) < spkthresh) & (filt_cur_data(2:alen) >= spkthresh) );
%         bad_spks = find(cur_data(cur_spkid) > cur_data(cur_spkid+25));
%         for t = 1:length(cur_spkid)
%             if min(filt_cur_data(cur_spkid(t)-10:cur_spkid(t))) < -20
%                 bad_spks = [bad_spks;t];
%             end
%         end
%         cur_spkid(unique(bad_spks)) = [];
        
        spkid = [spkid;(sbeg_times(i)+cur_spkid)];
        artifact_times = [artifact_times (sbeg_times(i)+wait_samps-25):(sbeg_times(i)+wait_samps+25) ... 
            (sbeg_times(i)+wait_samps+ci_dur):(sbeg_times(i)+wait_samps+ci_dur+25)];
        stim_times = [stim_times (sbeg_times(i)+wait_samps):(sbeg_times(i)+wait_samps+ci_dur)];
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
%     pause
%     close all
    
        stim_spikes = ismember(spkid,stim_times);
    num_stim_spikes = sum(stim_spikes);
    num_non_stim_spikes = length(spkid)-num_stim_spikes;
    stim_time = 60;
    non_stim_time = 60;
    stim_rate(d) = num_stim_spikes/stim_time
    non_stim_rate(d) = num_non_stim_spikes/non_stim_time
    
d
end
% 
% cd C:\WC_Germany\current_injection\persistent_ci
% save seprec_1s_ci_data tvec mean_* avg_amp* norm_mean* prior*


