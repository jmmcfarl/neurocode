clear all

depol_dir{1} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-05-24_CWC_LFP_B';
depol_dir{2} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-05-31_CWC_LFP';
depol_dir{3} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-06-04_CWC_LFP_B';
depol_dir{4} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-03_CWC_LFP_A';
depol_dir{5} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-03_CWC_LFP_B';
depol_dir{6} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-04_CWC_LFP_B';
depol_dir{7} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-04_CWC_LFP_D';
depol_dir{8} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-05_CWC_LFP';

prior_dir{1} = 'C:\WC_Germany\EC_MPasci\A2007_05_24_CWC_LFP_B';
prior_dir{2} = 'C:\WC_Germany\EC_MPasci\A2007_05_31_CWC_LFP';
prior_dir{3} = 'C:\WC_Germany\EC_MPasci\A2007_06_04_CWC_LFP_B';
prior_dir{4} = 'C:\WC_Germany\EC_MPasci\A2007_07_03_CWC_LFP_A';
prior_dir{5} = 'C:\WC_Germany\EC_MPasci\A2007_07_03_CWC_LFP_B';
prior_dir{6} = 'C:\WC_Germany\EC_MPasci\A2007_07_04_CWC_LFP_B';
prior_dir{7} = 'C:\WC_Germany\EC_MPasci\A2007_07_04_CWC_LFP_D';
prior_dir{8} = 'C:\WC_Germany\EC_MPasci\A2007_07_05_CWC_LFP';

save_name{1} = '2007-5-24';
save_name{2} = '2007-5-31';
save_name{3} = '2007-6-4';
save_name{4} = '2007-7-3_A';
save_name{5} = '2007-7-3_B';
save_name{6} = '2007-7-4_B';
save_name{7} = '2007-7-4_D';
save_name{8} = '2007-7-5';





mp_min = -1.0;
mp_max = 0.1;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

Fs = 2e4;
% dsf = 20;
% Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[100/niqf 5e3/niqf]);
spkthresh = 0.01*2;

wait_samps = round(0.5*Fs);
ci_dur = Fs;

load old_1s_dep_data_minus_spike

for d = 1:8
    
    load(depol_dir{d}) 
    data = data_minus_spike{d};
    
    if length(data) ~= length(time)
        time(length(data)+1:end) = [];
    end
    
    send_times = find(time==1.9999);
    sbeg_times = find(time==0);

    
    
    num_sweeps(d) = length(send_times);
    stim_times = [];
    spkid = [];
    for i = 1:num_sweeps(d)
        cur_data = data(sbeg_times(i):send_times(i));
        diff_data = [0;diff(cur_data)];
        alen = length(cur_data);
        filt_cur_data = filtfilt(b,a,diff_data);
        cur_spkid = find((filt_cur_data(1:alen-1) < spkthresh) & (filt_cur_data(2:alen) >= spkthresh) );
        bad_spks = find(data(cur_spkid) > data(cur_spkid+30));
        cur_spkid(bad_spks) = [];
        spkid = [spkid;(sbeg_times(i)+cur_spkid)];
        stim_times = [stim_times (sbeg_times(i)+wait_samps):(sbeg_times(i)+wait_samps+ci_dur)];
    end
    stim_spikes = ismember(spkid,stim_times);
    num_stim_spikes = sum(stim_spikes);
    num_non_stim_spikes = length(spkid)-num_stim_spikes;
    stim_time = num_sweeps(d)*1;
    non_stim_time = length(data)/Fs - stim_time;
    stim_rate(d) = num_stim_spikes/stim_time;
    non_stim_rate(d) = num_non_stim_spikes/non_stim_time; 
    
d
end
% 
cd C:\WC_Germany\current_injection\persistent_ci
% save seprec_1s_ci_data tvec mean_* avg_amp* norm_mean* prior*


