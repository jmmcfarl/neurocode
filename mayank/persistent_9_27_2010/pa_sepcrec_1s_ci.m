clear all

depol_dir{1} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-05-24_CWC_LFP_B';
depol_dir{2} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-05-31_CWC_LFP';
depol_dir{3} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-06-04_CWC_LFP_B';
depol_dir{4} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-03_CWC_LFP_A';
depol_dir{5} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-03_CWC_LFP_B';
depol_dir{6} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-04_CWC_LFP_B';
depol_dir{7} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-04_CWC_LFP_D';
depol_dir{8} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-05_CWC_LFP';
depol_dir{9} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-05-25_CWC_LFP_B';

% prior_dir{1} = 'C:\wc_data\MPasci\A2007_05_24_CWC_LFP_B';
% prior_dir{2} = 'C:\wc_data\MPasci\A2007_05_31_CWC_LFP';
% prior_dir{3} = 'C:\wc_data\MPasci\A2007_06_04_CWC_LFP_B';
% prior_dir{4} = 'C:\wc_data\MPasci\A2007_07_03_CWC_LFP_A';
% prior_dir{5} = 'C:\wc_data\MPasci\A2007_07_03_CWC_LFP_B';
% prior_dir{6} = 'C:\wc_data\MPasci\A2007_07_04_CWC_LFP_B';
% prior_dir{7} = 'C:\wc_data\MPasci\A2007_07_04_CWC_LFP_D';
% prior_dir{8} = 'C:\wc_data\MPasci\A2007_07_05_CWC_LFP';
% prior_dir{9} = 'C:\wc_data\MPasci\A2007_05_25_CWC_LFP_B';

save_name{1} = '2007-5-24'; %4
save_name{2} = '2007-5-31'; %6
save_name{3} = '2007-6-4'; %8
save_name{4} = '2007-7-3_A'; %11
save_name{5} = '2007-7-3_B'; %12
save_name{6} = '2007-7-4_B'; %13
save_name{7} = '2007-7-4_D'; %15
save_name{8} = '2007-7-5'; %16
save_name{9} = '2007-5-25_B'; %22

mp_min = -1.0;
mp_max = 0.1;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[100/niqf 5e3/niqf]);
spkthresh = 5;

wait_samps = round(0.5*Fs);
ci_dur = Fs;

tvec = 0:1/Fsd:(2-1/Fsd);
dep_end = find(tvec > 1.6,1,'first');
dep_start = find(tvec > 0.5,1,'first');
stim_time = 60;
non_stim_time = 60;

for d = 1:9
    
    load(depol_dir{d})
    fdata = filtfilt(b,a,[0;diff(data(:))]);
    f_std = std(fdata);

    %if there was no time variable in the data set, create a uniform time
    %grid
    if ~exist('time','var')
        time = 0:(length(data));
        time = mod(time/Fs,2-1e-13);
        time = time(:);
        if time(end) < 2-1/Fs
            time = [time; (time(end)+1/Fs:1/Fs:(2-1/Fs))'];
        end
    end
    time = time(:);
    data = data(:);

    %locate the times corresponding to start and stop of sweeps
    dtimes = [Inf; diff(time)];
    sbeg_times = find(dtimes<0);
    send_times = sbeg_times - 1;
    sbeg_times = [1; sbeg_times];
    send_times = [send_times; length(time)];
    
    %in the data and time vectors are mismatched, adjust accordingly
    if length(data) < length(time)
        data = [data; nan(length(time)-length(data),1)];
    elseif length(data) > length(time)
        time(length(data)+1:end) = [];
    end

    num_sweeps(d) = length(send_times);
    if d < 9
        dep_sweeps = 61:120;
        hyp_sweeps = 1:60;
    else
        dep_sweeps = 1:60;
        hyp_sweeps = [];
    end
    
    stim_times = [];
    artifact_times = [];
    spkid = [];
    for i = dep_sweeps
        cur_data = data(sbeg_times(i):send_times(i));
        diff_data = [0;diff(cur_data)];
        alen = length(cur_data);
        filt_cur_data = filtfilt(b,a,diff_data)/f_std;
        cur_spkid = find((filt_cur_data(1:alen-1) < spkthresh) & (filt_cur_data(2:alen) >= spkthresh) );
%         bad_spks = find(cur_data(cur_spkid) > cur_data(cur_spkid+30)); %the stim artifacts have a big backswing, very diff from real spikes
%         cur_spkid(bad_spks) = [];
        spkid = [spkid;(sbeg_times(i)+cur_spkid)];
        artifact_times = [artifact_times (sbeg_times(i)+wait_samps-25):(sbeg_times(i)+wait_samps+25) ...
            (sbeg_times(i)+wait_samps+ci_dur):(sbeg_times(i)+wait_samps+ci_dur+25)];
        stim_times = [stim_times (sbeg_times(i)+wait_samps):(sbeg_times(i)+wait_samps+ci_dur)];
    end
    spkid(ismember(spkid,artifact_times)) = [];
    bad_spks = find(diff(spkid) < 0.002*Fs);
    spkid(bad_spks) = [];

    stim_spikes = ismember(spkid,stim_times); %spikes occuring during depolarization
    num_stim_spikes = sum(stim_spikes);
    num_non_stim_spikes = length(spkid)-num_stim_spikes;
    stim_time = length(dep_sweeps)*1; %total stimulation duration
%     non_stim_time = length(data)/Fs - stim_time;
    stim_rate(d) = num_stim_spikes/stim_time; %firing rate during dep.
%     non_stim_rate(d) = num_non_stim_spikes/non_stim_time; %firing rate outside dep.
    
    %interpolate data to remove spikes
    spkRmForward = round(.003*Fs); %width of spike to remove
    spkRmBackward = round(0.0005*Fs); %remove this much preceding spike triggering
    data_minus_spike = data;
    during_spk = zeros(size(data));
    for i = -spkRmBackward:spkRmForward
        during_spk(spkid+i) = 1;
    end
    data_minus_spike = interp1(find(~during_spk),data(~during_spk),1:length(data));
 
%%
    %store the raw data across sweeps in a matrix
    sweep_mat = zeros(num_sweeps(d),length(tvec));
    for i = 1:num_sweeps(d)
        curdata = data_minus_spike(sbeg_times(i):send_times(i));
        sweep_mat(i,:) = downsample(curdata,dsf);
    end
    mean_afterdep_traj(d,:) = nanmean(sweep_mat(dep_sweeps,:));
    mean_afterhyp_traj(d,:) = nanmean(sweep_mat(hyp_sweeps,:));
    
    avg_amp_before_dep(d) = nanmean(nanmean(sweep_mat(dep_sweeps,1:dep_start-100)));
    avg_amp_during_dep(d) = nanmean(nanmean(sweep_mat(dep_sweeps,dep_start:dep_end)));
    avg_amp_before_hyp(d) = nanmean(nanmean(sweep_mat(hyp_sweeps,1:dep_start-100)));
    avg_amp_during_hyp(d) = nanmean(nanmean(sweep_mat(hyp_sweeps,dep_start:dep_end)));
        
%     imagesc(sweep_mat)
%     pause
%     clf
%     d
    
    clear data time
end
%
cd F:\WC_Germany\persistent_9_27_2010\
save seprec_1s_ci_data mean_* avg_amp* *rate


