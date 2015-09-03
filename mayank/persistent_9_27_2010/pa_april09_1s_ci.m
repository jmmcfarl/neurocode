clear all

depol_array{1} = 'C:\wc_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_1s_POS_CI';
depol_array{2} = 'C:\wc_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_1s_POS';
hypol_array{1} = 'C:\wc_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_1s_NEG_CI';
hypol_array{2} = 'C:\wc_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_1s_NEG';

ci_time = 1;
wait_time = 2;

first_pulse = 1; %time of first CI

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
    
    %% get depolarizing CI
    load(depol_array{d})  
    time = (1:length(data))/Fs;
    diff_data = [0;diff(data)];
    alen = length(data);
    filt_data = zscore(filtfilt(b,a,data));
    
    spkid = find( (filt_data(1:alen-1) < spkthresh) & (filt_data(2:alen) >= spkthresh) );
    bad_spks = find(diff(spkid) < 0.001*Fs);
    spkid(bad_spks) = [];
    
    num_cis = floor((max(time)-first_pulse)/(ci_time+wait_time));
    stim_times = [];
    artifact_times = [];
    for i = 1:num_cis     
        begpt = (i-1)*(wait_time+ci_time)*Fs+1+round(first_pulse*Fs);
        endpt = begpt+ci_time*Fs;
        stim_times = [stim_times begpt:endpt];
        artifact_times = [artifact_times (begpt-25):begpt endpt:(endpt+25)]; %dont count stim artifacts as spikes
    end
    spkid(ismember(spkid,artifact_times)) = [];

    %interpolate to remove spikes
    spkRmForward = round(.003*Fs); %width of spike to remove
    spkRmBackward = round(0.0005*Fs);
    during_spk = zeros(size(data));
    for i = -spkRmBackward:spkRmForward
        during_spk(spkid+i) = 1;
    end
    data_minus_spike = interp1(find(~during_spk),data(~during_spk),1:length(data));
    
    %compute firing rate during stimulation and outside stimulation
    stim_spikes = ismember(spkid,stim_times);
    num_stim_spikes = sum(stim_spikes);
    num_non_stim_spikes = length(spkid)-num_stim_spikes;
    stim_time = num_cis*ci_time;
    non_stim_time = max(time) - stim_time;
    stim_rate(d) = num_stim_spikes/stim_time;
    non_stim_rate(d) = num_non_stim_spikes/non_stim_time;
   
    %%
    data = downsample(data,dsf);
    time = downsample(time,dsf);
    
    %compile stimulations into a matrix
    num_dep_cis(d) = floor((max(time)-first_pulse)/(ci_time+wait_time));
    sweep_time = (ci_time+wait_time)*Fsd;
    depol_mat = zeros(num_dep_cis(d),sweep_time);
    for i = 1:num_dep_cis(d)
        begpt = (i-1)*sweep_time+1;
        endpt = begpt+sweep_time;
        depol_mat(i,:) = data(begpt:endpt-1);
    end
    tvec = (1:sweep_time)/Fsd;
    
%     dep_end = find(tvec > 2.0,1,'first');
%     dep_start = find(tvec > 1.0,1,'first');
    mean_afterdep_traj(d,:) = mean(depol_mat);
    
    clear data time
    
    %% get hyperpolarizing CI
    load(hypol_array{d})  
    time = (1:length(data))/Fs;
    diff_data = [0;diff(data)];
    alen = length(data);
    filt_data = zscore(filtfilt(b,a,data));
    
    spkid = find( (filt_data(1:alen-1) < spkthresh) & (filt_data(2:alen) >= spkthresh) );
     bad_spks = find(diff(spkid) < 0.001*Fs);
    spkid(bad_spks) = [];
   
    num_cis = floor((max(time)-first_pulse)/(ci_time+wait_time));
    stim_times = [];
    artifact_times = [];
    for i = 1:num_cis     
        begpt = (i-1)*(ci_time+wait_time)*Fs+1+round(first_pulse*Fs);
        endpt = begpt+ci_time*Fs;
        stim_times = [stim_times begpt:endpt];
        artifact_times = [artifact_times (begpt-25):begpt endpt:(endpt+25)];
    end
    spkid(ismember(spkid,artifact_times)) = [];

    spkRmForward = round(.003*Fs); %width of spike to remove
    spkRmBackward = round(0.0005*Fs);
    during_spk = zeros(size(data));
    for i = -spkRmBackward:spkRmForward
        during_spk(spkid+i) = 1;
    end
    data_minus_spike = interp1(find(~during_spk),data(~during_spk),1:length(data));
    
    stim_spikes = ismember(spkid,stim_times);
    num_stim_spikes = sum(stim_spikes);
    num_non_stim_spikes = length(spkid)-num_stim_spikes;
    stim_time = num_cis*ci_time;
    non_stim_time = max(time) - stim_time;
    hstim_rate(d) = num_stim_spikes/stim_time;
    non_hstim_rate(d) = num_non_stim_spikes/non_stim_time;
   
    %%
    data = downsample(data,dsf);
    time = downsample(time,dsf);
    
    num_hcis(d) = floor((max(time)-first_pulse)/(ci_time+wait_time));
    sweep_time = (ci_time+wait_time)*Fsd;
    hypol_mat = zeros(num_hcis(d),sweep_time);
    for i = 1:num_hcis(d)
        begpt = (i-1)*sweep_time+1;
        endpt = begpt+sweep_time;
        hypol_mat(i,:) = data(begpt:endpt-1);
    end
    tvec = (1:sweep_time)/Fsd;
    dep_end = find(tvec > 2.0,1,'first');
    dep_start = find(tvec > 1.0,1,'first');
    mean_afterhyp_traj(d,:) = mean(hypol_mat);
    
    clear data time
end

cd C:\WC_Germany\persistent_9_27_2010\
save april09_1s_ci_data mean_* *rate num_cis num_hcis Fsd
