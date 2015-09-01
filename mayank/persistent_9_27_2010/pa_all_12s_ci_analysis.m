clear all

mp_min = -1.0;
mp_max = 0.1;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

Fs = 500;
niqf = Fs/2;
[b,a] = butter(2,100/niqf,'high');
spkthresh = 1.5;


%% FOR OLD DATA
old_depol_array{1} = 'G:\WC_Germany\longdepol\2007-05-31_CWC_LFP\2007-05-31_CWC_LFP_Hekadata\sweep_data'; %6
old_depol_array{2} = 'G:\WC_Germany\longdepol\2007-06-04_CWC_LFP_B\2007-06-04_CWC_LFP_B_Hekadata\sweep_data'; %8
stim_beg = round(1*Fs); %index of stimulation start within each sweep
stim_end = round(13*Fs); %index of stimulation stop

for d = 1:2
    spkid = [];
    stim_times = [];
    artifact_times = [];
    load(old_depol_array{d})
    temp = sweep_data';
    data = temp(:);
    fdata = filtfilt(b,a,[0;diff(data)]);
    f_std = std(fdata);

    [num_sweeps,sweep_dur] = size(sweep_data);
    for s = 1:num_sweeps
        cur_data = sweep_data(s,:);
        diff_data = [0 diff(cur_data)];
        alen = length(cur_data);
        filt_cur_data = zscore(filtfilt(b,a,diff_data));
        cur_spkid = find((filt_cur_data(1:alen-1) < spkthresh) & (filt_cur_data(2:alen) >= spkthresh) );
        spkid = [spkid ((s-1)*size(sweep_data,2)+cur_spkid)];
        stim_times = [stim_times ((s-1)*sweep_dur+stim_beg):((s-1)*sweep_dur+stim_end)];
        artifact_times = [artifact_times ((s-1)*sweep_dur+stim_beg-4):((s-1)*sweep_dur+stim_beg+4) ...
            ((s-1)*size(sweep_data,2)+stim_end-4):((s-1)*size(sweep_data,2)+stim_end+4)]; %identify then remove stimulus artifact times so as not to count as spikes
    end
    spkid(ismember(spkid,artifact_times)) = [];
    bad_spks = find(diff(spkid) < 0.002*Fs);
    spkid(bad_spks) = [];
    
    data = reshape(sweep_data',numel(sweep_data),1);
    
    spkRmForward = round(.003*Fs); %width of spike to remove
    spkRmBackward = round(0.0005*Fs);
    during_spk = zeros(size(data));
    for i = -spkRmBackward:spkRmForward
        during_spk(spkid+i) = 1;
    end
    data_minus_spike = interp1(find(~during_spk),data(~during_spk),1:length(data));
    data_minus_spike = reshape(data_minus_spike,size(sweep_data'));
    data_minus_spike = data_minus_spike';
    
    stim_spikes = ismember(spkid,stim_times);
    num_stim_spikes = sum(stim_spikes);
    num_non_stim_spikes = length(spkid)-num_stim_spikes;
    stim_time = size(sweep_data,1)*12;
    non_stim_time = size(sweep_data,1)*4;
    stim_rate(d) = num_stim_spikes/stim_time;
    non_stim_rate(d) = num_non_stim_spikes/non_stim_time;
    
%     mean_sweep(d,:) = mean(sweep_data);
    mean_sweep(d,:) = mean(data_minus_spike);
    num_cis(d) = size(sweep_data,1);
    
end


%% FOR NEW DATA
Fs = 2e4;
dsf = 40;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[100/niqf 5e3/niqf]);
spkthresh = 3;

stim_dur = 12;
pause_dur = 18;
first_stim = 9;

depol_array{1} = 'G:\wc_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_12s_CI'; %18
depol_array{2} = 'G:\wc_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_12s_CI'; %19
depol_array{3} = 'G:\wc_data\2009-04-13_B\2009-04-13_CWC_LFP_B_12s_CI_smaller'; %20

for d = 1:3
    
    load(depol_array{d})
    time = (1:length(data))/Fs;
    diff_data = [0;diff(data)];
    alen = length(data);
    filt_data = zscore(filtfilt(b,a,data));
    
    spkid = find( (filt_data(1:alen-1) < spkthresh) & (filt_data(2:alen) >= spkthresh) );
    % bad_spks = find(data(spkid) > data(spkid+30));
    % spkid(bad_spks) = [];
    
    %some broad spike bursts can trigger multiple spikes, eliminate these
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
    bad_spks = find(diff(spkid) < 0.001*Fs);
    spkid(bad_spks) = [];

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
    stim_time = num_cis*12;
    non_stim_time = num_cis*18;
    stim_rate(d+2) = num_stim_spikes/stim_time;
    non_stim_rate(d+2) = num_non_stim_spikes/non_stim_time;
    
%     data = downsample(data,dsf);
    data = downsample(data_minus_spike,dsf);
    time = downsample(time,dsf);
    
    num_cis(d+2) = floor((max(time)-2*first_stim)/(stim_dur+pause_dur))+1;
    
    sweep_time = (stim_dur+pause_dur)*Fsd;
    sweep_mat = zeros(num_cis(d+2),sweep_time);
    
    for i = 1:num_cis(d+2)
        begpt = (i-1)*30*Fsd+1;
        endpt = begpt+30*Fsd;
        if endpt > length(data)
            endpt = length(data);
            cur_length = endpt-begpt+1;
            sweep_mat(i,1:cur_length) = data(begpt:endpt);
            sweep_mat(i,cur_length+1:end) = nan;
        else
            sweep_mat(i,:) = data(begpt:endpt-1);
        end
    end
    
    new_mean_sweep(d,:) = nanmean(sweep_mat);
    t_axis = (1:sweep_time)/Fsd;
end

cd G:\WC_Germany\persistent_9_27_2010\
% save all_12s_CI_data_minus_spk *mean_sweep num_cis *rate t_axis Fsd
