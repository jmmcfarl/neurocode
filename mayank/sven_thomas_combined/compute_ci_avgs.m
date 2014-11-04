
function [mean_sweep,sweep_mat,t_axis,stim_rate,non_stim_rate,ci_inds] = compute_ci_avgs(data,time)

Fs = 2e4;
dsf = 40;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[100/niqf 5e3/niqf]);
spkthresh = 3;

stim_dur = 12;
pause_dur = 18;
first_stim = 9;

diff_data = [0;diff(data)];
alen = length(data);
filt_data = zscore(filtfilt(b,a,data));

spkid = find( (filt_data(1:alen-1) < spkthresh) & (filt_data(2:alen) >= spkthresh) );

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
stim_rate = num_stim_spikes/stim_time;
non_stim_rate = num_non_stim_spikes/non_stim_time;

%     data = downsample(data,dsf);
data = downsample(data_minus_spike,dsf);
time = downsample(time,dsf);

num_cis = floor((max(time)-2*first_stim)/(stim_dur+pause_dur))+1;
fprintf('%d CIs detected\n',num_cis);

sweep_time = (stim_dur+pause_dur)*Fsd;
sweep_mat = zeros(num_cis,sweep_time);
for i = 1:num_cis
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

mean_sweep = nanmean(sweep_mat);
t_axis = (1:sweep_time)/Fsd;
