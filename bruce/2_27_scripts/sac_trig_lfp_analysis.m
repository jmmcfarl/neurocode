clear all
close all

cd ~/Data/bruce/2_27_12
load Blocks
cd /Users/James/Data/bruce/2_27_12/stimrecon
load saccade_data

%%
cd ~/Data/bruce/2_27_12
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    
    %%
    %get start times of each LFP trial
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    % lfp_trial_stop = lfp_trial_start+lfp_dur';
    
    %get block start times
    block_trial_times = Blocks{blockid}.blocktimes;
    
    %%
    ov_t = block_trial_times(1,1):.001:block_trial_times(2,end);
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
    end
    % lfp_samps = zscore(lfp_samps);
    
    Fs = 1000;
    niqf = Fs/2;
    [b,a] = butter(2,[1 50]/niqf);
    lfp_samps = filtfilt(b,a,lfp_samps);
    
    dsf = 6;
    lfp_sampsd = downsample(lfp_samps,dsf);
    lfp_timed = downsample(lfp_time,dsf);
    lfp_sampsd = zscore(lfp_sampsd);
    Fsd = Fs/dsf;
    
    %%
    n_channels = size(lfp_samps,2);
    cur_nsacs = length(sac_stats(blockid).sac_peakvel);
    maxlag = round(Fsd*0.3);
    sac_trig_avgs = zeros(n_channels,2*maxlag+1);
    sac_trig_peaks = zeros(length(cur_nsacs),n_channels);
    
    cur_sac_start_times = sac_stats(blockid).sac_start_time;
    sac_start_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_sac_start_times));
    used_sacs = find(sac_start_inds > maxlag & sac_start_inds < (length(lfp_timed) - maxlag));
    for i = 1:length(used_sacs)
        cur_set = (sac_start_inds(used_sacs(i))-maxlag):(sac_start_inds(used_sacs(i))+maxlag);
        sac_trig_avgs = sac_trig_avgs + lfp_sampsd(cur_set,:)';
        sac_trig_peaks(used_sacs(i),:) = max(lfp_sampsd(sac_start_inds(used_sacs(i)):sac_start_inds(used_sacs(i))+maxlag,:));
    end
    sac_trig_avgs = sac_trig_avgs/length(used_sacs);
    lags = (-maxlag:maxlag)/Fsd;
    
end
