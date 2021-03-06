clear all
close all

cd ~/Data/bruce/2_27_12
    load Blocks
all_S = [];
for blockid = 1:5;
    fprintf('Block %d of %d\n',blockid,5);
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
    [b,a] = butter(2,[1 200]/niqf);
    lfp_samps = filtfilt(b,a,lfp_samps);
    
    dsf = 2;
    lfp_sampsd = downsample(lfp_samps,dsf);
    Fsd = Fs/dsf;
    
    %%
%     win = 20;
    params.Fs = Fsd;
    params.tapers = [3 5];
    movingwin = [20 10];
    n_lfps = size(lfp_samps,2);
%     for i = 1:n_lfps
        [S{blockid},t,f]=mtspecgramc(lfp_sampsd,movingwin,params);
        n_pt(blockid) = length(t);
%     end
    all_S = cat(1,all_S, S{blockid});
end
%%
S_n = log(all_S);
% S_n = all_S;
S_n = bsxfun(@minus,S_n,nanmean(S_n));
S_n = bsxfun(@rdivide,S_n,nanstd(S_n));