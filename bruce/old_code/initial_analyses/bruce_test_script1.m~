clear all
close all

cd ~/Data/bruce/2_27_12

load Blocks
load lemM232A.51.lfp.mat

cur_block = 1;
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
block_trial_times = Blocks{cur_block}.blocktimes;

%%
cd saccades/
load lemM232.51.em.sac.mat

sac_times =  Expt.Trials.EyeMsac.sacT;
sac_endtimes = sac_times + Expt.Trials.EyeMsac.duration/1e3;

sac_times = sac_endtimes; %use end times
%% Find used saccades 
forward_lag = 0.75;
backward_lag = 0.25;

is_sac_used = zeros(size(sac_times));
for i = 1:length(sac_times)
   prev_lfp_start = find(lfp_trial_start < sac_times(i),1,'last');
   if sac_times(i) < lfp_trial_stop(prev_lfp_start) - forward_lag & sac_times(i) > lfp_trial_start(prev_lfp_start) + backward_lag
       is_sac_used(i) = 1;
   end
end
    
used_sac_times = sac_times(is_sac_used==1);
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
[b,a] = butter(2,[1 20]/niqf);
% lfp_samps = filtfilt(b,a,lfp_samps);

%%

lags = -backward_lag:.001:forward_lag;
sac_trig_avgs = zeros(24,length(lags));
for i = 1:length(used_sac_times)
   [~,match_samp] = min(abs(lfp_time-used_sac_times(i)));
   bs = match_samp - round(backward_lag*1000);
   es = match_samp + round(forward_lag*1000);
   sac_trig_avgs = sac_trig_avgs + lfp_samps(bs:es,:)'; 
   
%    plot(lags,lfp_samps(bs:es,[1 15]))
%    pause(2)
   
end
sac_trig_avgs = sac_trig_avgs/length(used_sac_times);