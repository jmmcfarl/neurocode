clear all
close all

cd G:\WC_Germany\parietal_cortical_2010
load parietal_cortical_2010

raw_Fs = 2016;
dsf = 8;
niqf = raw_Fs/2;
Fsd = raw_Fs/dsf;

lf_lcf = 0.05;
lf_hcf = 2;
lcf = 0.05;
hcf = 20;
[b,a] = butter(2,[lcf/niqf hcf/niqf]);
[b_lf,a_lf] = butter(2,[lf_lcf/niqf lf_hcf/niqf]);

d = 16;

cd(sess_data(d).directory)
pwd

load used_data lf8 wcv_minus_spike
obs_lf = filtfilt(b_lf,a_lf,lf8);
obs_lf = zscore(downsample(obs_lf,dsf));
obs = filtfilt(b,a,lf8);
obs = zscore(downsample(obs,dsf));
obsw = filtfilt(b,a,wcv_minus_spike);
obsw = zscore(downsample(obsw,dsf));
time = (1:length(obs))/Fsd;

load hsmm_state_seq_lf
load hsmm_state_seq8_lf
timed = (1:length(hsmm_state_seq))/Fs;
timed(1) = time(1);
timed(end) = time(end);
hsmm_state_seq = round(interp1(timed,hsmm_state_seq,time));

seg = find(time >= 100 & time <= 300);

[hf_features,t_axis,Fs] = get_hf_features(lf8,raw_Fs,Fsd);
figure
plot(time(seg),obs(seg)), hold on
plot(time(seg),obs_lf(seg),'r')
plot(time(seg),hf_features(seg),'k')
plot(time(seg),hsmm_bbstate_seq(seg),'c','linewidth',2)

% figure
% plot(time,obsw), hold on
% plot(time,hsmm_bbstate_seq,'k','linewidth',2)
% plot(time,hsmm_state_seq-2,'c','linewidth',2)
% xlim([513,522])