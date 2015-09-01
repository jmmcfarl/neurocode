%% Initialize variables
clear all
close all

load used_data
load spike_time_jmm
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data

cur_up_trans = up_trans{1}(synch_ups{1});
cur_down_trans = down_trans{1}(synch_downs{1});

Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;

niqf = Fs/2;
lcf = 0.1/niqf;

[b,a] = butter(2,lcf,'high');

lf3_f = filtfilt(b,a,lf3);
lf3_f = downsample(lf3_f,dsf);
lf3_f = zscore(lf3_f);

lf5_f = filtfilt(b,a,lf5);
lf5_f = downsample(lf5_f,dsf);
lf5_f = zscore(lf5_f);

lf7_f = filtfilt(b,a,lf7);
lf7_f = downsample(lf7_f,dsf);
lf7_f = zscore(lf7_f);

lf8_f = filtfilt(b,a,lf8);
lf8_f = downsample(lf8_f,dsf);
lf8_f = zscore(lf8_f);

lfp_mat = [lf3_f lf5_f lf7_f lf8_f];

spkid = round(spkid/dsf);

[sta,U,sing_vals,V,lags] = get_event_trig_stats_v2(spkid,lfp_mat,1*Fsd,0,Fsd,5);
