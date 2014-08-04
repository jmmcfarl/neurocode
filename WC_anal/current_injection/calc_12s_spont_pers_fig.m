
clear all

cd C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19
dsf = 8;
Fs = 2016;
Fsd = Fs/dsf;

niqf = Fs/2;
load used_data wcv_minus_spike lf8
load spike_time_jmm
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

spike_ids = round(spkid/dsf);

lf8_f = filtfilt(b,a,lf8);
wcv_f = filtfilt(b,a,wcv_minus_spike);

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

t = (1:length(wcv_d))/2016*dsf;
loc = 690;
% loc = 302;


dur = 12;

load 2009-04-07_CWC_LFP_spontaneous
Fs = 2e4;
dsf = 5;
Fsd = Fs/dsf;
data = downsample(data,dsf);
time = downsample(time,dsf);
offset = 692.6-683.1
plot(time+offset-loc,data,'k')
xlim([0 dur])
ylim([-90 20])
load('C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\april_09_1s_dep_data')
d = 1
down_mp = n_pr_down_mp(d);
up_mp = n_pr_up_mp(d);
line([0 dur],[down_mp down_mp],'Color','k')
line([0 dur],[up_mp up_mp],'Color','k')


figure
plot(t-loc,lf8_d,'r','linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
xlim([0 dur])
ylim([-2 2.5])


