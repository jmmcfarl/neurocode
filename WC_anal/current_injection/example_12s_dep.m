clear all
close all

cd C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19
dsf = 8;
Fs = 2016;
Fsd = Fs/dsf;

niqf = Fs/2;
load raw_data wcv_minus_spike lf8
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

lf8_f = filtfilt(b,a,lf8);
wcv_f = filtfilt(b,a,wcv_minus_spike);

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

t = (1:length(wcv_d))/2016*dsf;

offset = 1112.3 - 1.031-305.5750;

plot(t-offset,lf8_d,'r')
xlim([0 30])
figure
depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_12s_CI';
load(depol_array{1})
Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;
    data = downsample(data,dsf);
    time = downsample(time,dsf);
    load('C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\april_09_1s_dep_data')
d = 1
% data = (data - mean_prior(d))/std_prior(d);
plot(time,data)
% line([9 9],[-3 -2],'Color','k')
% line([21 21],[-3 -2],'Color','k')
% line([9 9],[-90 -60],'Color','k')
% line([21 21],[-90 -60],'Color','k')
xlim([0 30])
% ylim([-3 5.5])
% ylim([-90 20])

% down_mp = n_pr_down_mp(d);
% up_mp = n_pr_up_mp(d);
% down_mp = (down_mp - mean_prior(d))/std_prior(d);
% up_mp = (up_mp - mean_prior(d))/std_prior(d);
% line([0 30],[down_mp down_mp],'Color','k')
% line([0 30],[up_mp up_mp],'Color','k')