clear all
close all

%% load heka data
cd C:\wc_data\2009-04-07\2009-4-7-19
dsf = 8;
Fs = 2016;
Fsd = Fs/dsf;
niqf = Fs/2;
load ./raw_data wcv lf8
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

lf8_f = filtfilt(b,a,lf8);
lf8_d = downsample(lf8_f,dsf);
lf8_d = zscore(lf8_d);

nlx_t_d = (1:length(lf8_d))/Fsd;
nlx_t = (1:length(wcv))/Fs;

base_nlxoffset = 800;

%% load Heka data
junction_pot = 7;
depol_array{1} = 'C:\wc_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_12s_CI';
load(depol_array{1})
Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;
data = downsample(data,dsf);
time = downsample(time,dsf);

%% align data
approx_nlx_times = find(nlx_t >= base_nlxoffset & nlx_t < base_nlxoffset + max(time));
[dc_offset,dc_maxcorr] = align_dc_ac_sigs_initial_v2(data,time,wcv(approx_nlx_times));
nlx_d = wcv(approx_nlx_times);
nlx_t = nlx_t(approx_nlx_times)- base_nlxoffset - dc_offset;

nlx_t_d = nlx_t_d - base_nlxoffset - dc_offset;

%% plot data
figure
plot(nlx_t_d,lf8_d,'r')
xlim([0 30])

figure
plot(time,data-junction_pot)
xlim([0 30])
x = 0:.001:30;
s = zeros(size(x));
s(x >= 9 & x <= 21) = 1;
hold on
plot(x,s*20-60,'k')