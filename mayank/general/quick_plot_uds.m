function [] = quick_plot_uds(mp,lfp)

Fs = 2016;
dsf = 8;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

mp_f = filtfilt(b,a,mp);
lfp_f = filtfilt(b,a,lfp);

mp_d = downsample(mp_f,dsf);
lfp_d = downsample(lfp_f,dsf);

Fsd = Fs/dsf;
t = (1:length(mp_d))/Fsd;

plot(t,mp_d)
hold on
plot(t,lfp_d,'r')