clear all
cd C:\WC_Germany\persistent_revised
load pers_revised_dir

d = 14;

cd(dir_array{d})
pwd

Fs = 2016;
niqf = Fs/2;
lcf = 0.1/niqf;
hcf = 4/niqf;
dsf = 8;
Fsd = Fs/dsf;

[b,a] = butter(2,[lcf hcf]);

load used_data wcv_minus_spike lf8 

wcv_f = filtfilt(b,a,wcv_minus_spike);
lf8_f = filtfilt(b,a,lf8);
% lf3_f = filtfilt(b,a,lf3);
% lf2_f = filtfilt(b,a,lf2);

wcv_f = zscore(downsample(wcv_f,dsf));
lf8_f = zscore(downsample(lf8_f,dsf));
% lf3_f = zscore(downsample(lf3_f,dsf));
% lf2_f = zscore(downsample(lf2_f,dsf));

t = (1:length(wcv_f))/Fsd;

figure
plot(t,wcv_f), hold on
plot(t,lf8_f,'r')
% plot(t,lf3_f,'k')
% xlim([710 730])