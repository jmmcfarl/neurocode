%% Cell 1
clear all, close all
cd C:\WC_Germany\persistent_revised
load pers_revised_dir

%PA
% sp = 713;
% ep = 733;
% line_bottom = 1.5;
% line_top = 2.0;


%NPA
sp = 35;
ep = 47;
line_bottom = 2.3;
line_top = 2.8;

d = 19;

cd(dir_array{d})
pwd

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
dsf = 8;
Fsd = Fs/dsf;

[b,a] = butter(2,[lcf hcf]);

load used_data wcv_minus_spike lf8 wcv
load spike_time_jmm
spike_ids = round(spkid/dsf);

wcv_f = filtfilt(b,a,wcv_minus_spike);
% wcv_f = wcv;
lf8_f = filtfilt(b,a,lf8);
% lf3_f = filtfilt(b,a,lf3);
% lf2_f = filtfilt(b,a,lf2);

wcv_f = zscore(downsample(wcv_f,dsf));
lf8_f = zscore(downsample(lf8_f,dsf));
% lf3_f = zscore(downsample(lf3_f,dsf));
% lf2_f = zscore(downsample(lf2_f,dsf));
t = (1:length(wcv_f))/Fsd;

used_spikes = find(spike_ids/Fsd > sp & spike_ids/Fsd < ep);


figure
plot(t-sp,wcv_f)
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i)))-sp t(spike_ids(used_spikes(i)))-sp],...
       [line_bottom line_top],'Color','k','linewidth',1) 
end
figure(2)
plot(t-sp,lf8_f,'r')
% plot(t,lf3_f,'k')
xlim([0 ep-sp]), figure(1), xlim([0 ep-sp])

cd C:\WC_Germany\persistent_revised\figs\
%% Cell 2
clear all, close all
cd C:\WC_Germany\persistent_revised
load pers_revised_dir

%PA
% sp = 424.5;
% ep = 449.8;
% line_bottom = 1.5;
% line_top = 2.0;


%NPA
sp = 854.5;
ep = sp+8.5;
line_bottom = 1.3;
line_top = 1.8;
% 
d = 14;

cd(dir_array{d})
pwd

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
dsf = 8;
Fsd = Fs/dsf;

[b,a] = butter(2,[lcf hcf]);

load used_data wcv_minus_spike lf8 wcv
load spike_time_jmm
spike_ids = round(spkid/dsf);

wcv_f = filtfilt(b,a,wcv_minus_spike);
% wcv_f = wcv;
lf8_f = filtfilt(b,a,lf8);
% lf3_f = filtfilt(b,a,lf3);
% lf2_f = filtfilt(b,a,lf2);

wcv_f = zscore(downsample(wcv_f,dsf));
% wcv_f(wcv_f >1.7) = 1.7;

lf8_f = zscore(downsample(lf8_f,dsf));
% lf3_f = zscore(downsample(lf3_f,dsf));
% lf2_f = zscore(downsample(lf2_f,dsf));
t = (1:length(wcv_f))/Fsd;

used_spikes = find(spike_ids/Fsd > sp & spike_ids/Fsd < ep);


figure
plot(t,wcv_f)
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
       [line_bottom line_top],'Color','k','linewidth',1) 
end
figure(2)
plot(t,lf8_f,'r')
% plot(t,lf3_f,'k')
xlim([sp ep]), figure(1), xlim([sp ep])

cd C:\WC_Germany\persistent_revised\figs\
%% Cell 3
clear all, close all
cd C:\WC_Germany\persistent_revised
load pers_revised_dir

%PA
% sp = 213;
% ep = 280;
% line_bottom = 1.7;
% line_top = 2.2;


%NPA
sp = 448;
ep = 465;
line_bottom = 2.5;
line_top = 3.0;

d = 5;

cd(dir_array{d})
pwd

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
dsf = 8;
Fsd = Fs/dsf;

[b,a] = butter(2,[lcf hcf]);

load used_data wcv_minus_spike lf8 wcv
load spike_time_jmm
spike_ids = round(spkid/dsf);

wcv_f = filtfilt(b,a,wcv_minus_spike);
% wcv_f = wcv;
lf8_f = filtfilt(b,a,lf8);
% lf3_f = filtfilt(b,a,lf3);
% lf2_f = filtfilt(b,a,lf2);

wcv_f = zscore(downsample(wcv_f,dsf));
% wcv_f(wcv_f > 1.) = 1.;
lf8_f = zscore(downsample(lf8_f,dsf));
% lf3_f = zscore(downsample(lf3_f,dsf));
% lf2_f = zscore(downsample(lf2_f,dsf));
t = (1:length(wcv_f))/Fsd;

used_spikes = find(spike_ids/Fsd > sp & spike_ids/Fsd < ep);


figure
plot(t-sp,wcv_f)
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i)))-sp t(spike_ids(used_spikes(i)))-sp],...
       [line_bottom line_top],'Color','k','linewidth',1) 
end
figure(2)
plot(t-sp,lf8_f,'r')
% plot(t,lf3_f,'k')
xlim([0 ep-sp]), figure(1), xlim([0 ep-sp])

cd C:\WC_Germany\persistent_revised\figs\