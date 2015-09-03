clear all
close all


drive_letter = 'F';
cd F:\WC_Germany\overall_EC
load overall_EC_dir
%%
dsf = 16;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.05/niqf 4/niqf]);

d = find_struct_field_vals(sess_data,'name','2006-04-07_C')
% d = find_struct_field_vals(sess_data,'name','2009_05-10_A')
% d = find_struct_field_vals(sess_data,'name','2005_12-09_A')

cdir = sess_data(d).directory;
cdir(1) = 'F';
disp(sprintf('session %d',d))
cd(cdir);
s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
load used_data lf2 lf5 lf8 wcv_minus_spike
load spike_time_jmm
spk_id = round(spkid/dsf);

wcv = filtfilt(b,a,wcv_minus_spike);
wcv = zscore(downsample(wcv,dsf));
lf8 = filtfilt(b,a,lf8);
lf8 = zscore(downsample(lf8,dsf));
lf2 = lf2/sess_data(d).gains(2);
lf5 = lf5/sess_data(d).gains(5);
lf2_r = lf2-lf5;
lf2_r = filtfilt(b,a,lf2_r);
lf2_r = zscore(downsample(lf2_r,dsf));

t_axis = (1:length(lf8))/Fsd;

plot(t_axis,wcv+2), hold on
plot(t_axis,lf8,'r')
plot(t_axis,lf2_r*1.2-2,'k')
% xlim([560 576]) %for 2007-06-04_B
% xlim([617 630]) %for 2007-06-04_B
% xlim([888 912])
xlabel('Time (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)
legend('MP','Cortical LFP','HPC LFP')

% plot(t_axis,wcv,'g'), hold on
% plot(t_axis,lf8,'r')
% plot(t_axis,lf2_r,'k')
% % xlim([198 220]) % for 2009-05-10_B
% xlim([996 1016]) %for 2009-05-10_A
% xlabel('Time (s)','fontsize',14)
% ylabel('Amplitude (z)','fontsize',14)
% legend('MP','Cortical LFP','HPC LFP')
