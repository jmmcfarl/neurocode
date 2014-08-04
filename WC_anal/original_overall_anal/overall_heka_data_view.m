load('C:\WC_Germany\JMM_analysis_ste\stellate_heka_dir.mat')
% load('C:\WC_Germany\JMM_analysis_pyr\pyr_heka_dir.mat')

d = 5;
load(f_loc{d})
dat_name = f_loc{d}(25:end);

eval(['data = ' dat_name '_MP;'])
eval(['time = ' dat_name '_sampletimes;'])
Fs = 5000;
niqf = Fs/2;
[b,a] = butter(2,2/niqf,'low');
[b2,a2] = butter(2,[7/niqf 20/niqf]);

sweep_t = find(diff(time) < 0);
low = min(data);
high = max(data)-0.3;
figure
dataseg = data(1:sweep_t(1));
filtseg = filtfilt(b,a,dataseg);
% plot(time(1:sweep_t(1)),data(1:sweep_t(1)))
% hold on
% plot(time(1:sweep_t(1)),filtseg,'r')
% ylim([low high])
% pause
% clf
for i = 1:length(sweep_t)-1
    dataseg = data(sweep_t(i)+1:sweep_t(i+1));
filtseg = filtfilt(b,a,dataseg);
% filtseg2 = filtfilt(b2,a2,dataseg);
dfilt = [0;diff(filtseg)];
dfilt = jmm_smooth_1d(dfilt,100);
subplot(2,1,1)
   plot(time(sweep_t(i)+1:sweep_t(i+1)),data(sweep_t(i)+1:sweep_t(i+1)))
   hold on
      plot(time(sweep_t(i)+1:sweep_t(i+1)),filtseg,'r','linewidth',2)
%       plot(time(sweep_t(i)+1:sweep_t(i+1)),filtseg2*3+mean(dataseg),'g','linewidth',2)
   ylim([low high])
      subplot(2,1,2)
plot(time(sweep_t(i)+1:sweep_t(i+1)),dfilt,'k')
line([time(sweep_t(i)+1) time(sweep_t(i+1))],[5e-5 5e-5],'Color','r')
line([time(sweep_t(i)+1) time(sweep_t(i+1))],[1e-5 1e-5],'Color','b')
line([time(sweep_t(i)+1) time(sweep_t(i+1))],[-1e-5 -1e-5],'Color','b')

   grid
   pause
   clf
   
end