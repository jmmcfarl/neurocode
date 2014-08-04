  
clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')

Fs = 2016;

maxlag = 1*Fs;
lags = -maxlag:maxlag;
niqf = 2016/2;
lcf = 10/niqf;
[b,a] = butter(2,lcf,'high');
for d = 1:length(over_dir)
    d=48
    cd(over_dir{d})
    pwd

    load used_data wcv lf8 lf2 lf3 lf5
% wcv = filtfilt(b,a,wcv);
% lf8 = filtfilt(b,a,lf8);
% lf3 = filtfilt(b,a,lf3);
lf2 = filtfilt(b,a,lf2);
% lf5 = filtfilt(b,a,lf5);
% wcv = zscore(wcv);
% lf8 = zscore(lf8);
% lf3 = zscore(lf3);
lf2 = zscore(lf2);
% lf5 = zscore(lf5);
if exist('spike_time.mat')
        load spike_time   
    else 
        load spike_time_br
end

rec_dur = length(wcv)/Fs;
num_spikes(d) = length(spkid);
mean_rate(d) = num_spikes(d)/rec_dur;

%       spk_trig_wcv = zeros(length(spkid),length(lags));
%       spk_trig_lf8 = zeros(length(spkid),length(lags));
%       spk_trig_lf3 = zeros(length(spkid),length(lags));
      spk_trig_lf2 = zeros(length(spkid),length(lags));
%        spk_trig_lf5 = zeros(length(spkid),length(lags));
     
      for i = 1:length(spkid)
          
         if spkid(i) > maxlag && length(wcv)-spkid(i) > maxlag 
%              spk_trig_wcv(i,:) = wcv(spkid(i)-maxlag:spkid(i)+maxlag);
%              spk_trig_lf8(i,:) = lf8(spkid(i)-maxlag:spkid(i)+maxlag);
%              spk_trig_lf3(i,:) = lf3(spkid(i)-maxlag:spkid(i)+maxlag);
             spk_trig_lf2(i,:) = lf2(spkid(i)-maxlag:spkid(i)+maxlag);
%              spk_trig_lf5(i,:) = lf5(spkid(i)-maxlag:spkid(i)+maxlag);

         else
%              spk_trig_wcv(i,:) = nan;
%              spk_trig_lf8(i,:) = nan;
%              spk_trig_lf3(i,:) = nan;
             spk_trig_lf2(i,:) = nan;
%               spk_trig_lf5(i,:) = nan;
         
         end
          
      end

%       sta_wcv(d,:) = nanmean(spk_trig_wcv);
%       sta_lf8(d,:) = nanmean(spk_trig_lf8);
%       sta_lf3(d,:) = nanmean(spk_trig_lf3);
      sta_lf2(d,:) = nanmean(spk_trig_lf2);
%       sta_lf5(d,:) = nanmean(spk_trig_lf5);
      
%       plot(lags/Fs,sta_wcv(d,:))
%       xlim([-0.2 0.2]);grid
%       title(['total spikes: ' num2str(num_spikes(d)) '  mean rate: ' num2str(mean_rate(d))])
%       t_names = ['C:\WC_Germany\overall_calcs\sta\sta_mp_' num2str(cell_type(d)) '_' over_names{d}];
%         print('-dpng',t_names);
%         close
%         
%         plot(lags/Fs,sta_lf8(d,:),'r')
%         hold on
%         plot(lags/Fs,sta_lf3(d,:),'g')
%         plot(lags/Fs,sta_lf2(d,:),'k')
%         legend('LF8','LF3','LF2')
%               xlim([-0.2 0.2]);grid

        plot(lags/Fs,sta_lf2(d,:),'m')
              xlim([-0.2 0.2]);grid

%       title(['total spikes: ' num2str(num_spikes(d)) '  mean rate: ' num2str(mean_rate(d))])
%       t_names = ['C:\WC_Germany\overall_calcs\sta\sta_lf5_' num2str(cell_type(d)) '_' over_names{d}];
%         print('-dpng',t_names);
%         close      
      
end

% save C:\WC_Germany\overall_calcs\sta\sta_data5 lags Fs sta* mean_rate num_spikes