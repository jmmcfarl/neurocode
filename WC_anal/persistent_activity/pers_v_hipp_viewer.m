clear all
close all
load C:\WC_Germany\Persistent_activity\dir_tree_update

Fs = 2016;
niqf = Fs/2;
dsf = 4;
Fsd = Fs/dsf;

lcf = 0.05/niqf;
hcf = 8/niqf;
[blow,alow] = butter(2,[lcf hcf]);

lcf = 2/niqf;
hcf = 10/niqf;
[bmid,amid] = butter(2,[lcf hcf]);

lcf = 100/niqf;
hcf = 250/niqf;
[bhigh,ahigh] = butter(2,[lcf hcf]);


win = 20;

for d = 1:length(dir_array)
    
   cd(dir_array{d})
   pwd
   
   load used_data lf2 lf3 lf8 wcv_minus_spike
   if exist('spike_time.mat') > 0
       load spike_time
   else
    load spike_time_br
   end
   lf3_theta = filtfilt(bmid,amid,lf3);
   lf3_ripple = filtfilt(bhigh,ahigh,lf3);
   lf2_uds = filtfilt(blow,alow,lf2);
   
   lf3_uds = filtfilt(blow,alow,lf3);

   lf8_uds = filtfilt(blow,alow,lf8);
   wcv = filtfilt(blow,alow,wcv_minus_spike);
   
   lf3_theta = downsample(lf3_theta,dsf);
   lf3_ripple = downsample(lf3_ripple,dsf);
   lf2_uds = downsample(lf2_uds,dsf);
   lf3_uds = downsample(lf3_uds,dsf);
   lf8_uds = downsample(lf8_uds,dsf);
   wcv = downsample(wcv,dsf);
   
   spike_id = round(spkid/dsf);
   
   lf3_theta = zscore(lf3_theta);
   lf3_ripple = zscore(lf3_ripple);
   lf2_uds = zscore(lf2_uds);
   lf8_uds = zscore(lf8_uds);
   lf3_uds = zscore(lf3_uds);
   wcv = zscore(wcv);
   
   t_axis = (1:length(wcv))/Fsd;
   
%     plot(t_axis,lf2_theta,'k')
      hold on
    plot(t_axis,wcv,'linewidth',4)
   plot(t_axis,lf8_uds-1,'r','linewidth',4)
%    plot(t_axis,lf2_uds-1,'c','linewidth',2)
%    plot(t_axis,lf3_uds-1,'k','linewidth',2)
   plot(t_axis,lf3_ripple,'k','linewidth',2)
    num_wins = ceil(max(t_axis)/win);
    plot(t_axis(spike_id),ones(size(spike_id))*1,'g*')
    for n = 1:num_wins
       xlim([1+(n-1)*win n*win])
       ylim([-5 5])
        pause
    end
    
end