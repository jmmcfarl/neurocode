clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update

niqf = 2016/2;
% lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,hcf,'low');
dsf = 8;

for d = 1:17
    
   cd(dir_array{d})
   pwd
   
   load used_data lf8 wcv_minus_spike
   
   lf8_f = filtfilt(b,a,lf8);
   wcv_f = filtfilt(b,a,wcv_minus_spike);
   
   lf8_d = downsample(lf8_f,dsf);
   wcv_d = downsample(wcv_f,dsf);
   
   lf8_d = zscore(lf8_d);
   wcv_d = zscore(wcv_d);
   
    [lf8_dist(d,:),gridv] = gpkde(lf8_d,-3,[-4;4;400]);
    [wcv_dist(d,:),gridv] = gpkde(wcv_d,-3,[-4;4;400]);
    
 plot(gridv,wcv_dist(d,:),'linewidth',2)
 hold on
 plot(gridv,lf8_dist(d,:),'r','linewidth',2)
 legend('MP','LFP')
    tname = ['C:\WC_Germany\Persistent_activity\neural_amp\' f_names{d}];
    print('-dpng',tname);
    close
   
end

save C:\WC_Germany\Persistent_activity\neural_amp\low_pass_data_3 lf8_dist gridv wcv_dist