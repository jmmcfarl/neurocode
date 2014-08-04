clear all
close all

cd C:\WC_Germany\Cortical_analysis
load cortical_dir
load C:\WC_Germany\Cortical_analysis\sigmoid_fit\sig_fit_alldata

Fs = 2016;
nbins = 100;
lag_range = linspace(-0.5,0.5,nbins);
slope_range = linspace(0,100,nbins);
sig_width_range = linspace(0,0.4,nbins);

for d = 1:length(sess_data)
    
   num_lf8_ups = length(rlid_lup{d});
   
   %mp - lfp
   t_10_10_delay{d} = zeros(num_lf8_ups,1);
   t_50_10_delay{d} = zeros(num_lf8_ups,1);
   t_90_10_delay{d} = zeros(num_lf8_ups,1);
   t_10_50_delay{d} = zeros(num_lf8_ups,1);
   t_50_50_delay{d} = zeros(num_lf8_ups,1);
   t_90_50_delay{d} = zeros(num_lf8_ups,1);
   t_10_90_delay{d} = zeros(num_lf8_ups,1);
   t_50_90_delay{d} = zeros(num_lf8_ups,1);
   t_90_90_delay{d} = zeros(num_lf8_ups,1);
  
   
   for i = 1:num_lf8_ups
      
       [dummy, nearest_mp_up] = min(abs(rlid_lup{d}(i) - rlid_wup{d}));
       t_10_10_delay{d}(i) = t_10_wup{d}(nearest_mp_up) - t_10_lup{d}(i);
       t_50_10_delay{d}(i) = rlid_wup{d}(nearest_mp_up) - t_10_lup{d}(i);
       t_90_10_delay{d}(i) = t_90_wup{d}(nearest_mp_up) - t_10_lup{d}(i);
       t_10_50_delay{d}(i) = t_10_wup{d}(nearest_mp_up) - rlid_lup{d}(i);
       t_50_50_delay{d}(i) = rlid_wup{d}(nearest_mp_up) - rlid_lup{d}(i);
       t_90_50_delay{d}(i) = t_90_wup{d}(nearest_mp_up) - rlid_lup{d}(i);
       t_10_90_delay{d}(i) = t_10_wup{d}(nearest_mp_up) - t_90_lup{d}(i);
       t_50_90_delay{d}(i) = rlid_wup{d}(nearest_mp_up) - t_90_lup{d}(i);
       t_90_90_delay{d}(i) = t_90_wup{d}(nearest_mp_up) - t_90_lup{d}(i);
       
   end
    
   t_10_10_delay{d} = t_10_10_delay{d}/Fs;
   t_50_10_delay{d} = t_50_10_delay{d}/Fs;
   t_90_10_delay{d} = t_90_10_delay{d}/Fs;
   t_10_50_delay{d} = t_10_50_delay{d}/Fs;
   t_50_50_delay{d} = t_50_50_delay{d}/Fs;
   t_90_50_delay{d} = t_90_50_delay{d}/Fs;
   t_10_90_delay{d} = t_10_90_delay{d}/Fs;
   t_50_90_delay{d} = t_50_90_delay{d}/Fs;
   t_90_90_delay{d} = t_90_90_delay{d}/Fs;
   
   t_10_10_hist(d,:) = histc(t_10_10_delay{d},lag_range);
   t_50_10_hist(d,:) = histc(t_50_10_delay{d},lag_range);
   t_90_10_hist(d,:) = histc(t_90_10_delay{d},lag_range);
   t_10_50_hist(d,:) = histc(t_10_50_delay{d},lag_range);
   t_50_50_hist(d,:) = histc(t_50_50_delay{d},lag_range);
   t_90_50_hist(d,:) = histc(t_90_50_delay{d},lag_range);
   t_10_90_hist(d,:) = histc(t_10_90_delay{d},lag_range);
   t_50_90_hist(d,:) = histc(t_50_90_delay{d},lag_range);
   t_90_90_hist(d,:) = histc(t_90_90_delay{d},lag_range);
   
   t_10_10_hist(d,:) = t_10_10_hist(d,:)/sum(t_10_10_hist(d,:));
   t_50_10_hist(d,:) = t_50_10_hist(d,:)/sum(t_50_10_hist(d,:));
   t_90_10_hist(d,:) = t_90_10_hist(d,:)/sum(t_90_10_hist(d,:));
   t_10_50_hist(d,:) = t_10_50_hist(d,:)/sum(t_10_50_hist(d,:));
   t_50_50_hist(d,:) = t_50_50_hist(d,:)/sum(t_50_50_hist(d,:));
   t_90_50_hist(d,:) = t_90_50_hist(d,:)/sum(t_90_50_hist(d,:));
   t_10_90_hist(d,:) = t_10_90_hist(d,:)/sum(t_10_90_hist(d,:));
   t_50_90_hist(d,:) = t_50_90_hist(d,:)/sum(t_50_90_hist(d,:));
   t_90_90_hist(d,:) = t_90_90_hist(d,:)/sum(t_90_90_hist(d,:));

   t_10_10_cumdist(d,:) = cumsum(t_10_10_hist(d,:));
   t_50_10_cumdist(d,:) = cumsum(t_50_10_hist(d,:));
   t_90_10_cumdist(d,:) = cumsum(t_90_10_hist(d,:));
   t_10_50_cumdist(d,:) = cumsum(t_10_50_hist(d,:));
   t_50_50_cumdist(d,:) = cumsum(t_50_50_hist(d,:));
   t_90_50_cumdist(d,:) = cumsum(t_90_50_hist(d,:));
   t_10_90_cumdist(d,:) = cumsum(t_10_90_hist(d,:));
   t_50_90_cumdist(d,:) = cumsum(t_50_90_hist(d,:));
   t_90_90_cumdist(d,:) = cumsum(t_90_90_hist(d,:));

   
   mean_10_10_delay(d) = nanmean(t_10_10_delay{d});
   mean_50_10_delay(d) = nanmean(t_50_10_delay{d});
   mean_90_10_delay(d) = nanmean(t_90_10_delay{d});
   mean_10_50_delay(d) = nanmean(t_10_50_delay{d});
   mean_50_50_delay(d) = nanmean(t_50_50_delay{d});
   mean_90_50_delay(d) = nanmean(t_90_50_delay{d});
   mean_10_90_delay(d) = nanmean(t_10_90_delay{d});
   mean_50_90_delay(d) = nanmean(t_50_90_delay{d});
   mean_90_90_delay(d) = nanmean(t_90_90_delay{d});
   
   median_10_10_delay(d) = nanmedian(t_10_10_delay{d});
   median_50_10_delay(d) = nanmedian(t_50_10_delay{d});
   median_90_10_delay(d) = nanmedian(t_90_10_delay{d});
   median_10_50_delay(d) = nanmedian(t_10_50_delay{d});
   median_50_50_delay(d) = nanmedian(t_50_50_delay{d});
   median_90_50_delay(d) = nanmedian(t_90_50_delay{d});
   median_10_90_delay(d) = nanmedian(t_10_90_delay{d});
   median_50_90_delay(d) = nanmedian(t_50_90_delay{d});
   median_90_90_delay(d) = nanmedian(t_90_90_delay{d});

   std_10_10_delay(d) = nanstd(t_10_10_delay{d});
   std_50_10_delay(d) = nanstd(t_50_10_delay{d});
   std_90_10_delay(d) = nanstd(t_90_10_delay{d});
   std_10_50_delay(d) = nanstd(t_10_50_delay{d});
   std_50_50_delay(d) = nanstd(t_50_50_delay{d});
   std_90_50_delay(d) = nanstd(t_90_50_delay{d});
   std_10_90_delay(d) = nanstd(t_10_90_delay{d});
   std_50_90_delay(d) = nanstd(t_50_90_delay{d});
   std_90_90_delay(d) = nanstd(t_90_90_delay{d});

   iqr_10_10_delay(d) = iqr(t_10_10_delay{d});
   iqr_50_10_delay(d) = iqr(t_50_10_delay{d});
   iqr_90_10_delay(d) = iqr(t_90_10_delay{d});
   iqr_10_50_delay(d) = iqr(t_10_50_delay{d});
   iqr_50_50_delay(d) = iqr(t_50_50_delay{d});
   iqr_90_50_delay(d) = iqr(t_90_50_delay{d});
   iqr_10_90_delay(d) = iqr(t_10_90_delay{d});
   iqr_50_90_delay(d) = iqr(t_50_90_delay{d});
   iqr_90_90_delay(d) = iqr(t_90_90_delay{d});
  
   %    plot(lag_range,t_10_hist(d,:))
%    hold on
%    plot(lag_range,t_50_hist(d,:),'r')
%    plot(lag_range,t_90_hist(d,:),'k')
%    plot(lag_range,t_90_10_hist(d,:),'c')
%    legend('10','50','90','90-10')
%    grid
%   cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
%     t_names = ['C:\WC_Germany\Cortical_analysis\sigmoid_fit\rel_delay_dist_' cell_name];
%    print('-dpng',t_names)
%    close all
   
   rltau_wup{d} = rltau_wup{d}/Fs;
   rltau_lup{d} = rltau_lup{d}/Fs;
   
   sigwidth_wup = (t_90_wup{d} - t_10_wup{d})/Fs;
    sigwidth_lup = (t_90_lup{d} - t_10_lup{d})/Fs;
 
%    useable_widths = find(~isnan(sigwidth_wup(nearest_mp_up)) & ~isnan(sigwidth_lup));
%    temp = corrcoef(sigwidth_wup(useable_widths),sigwidth_lup(nearest_lf8_up(useable_widths)));
%    width_corr(d) = temp(2,1);
   
   clear nearest_mp_up
   
%    w_slope_hist(d,:) = jmm_smooth_1d_cor(histc(rltau_wup{d},slope_range),2);
%    l_slope_hist(d,:) = jmm_smooth_1d_cor(histc(rltau_lup{d},slope_range),2);
   
w_sigwidth_hist(d,:) = histc(sigwidth_wup,sig_width_range);
l_sigwidth_hist(d,:) = histc(sigwidth_lup,sig_width_range);
w_sigwidth_cumdist(d,:) = cumsum(w_sigwidth_hist(d,:)/sum(w_sigwidth_hist(d,:)));
l_sigwidth_cumdist(d,:) = cumsum(l_sigwidth_hist(d,:)/sum(l_sigwidth_hist(d,:)));

%    plot(slope_range,w_slope_hist(d,:))
%    hold on
%    plot(slope_range,l_slope_hist(d,:),'r')
%    grid
%    legend('MP','LF8')
%               cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\sigmoid_fit\slope_dist_' cell_name];
%    print('-dpng',t_names)
%    close all
% 
%       plot(sig_width_range,w_sigwidth_hist(d,:))
%    hold on
%    plot(sig_width_range,l_sigwidth_hist(d,:),'r')
%    grid
%    legend('MP','LF8')
%               cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\sigmoid_fit\sigwidth_dist_' cell_name];
%    print('-dpng',t_names)
%    close all

   
%    mean_w_slope(d) = nanmean(rltau_wup{d});
%    mean_l_slope(d) = nanmean(rltau_lup{d});
%    median_w_slope(d) = nanmedian(rltau_wup{d});
%    median_l_slope(d) = nanmedian(rltau_lup{d});

  
   
end

save C:\WC_Germany\Cortical_analysis\sigmoid_fit\sigfitcompare *range *cumdist *hist mean* std* median* iqr*