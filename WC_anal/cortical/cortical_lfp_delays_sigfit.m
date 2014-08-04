clear all
close all

cd C:\WC_Germany\Cortical_analysis
load cortical_dir
load C:\WC_Germany\Cortical_analysis\sigmoid_fit\sig_fit_alldata

Fs = 2016;
nbins = 100;
lag_range = linspace(-0.1,0.1,nbins);
slope_range = linspace(0,100,nbins);

for d = 1:length(sess_data)
    
   num_mp_ups = length(rlid_wup{d});
   
%    t_10_delay87{d} = zeros(num_lf8_ups,1);
%    t_50_delay87{d} = zeros(num_lf8_ups,1);
%    t_90_delay87{d} = zeros(num_lf8_ups,1);
%    t_10_delay85{d} = zeros(num_lf8_ups,1);
%    t_50_delay85{d} = zeros(num_lf8_ups,1);
%    t_90_delay85{d} = zeros(num_lf8_ups,1);
   
t_90_10_delay{d} = zeros(num_mp_ups,1);

   for i = 1:num_mp_ups
      
%        [dummy, nearest_lf7_up] = min(abs(rlid_lup{d}(i) - rlid_lup7{d}));
%        t_10_delay87{d}(i) = t_10_lup7{d}(nearest_lf7_up) - t_10_lup{d}(i);
%        t_50_delay87{d}(i) = rlid_lup7{d}(nearest_lf7_up) - rlid_lup{d}(i);
%        t_90_delay87{d}(i) = t_90_lup7{d}(nearest_lf7_up) - t_90_lup{d}(i);
%        [dummy, nearest_lf5_up] = min(abs(rlid_lup{d}(i) - rlid_lup5{d}));
%        t_10_delay85{d}(i) = t_10_lup5{d}(nearest_lf5_up) - t_10_lup{d}(i);
%        t_50_delay85{d}(i) = rlid_lup5{d}(nearest_lf5_up) - rlid_lup{d}(i);
%        t_90_delay85{d}(i) = t_90_lup5{d}(nearest_lf5_up) - t_90_lup{d}(i);

        t_90_10_delay{d}(i) = t_90_wup{d}(i) - t_10_wup{d}(i);


   end
    
%    t_10_delay87{d} = t_10_delay87{d}/Fs;
%    t_50_delay87{d} = t_50_delay87{d}/Fs;
%    t_90_delay87{d} = t_90_delay87{d}/Fs;
%    t_10_delay85{d} = t_10_delay85{d}/Fs;
%    t_50_delay85{d} = t_50_delay85{d}/Fs;
%    t_90_delay85{d} = t_90_delay85{d}/Fs;
t_90_10_delay{d} = t_90_10_delay{d}/Fs;  

median_90_10_delay(d) = nanmedian(t_90_10_delay{d});
mean_90_10_delay(d) = nanmean(t_90_10_delay{d});

%    t_10_hist87(d,:) = jmm_smooth_1d_cor(histc(t_10_delay87{d},lag_range),2);
%    t_50_hist87(d,:) = jmm_smooth_1d_cor(histc(t_50_delay87{d},lag_range),2);
%    t_90_hist87(d,:) = jmm_smooth_1d_cor(histc(t_90_delay87{d},lag_range),2);
%     t_10_hist85(d,:) = jmm_smooth_1d_cor(histc(t_10_delay85{d},lag_range),2);
%    t_50_hist85(d,:) = jmm_smooth_1d_cor(histc(t_50_delay85{d},lag_range),2);
%    t_90_hist85(d,:) = jmm_smooth_1d_cor(histc(t_90_delay85{d},lag_range),2);
%   
%    mean_10_delay87(d) = nanmean(t_10_delay87{d});
%    mean_50_delay87(d) = nanmean(t_50_delay87{d});
%    mean_90_delay87(d) = nanmean(t_90_delay87{d});
%    mean_10_delay85(d) = nanmean(t_10_delay85{d});
%    mean_50_delay85(d) = nanmean(t_50_delay85{d});
%    mean_90_delay85(d) = nanmean(t_90_delay85{d});
%    
%    median_10_delay87(d) = nanmedian(t_10_delay87{d});
%    median_50_delay87(d) = nanmedian(t_50_delay87{d});
%    median_90_delay87(d) = nanmedian(t_90_delay87{d});
%     median_10_delay85(d) = nanmedian(t_10_delay85{d});
%    median_50_delay85(d) = nanmedian(t_50_delay85{d});
%    median_90_delay85(d) = nanmedian(t_90_delay85{d});
  
%    stairs(lag_range,t_10_hist)
%    hold on
%    stairs(lag_range,t_50_hist,'r')
%    stairs(lag_range,t_90_hist,'k')
%    legend('10','50','90')
%    grid

%    plot(lag_range,t_10_hist87(d,:))
%    hold on
%    plot(lag_range,t_50_hist87(d,:),'r')
%    plot(lag_range,t_90_hist87(d,:),'k')
%       plot(lag_range,t_10_hist85(d,:),'--')
%       plot(lag_range,t_50_hist85(d,:),'r--')
%    plot(lag_range,t_90_hist85(d,:),'k--')
% 
%    legend('10_87','50_87','90_87')
%    grid
%    t_names = ['C:\WC_Germany\Cortical_analysis\sigmoid_fit\lfp_rel_delay_dist_' sess_data(d).name];
%    print('-dpng',t_names)
%    close all
   
%    rltau_lup{d} = rltau_lup{d}/Fs;
%    rltau_lup7{d} = rltau_lup7{d}/Fs;
%    rltau_lup5{d} = rltau_lup5{d}/Fs;
%   
%    l_slope_hist(d,:) = jmm_smooth_1d_cor(histc(rltau_lup{d},slope_range),2);
%    l_slope_hist7(d,:) = jmm_smooth_1d_cor(histc(rltau_lup7{d},slope_range),2);
%    l_slope_hist5(d,:) = jmm_smooth_1d_cor(histc(rltau_lup5{d},slope_range),2);
   
%    plot(slope_range,l_slope_hist(d,:))
%    hold on
%    plot(slope_range,l_slope_hist7(d,:),'r')
%       plot(slope_range,l_slope_hist5(d,:),'k')
%    grid
%    legend('LF8','LF7','LF5')
%    t_names = ['C:\WC_Germany\Cortical_analysis\sigmoid_fit\lfp_slope_dist_' sess_data(d).name];
%    print('-dpng',t_names)
%    close all

   
end


l_23_pyr_par = 1:11;

int_par = 12:14;

l_56_pyr_par = 15:24;

l_23_pyr_fro = 25:28;

l_5_pyr_fro = 29:31;

l_23_pyr_pre = 32:36;

l_5_pyr_pre = 37:39;

l_23_pyr_fro = [l_23_pyr_fro l_23_pyr_pre];
l_5_pyr_fro = [l_5_pyr_fro l_5_pyr_pre];

thal = 40:43;

barrel = 44:47;



hist_range = linspace(-0.3, 0.3, 30);

hmed_10_23_par = histc(median_10_delay(l_23_pyr_par),hist_range);
hmed_10_56_par = histc(median_10_delay(l_56_pyr_par),hist_range);
hmed_10_23_fro = histc(median_10_delay(l_23_pyr_fro),hist_range);
hmed_10_5_fro = histc(median_10_delay(l_5_pyr_fro), hist_range);
hmed_10_23_int = histc(median_10_delay(int_par),hist_range);
hmed_10_thal = histc(median_10_delay(thal),hist_range);
hmed_10_barrel = histc(median_10_delay(barrel), hist_range);
hmed_10_23_pre = histc(median_10_delay(l_23_pyr_pre), hist_range);


stairs(hist_range, hmed_10_23_par)
hold on
stairs(hist_range, hmed_10_56_par,'r')
stairs(hist_range, hmed_10_23_fro,'k')
stairs(hist_range, hmed_10_5_fro,'g')
stairs(hist_range, hmed_10_barrel, 'c')
stairs(hist_range, hmed_10_thal,'Color',[0.2 0.3 0.4])
% legend('23 par','56 par','23 fro','5 fro','barrel','thal')
legend('23 par','56 par','23 fro','5 fro')



hmed_50_23_par = histc(median_50_delay(l_23_pyr_par),hist_range);
hmed_50_56_par = histc(median_50_delay(l_56_pyr_par),hist_range);
hmed_50_23_fro = histc(median_50_delay(l_23_pyr_fro),hist_range);
hmed_50_5_fro = histc(median_50_delay(l_5_pyr_fro), hist_range);
hmed_50_23_int = histc(median_50_delay(int_par),hist_range);
hmed_50_thal = histc(median_50_delay(thal),hist_range);
hmed_50_barrel = histc(median_50_delay(barrel), hist_range);
hmed_50_23_pre = histc(median_50_delay(l_23_pyr_pre), hist_range);


stairs(hist_range, hmed_50_23_par)
hold on
stairs(hist_range, hmed_50_56_par,'r')
stairs(hist_range, hmed_50_23_fro,'k')
stairs(hist_range, hmed_50_5_fro,'g')
stairs(hist_range, hmed_50_barrel, 'c')
stairs(hist_range, hmed_50_thal,'Color',[0.2 0.3 0.4])
% legend('23 par','56 par','23 fro','5 fro','barrel','thal')
legend('23 par','56 par','23 fro','5 fro')



hmed_90_23_par = histc(median_90_delay(l_23_pyr_par),hist_range);
hmed_90_56_par = histc(median_90_delay(l_56_pyr_par),hist_range);
hmed_90_23_fro = histc(median_90_delay(l_23_pyr_fro),hist_range);
hmed_90_5_fro = histc(median_90_delay(l_5_pyr_fro), hist_range);
hmed_90_23_int = histc(median_90_delay(int_par),hist_range);
hmed_90_thal = histc(median_90_delay(thal),hist_range);
hmed_90_barrel = histc(median_90_delay(barrel), hist_range);
hmed_90_23_pre = histc(median_90_delay(l_23_pyr_pre), hist_range);


stairs(hist_range, hmed_90_23_par)
hold on
stairs(hist_range, hmed_90_56_par,'r')
stairs(hist_range, hmed_90_23_fro,'k')
stairs(hist_range, hmed_90_5_fro,'g')
stairs(hist_range, hmed_90_barrel, 'c')
stairs(hist_range, hmed_90_thal,'Color',[0.2 0.3 0.4])
% legend('23 par','56 par','23 fro','5 fro','barrel','thal')
legend('23 par','56 par','23 fro','5 fro')
