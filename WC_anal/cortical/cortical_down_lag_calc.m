clear all
close all

cd C:\WC_Germany\Cortical_analysis
load cortical_dir
load C:\WC_Germany\Cortical_analysis\trig_avgs\trig_avg_data_wideband
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data

Fs = 2016;
dsf = 8;
nbins = 100;
lag_range = linspace(-0.6,0.6,nbins);
Fsd = Fs/dsf;

for d = 1:length(sess_data)
    
   num_mp_ups = length(synch_ups{d});
   
   mp_down_lag{d} = zeros(1,num_mp_ups);
   
   for i = 1:num_mp_ups
      
       [dummy, nearest_lf8_down] = min(abs(down_trans{d}(synch_ups{d}(i)) - down_trans8{d}));
        mp_down_lag{d}(i) = down_trans8{d}(nearest_lf8_down) - down_trans{d}(synch_ups{d}(i));
        
   end
    
    mp_down_lag{d} = mp_down_lag{d}/Fsd;
   
    mp_down_lag_hist(d,:) = histc(mp_down_lag{d},lag_range);
    mp_down_lag_hist(d,:) = mp_down_lag_hist(d,:)/sum(mp_down_lag_hist(d,:));
    
   plot(lag_range,mp_down_lag_hist(d,:))
   grid
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    t_names = ['C:\WC_Germany\Cortical_analysis\down_lag\down_lag_dist_' cell_name];
   print('-dpng',t_names)
   close all
   
   
end

save C:\WC_Germany\Cortical_analysis\down_lag\down_lag_dist_data lag_range mp_down_lag_hist