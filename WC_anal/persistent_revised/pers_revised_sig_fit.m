%% Fit sigmoid to transition points

clear all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data

Fs = 2016;
niqf = Fs/2;
dsf = 8;

lcf = 0.05/niqf;
hcf = 40/niqf;
hcf2 = 5/niqf;

[b,a] = butter(2,[lcf hcf]);
[b2,a2] = butter(2,[lcf hcf2]);

% t_axis = -0.25*Fs:0.25*Fs;

for d = 1:length(dir_array)
    
   cd(dir_array{d});
   pwd

   load sync_times
   load used_data wcv_minus_spike lf8
   
   wcv = filtfilt(b,a,wcv_minus_spike);
   wcv = zscore(wcv);
   wcv = wcv';
   
   wcv_filt = filtfilt(b2,a2,wcv_minus_spike);
   wcv_slope = [0;diff(wcv_filt)];
   wcv_slope = wcv_slope/std(wcv_slope);
   
   lf8 = filtfilt(b,a,lf8);
   lf8 = zscore(lf8);
   lf8 = lf8';
   lf8_filt = filtfilt(b2,a2,lf8);
   lf8_slope = [0 diff(lf8_filt)];
   lf8_slope = lf8_slope/std(lf8_slope);

   cur_wcv_up_trans = up_trans{d}(synch_ups{d})*dsf;
   cur_wcv_down_trans = down_trans{d}(synch_downs{d})*dsf;
   cur_lf8_up_trans = up_trans8{d}(synch_ups8{d})*dsf;
   cur_lf8_down_trans = down_trans8{d}(synch_downs8{d})*dsf;
   
%% fit wcv state transitions
   [rlid_wup{d},rltime_wup{d},rlamp_wup{d},rlshift_wup{d},rltau_wup{d},rlerror_wup{d},t_90_wup{d},t_10_wup{d}] = ...
       get_lfp_wcv_sigmoid_fit_ut(cur_wcv_up_trans,wcv,wcv_slope,synct);
%      
  [rlid_wdown{d},rltime_wdown{d},rlamp_wdown{d},rlshift_wdown{d},rltau_wdown{d},rlerror_wdown{d},t_10_wdown{d}] = ...
      get_lfp_wcv_sigmoid_fit_dt(cur_wcv_down_trans,-wcv,-wcv_slope,synct);

%% fit lf8 state transitions
   [rlid_lup{d},rltime_lup{d},rlamp_lup{d},rlshift_lup{d},rltau_lup{d},rlerror_lup{d},t_90_lup{d},t_10_lup{d}] = ...
       get_lfp_wcv_sigmoid_fit_ut(cur_lf8_up_trans,lf8,lf8_slope,synct);
     
  [rlid_ldown{d},rltime_ldown{d},rlamp_ldown{d},rlshift_ldown{d},rltau_ldown{d},rlerror_ldown{d},t_10_ldown{d}] = ...
      get_lfp_wcv_sigmoid_fit_dt(cur_lf8_down_trans,-lf8,-lf8_slope,synct);


%% eliminate bad up transitions
  % eliminate bad wup fits
   z_errors_wup = rlerror_wup{d} - nanmean(rlerror_wup{d});
   z_errors_wup = z_errors_wup./nanstd(z_errors_wup);
   bad_errors_wup = find(z_errors_wup > 2); 
   rlid_wup{d}(bad_errors_wup) = nan;
   rltime_wup{d}(bad_errors_wup) = nan;
   rlamp_wup{d}(bad_errors_wup) = nan;
   rlshift_wup{d}(bad_errors_wup) = nan;
   rltau_wup{d}(bad_errors_wup) = nan;
   rlerror_wup{d}(bad_errors_wup) = nan;  
   t_90_wup{d}(bad_errors_wup) = nan;

   mean_tau_wup(d) = nanmean(rltau_wup{d});
   
     % eliminate bad wdown fits
   z_errors_wdown = rlerror_wdown{d} - nanmean(rlerror_wdown{d});
   z_errors_wdown = z_errors_wdown./nanstd(z_errors_wdown);
   bad_errors_wdown = find(z_errors_wdown > 2); 
   rlid_wdown{d}(bad_errors_wdown) = nan;
   rltime_wdown{d}(bad_errors_wdown) = nan;
   rlamp_wdown{d}(bad_errors_wdown) = nan;
   rlshift_wdown{d}(bad_errors_wdown) = nan;
   rltau_wdown{d}(bad_errors_wdown) = nan;
   rlerror_wdown{d}(bad_errors_wdown) = nan;  
   t_90_wdown{d}(bad_errors_wdown) = nan;

   mean_tau_wdown(d) = nanmean(rltau_wdown{d});
   
   %eliminate bad lup fits
   z_errors_lup = rlerror_lup{d} - nanmean(rlerror_lup{d});
   z_errors_lup = z_errors_lup./nanstd(z_errors_lup);
   bad_errors_lup = find(z_errors_lup > 2); 
   rlid_lup{d}(bad_errors_lup) = nan;
   rltime_lup{d}(bad_errors_lup) = nan;
   rlamp_lup{d}(bad_errors_lup) = nan;
   rlshift_lup{d}(bad_errors_lup) = nan;
   rltau_lup{d}(bad_errors_lup) = nan;
   rlerror_lup{d}(bad_errors_lup) = nan;  
   t_90_lup{d}(bad_errors_lup) = nan;
   
   mean_tau_lup(d) = nanmean(rltau_lup{d});
   
      %eliminate bad lup fits
   z_errors_ldown = rlerror_ldown{d} - nanmean(rlerror_ldown{d});
   z_errors_ldown = z_errors_ldown./nanstd(z_errors_ldown);
   bad_errors_ldown = find(z_errors_ldown > 2); 
   rlid_ldown{d}(bad_errors_ldown) = nan;
   rltime_ldown{d}(bad_errors_ldown) = nan;
   rlamp_ldown{d}(bad_errors_ldown) = nan;
   rlshift_ldown{d}(bad_errors_ldown) = nan;
   rltau_ldown{d}(bad_errors_ldown) = nan;
   rlerror_ldown{d}(bad_errors_ldown) = nan;  
   t_90_ldown{d}(bad_errors_ldown) = nan;

   mean_tau_ldown(d) = nanmean(rltau_ldown{d});
   
end

save C:\WC_Germany\persistent_revised\sig_fit_alldata t_90* t_10* rl*
