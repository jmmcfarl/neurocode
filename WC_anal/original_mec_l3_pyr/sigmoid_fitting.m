%% Fit sigmoid to transition points

clear all

load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28

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
   wcv(wcv > 3) = 3;
   
   wcv_filt = filtfilt(b2,a2,wcv_minus_spike);
   wcv_slope = [0;diff(wcv_filt)];
   wcv_slope = wcv_slope/std(wcv_slope);
   
   lf8 = filtfilt(b,a,lf8);
   lf8 = zscore(lf8);
   lf8 = lf8';
   lf8_filt = filtfilt(b2,a2,lf8);
   lf8_slope = [0 diff(lf8_filt)];
   lf8_slope = lf8_slope/std(lf8_slope);

   cur_wcv_up_trans = up_trans{d}*dsf;
   cur_wcv_down_trans = down_trans{d}*dsf;
   cur_lf8_up_trans = up_trans8{d}*dsf;
   cur_lf8_down_trans = down_trans8{d}*dsf;
   
%% fit wcv state transitions
   [rlid_wup{d},rltime_wup{d},rlamp_wup{d},rlshift_wup{d},rltau_wup{d},rlerror_wup{d},t_90_wup{d},t_10_wup{d}] = ...
       get_lfp_wcv_sigmoid_fit_ut(cur_wcv_up_trans,wcv,wcv_slope,synct);
%      
%   [rlid_wdown{d},rltime_wdown{d},rlamp_wdown{d},rlshift_wdown{d},rltau_wdown{d},rlerror_wdown{d},t_10_wdown{d}] = ...
%       get_lfp_wcv_sigmoid_fit_dt(cur_wcv_down_trans,-wcv,-wcv_slope,synct);

%% fit lf8 state transitions
   [rlid_lup{d},rltime_lup{d},rlamp_lup{d},rlshift_lup{d},rltau_lup{d},rlerror_lup{d},t_90_lup{d},t_10_lup{d}] = ...
       get_lfp_wcv_sigmoid_fit_ut(cur_lf8_up_trans,lf8,lf8_slope,synct);
     
%   [rlid_ldown{d},rltime_ldown{d},rlamp_ldown{d},rlshift_ldown{d},rltau_ldown{d},rlerror_ldown{d},t_10_ldown{d}] = ...
%       get_lfp_wcv_sigmoid_fit_dt(cur_lf8_down_trans,-lf8,-lf8_slope,synct);


%% eliminate bad up transitions
rltau_wup{d}(rltau_wup{d} <= 0) = nan;
rltau_lup{d}(rltau_lup{d} <= 0) = nan;
t_90_wup{d}(rltau_wup{d} <=0) = nan;
t_90_lup{d}(rltau_lup{d} <=0) = nan;
t_10_wup{d}(rltau_wup{d} <=0) = nan;
t_10_lup{d}(rltau_lup{d} <=0) = nan;

z_w = rltau_wup{d}-nanmean(rltau_wup{d});
z_w = z_w/nanstd(z_w);
z_l = rltau_lup{d}-nanmean(rltau_lup{d});
z_l = z_l/nanstd(z_l);

rltau_wup{d}(abs(z_w) > 3) = nan;
rltau_lup{d}(abs(z_l) > 3) = nan;
t_90_wup{d}(abs(z_w) > 3) = nan;
t_90_lup{d}(abs(z_l) > 3) = nan;
t_10_wup{d}(abs(z_w) > 3) = nan;
t_10_lup{d}(abs(z_l) > 3) = nan;

  %eliminate bad up fits
%    z_errors = zscore(rlerror_wup{d});
   
%    bad_errors = find(z_errors > 2);
   
%    rlid_wup{d}(bad_errors) = nan;
%    rltime_wup{d}(bad_errors) = nan;
%    rlamp_wup{d}(bad_errors) = nan;
%    rlshift_wup{d}(bad_errors) = nan;
%    rltau_wup{d}(bad_errors) = nan;
%    rlerror_wup{d}(bad_errors) = nan;  
%    t_90_wup{d}(bad_errors) = nan;
   
%    %eliminate bad down fits
%        z_errors = zscore(rlerror_wdown{d});
%    
%    bad_errors = find(z_errors > 2);
%    
%    rlid_wdown{d}(bad_errors) = nan;
%    rltime_wdown{d}(bad_errors) = nan;
%    rlamp_wdown{d}(bad_errors) = nan;
%    rlshift_wdown{d}(bad_errors) = nan;
%    rltau_wdown{d}(bad_errors) = nan;
%    rlerror_wdown{d}(bad_errors) = nan;  
%     t_90_wdown{d}(bad_errors) = nan;
    
%find corrected times and 90% times



   clear synct
   
end


%% Calculate Slope statistics
% slope_axis = linspace(0,2e5,200);
% 
% for d = 1:17
%     
%    %get rid of instances of negative slopes
%    neg_slopes = find(rltau_wup{d} < 0);
%    rltau_wup{d}(neg_slopes) = nan;
%    neg_slopes = find(rltau_wdown{d} < 0);
%    rltau_wdown{d}(neg_slopes) = nan;
%    
%    mean_up_slope(d) = nanmean(rltau_wup{d});
%    median_up_slope(d) = nanmedian(rltau_wup{d});
%    mean_down_slope(d) = nanmean(rltau_wdown{d});
%    median_down_slope(d) = nanmedian(rltau_wdown{d});
%    
%    %estimate density of slopes
%    dens_est(d,:) = ksdensity(rltau_wup{d},slope_axis);
%    [max_val,max_loc] = max(dens_est(d,:));
%    mode_up_slope(d) = slope_axis(max_loc);
%     
%    dens_est(d,:) = ksdensity(rltau_wdown{d},slope_axis);
%    [max_val,max_loc] = max(dens_est(d,:));
%    mode_down_slope(d) = slope_axis(max_loc);
% 
%     
% end
% clear dens_est
% %% Calculate state duration statistics
% 
% dur_axis = linspace(0,20,200);
% 
% for d = 1:17
%        
%    mean_up_dur(d) = nanmean(up_state_dur{d});
%    median_up_dur(d) = nanmedian(up_state_dur{d});
%    
%    mean_down_dur(d) = nanmean(down_state_dur{d});
%    median_down_dur(d) = nanmedian(down_state_dur{d});
%    
%    
%    %estimate density of slopes
%    dens_est(d,:) = ksdensity(up_state_dur{d},dur_axis);
%    [max_val,max_loc] = max(dens_est(d,:));
%    mode_up_dur(d) = dur_axis(max_loc);
%     
%    dens_est(d,:) = ksdensity(down_state_dur{d},dur_axis);
%    [max_val,max_loc] = max(dens_est(d,:));
%    mode_down_dur(d) = dur_axis(max_loc);
%     
% end
