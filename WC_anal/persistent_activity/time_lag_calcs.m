clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

Fsd = 2016/8;
up_lag_range = [0:0.02:1];
down_lag_range = [0:0.04:2];
near_up_lag_range = [-0.5:0.03:1];

for d = 1:length(dir_array)
    
    up_lag{d} = zeros(size(synch_ups{d}));
    down_lag{d} = zeros(size(synch_ups{d}));
    
   for i = 1:length(synch_ups{d})
       
       prev_up = find(up_trans8{d} < up_trans{d}(synch_ups{d}(i)),1,'last');
       [dummy,near_up_loc] = min(abs(up_trans8{d}-up_trans{d}(synch_ups{d}(i))));
       prev_down = find(down_trans8{d} < down_trans{d}(synch_ups{d}(i)),1,'last');
       [dummy,near_down_loc] = min(abs(down_trans8{d}-down_trans{d}(synch_ups{d}(i))));
       
       if ~isempty(prev_up)
           up_lag_prev{d}(i) = (up_trans{d}(synch_ups{d}(i))-up_trans8{d}(prev_up))/Fsd;
           up_lag_near{d}(i) = (up_trans{d}(synch_ups{d}(i))-up_trans8{d}(near_up_loc))/Fsd;
       else
           up_lag_prev{d}(i) = nan;
           up_lag_near{d}(i) = nan;
       end
       
       if ~isempty(prev_down)
           down_lag{d}(i) = (down_trans{d}(synch_ups{d}(i))-down_trans8{d}(prev_down))/Fsd;
           down_lag_near{d}(i) = (down_trans{d}(synch_ups{d}(i)) - down_trans8{d}(near_down_loc))/Fsd;
       else
           down_lag{d}(i) = nan;
           down_lag_near{d}(i) = nan;
       end
       
   end
    
   up_lag_median(d) = nanmedian(up_lag_prev{d});
   up_lag_iqr(d) = iqr(up_lag_prev{d});
   down_lag_median(d) = nanmedian(down_lag{d});
   down_lag_iqr(d) = iqr(down_lag{d});
   up_lag_mean(d) = nanmean(up_lag_prev{d});
   up_lag_std(d) = nanstd(up_lag_prev{d});
   down_lag_mean(d) = nanmean(down_lag{d});
   down_lag_std(d) = nanstd(down_lag{d});
      
   up_near_median(d) = nanmedian(up_lag_near{d});
   up_near_iqr(d) = iqr(up_lag_near{d});
   down_near_median(d) = nanmedian(down_lag_near{d});
   down_near_iqr(d) = iqr(down_lag_near{d});
   up_near_mean(d) = nanmean(up_lag_near{d});
   up_near_std(d) = nanstd(up_lag_near{d});
   down_near_mean(d) = nanmean(down_lag_near{d});
   down_near_std(d) = nanstd(down_lag_near{d});

%    up_lag_prev_dist(d,:) = hist(up_lag_prev{d},up_lag_range);
%    up_lag_prev_dist(d,:) = up_lag_prev_dist(d,:)/sum(up_lag_prev_dist(d,:));
%       up_lag_near_dist(d,:) = hist(up_lag_near{d},near_up_lag_range);
%    up_lag_near_dist(d,:) = up_lag_near_dist(d,:)/sum(up_lag_near_dist(d,:));
%    down_lag_dist(d,:) = hist(down_lag{d},down_lag_range);
%    down_lag_dist(d,:) = down_lag_dist(d,:)/sum(down_lag_dist(d,:));
%    
%    subplot(2,1,1)
%    bar(up_lag_range,up_lag_prev_dist(d,:))
%    xlim([0 0.9])
%    subplot(2,1,2)
%    bar(near_up_lag_range,up_lag_near_dist(d,:),'r')
%    xlim([-0.45 .9])
%    t_names = ['C:\WC_Germany\Persistent_activity\time_lags\up_' f_names{d}];
%    print('-dpng',t_names)
%    close all
%    
%    bar(down_lag_range,down_lag_dist(d,:))
%    xlim([0 1.9])
%    t_names = ['C:\WC_Germany\Persistent_activity\time_lags\down_' f_names{d}];
%    print('-dpng',t_names)
%    close all

end