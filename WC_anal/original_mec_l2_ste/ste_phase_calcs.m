clear all

load C:\WC_Germany\JMM_analysis_ste\dir_tree_ste
load C:\WC_Germany\JMM_analysis_ste\UDS_dur_run_hist\data
load C:\WC_Germany\JMM_analysis_ste\lf8_period_f_data
clear up_dist down_dist
phase_range = [0:3:360];

for d = 1:length(dir_array)
    
    synch_ups{d} = 1:length(up_trans{d});
    
for i = 1:length(synch_ups{d})
    
    if up_trans{d}(synch_ups{d}(i)) > length(lf8_period_p{d})
        up_phase{d}(i) = nan;
    else
        up_phase{d}(i) = lf8_period_p{d}(up_trans{d}(synch_ups{d}(i)));
    end
    if down_trans{d}(synch_ups{d}(i)) > length(lf8_period_p{d})
        down_phase{d}(i) = nan;
    else
        down_phase{d}(i) = lf8_period_p{d}(down_trans{d}(synch_ups{d}(i)));
    end
    
end 

   up_dist(d,:) = hist(up_phase{d},phase_range);
   up_dist(d,:) = up_dist(d,:)/sum(up_dist(d,:));
   
   down_dist(d,:) = hist(down_phase{d},phase_range);
   down_dist(d,:) = down_dist(d,:)/sum(down_dist(d,:));

%    polar(phase_range/360*2*pi,up_dist(d,:))
%    hold on
%    polar(phase_range/360*2*pi,down_dist(d,:),'r')
%    
%       t_names = ['C:\WC_Germany\JMM_analysis_ste\time_lags\phase_' f_names{d}];
%    print('-dpng',t_names)
%    close all

   up_phase{d}(isnan(up_phase{d})) = [];
   down_phase{d}(isnan(down_phase{d})) = [];
   circ_up_mean(d) = circ_mean(up_phase{d}/360*2*pi);
     circ_down_mean(d) = circ_mean(down_phase{d}/360*2*pi);
 circ_up_var(d) = circ_var(up_phase{d}/360*2*pi);
 circ_down_var(d) = circ_var(down_phase{d}/360*2*pi);
  
   
end

mean_up = mean(up_dist);
mean_down = mean(down_dist);

u_up = mean_up+2*std(up_dist)/sqrt(length(dir_array));
l_up = mean_up-2*std(up_dist)/sqrt(length(dir_array));

u_down = mean_down+2*std(down_dist)/sqrt(length(dir_array));
l_down = mean_down-2*std(down_dist)/sqrt(length(dir_array));

polar(phase_range/360*2*pi,u_up,'--')
hold on
polar(phase_range/360*2*pi,mean_up)
polar(phase_range/360*2*pi,l_up,'--')

polar(phase_range/360*2*pi,mean_down,'r')
polar(phase_range/360*2*pi,u_down,'--r')
polar(phase_range/360*2*pi,l_down,'--r')