clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data_new_method
load C:\WC_Germany\Persistent_activity\lf8_period_f_data_new_method
load C:\WC_Germany\Persistent_activity\UDS_synch_state_dur\UDS_synch_state_dur_data_new_method

phase_range = linspace(0,360,100);
clear up_dist down_dist
for d = 1:length(dir_array)
    
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
bad_values = find(isnan(up_phase{d}) | isnan(down_phase{d}));
up_phase{d}(bad_values) = [];
down_phase{d}(bad_values) = [];
    [corr_coef(d),p_value(d)] = circ_corrcc(up_phase{d}*pi/180,down_phase{d}*pi/180)

    scatter(up_phase{d},down_phase{d})
    pause
    clf
    
   up_dist(d,:) = hist(up_phase{d},phase_range);
   up_dist(d,:) = up_dist(d,:)/sum(up_dist(d,:));
   up_dist(d,:) = fgsmooth(up_dist(d,:),2);
   
   down_dist(d,:) = hist(down_phase{d},phase_range);
   down_dist(d,:) = down_dist(d,:)/sum(down_dist(d,:));
   down_dist(d,:) = fgsmooth(down_dist(d,:),2);
   
%    polar(phase_range/360*2*pi,up_dist(d,:))
%    hold on
%    polar(phase_range/360*2*pi,down_dist(d,:),'r')
%    
%       t_names = ['C:\WC_Germany\Persistent_activity\time_lags\phase_' f_names{d}];
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

u_up = mean_up+1*std(up_dist)/sqrt(length(dir_array));
l_up = mean_up-1*std(up_dist)/sqrt(length(dir_array));

u_down = mean_down+1*std(down_dist)/sqrt(length(dir_array));
l_down = mean_down-1*std(down_dist)/sqrt(length(dir_array));
phase_rads = phase_range/360*2*pi;
X = [u_up.*cos(phase_rads) fliplr(l_up.*cos(phase_rads))];
Y = [u_up.*sin(phase_rads) fliplr(l_up.*sin(phase_rads))];
polar(phase_range/360*2*pi,mean_up)
% polar(phase_range/360*2*pi,l_up,'--')
hold on
polar(phase_range/360*2*pi,mean_down,'r')

fill(X,Y,'b')
hold on

X = [u_down.*cos(phase_rads) fliplr(l_down.*cos(phase_rads))];
Y = [u_down.*sin(phase_rads) fliplr(l_down.*sin(phase_rads))];
fill(X,Y,'r')
polar(phase_range/360*2*pi,u_up,'--')
hold on
polar(phase_range/360*2*pi,u_down,'--r')
polar(phase_range/360*2*pi,l_down,'--r')