clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\lag_phase_data

lag_range = linspace(-2,2,200);

for d = 1:length(over_dir)


    up_dist(d,:) = hist(lfp_up_lag_near{d},lag_range);
    up_dist(d,:) = up_dist(d,:)/sum(up_dist(d,:));
    up_dist(d,:) = fgsmooth(up_dist(d,:),2);

    down_dist(d,:) = hist(lfp_down_lag_near{d},lag_range);
    down_dist(d,:) = down_dist(d,:)/sum(down_dist(d,:));
    down_dist(d,:) = fgsmooth(down_dist(d,:),2);

%        plot(lag_range,up_dist(d,:))
%        hold on
%        plot(lag_range,down_dist(d,:),'r')
%     
%         t_names = ['C:\WC_Germany\overall_calcs\time_lag\' num2str(cell_type(d)) '_' over_names{d}];
%        print('-dpng',t_names)
%        close all

end


save C:\WC_Germany\overall_calcs\time_lag\time_lag_data up_dist down_dist lag_range