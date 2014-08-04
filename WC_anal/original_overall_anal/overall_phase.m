clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\lag_phase_data

phase_range = linspace(0,360,100);

for d = 1:length(over_dir)


    up_dist(d,:) = hist(mp_up_phase{d},phase_range);
    up_dist(d,:) = up_dist(d,:)/sum(up_dist(d,:));
    up_dist(d,:) = fgsmooth(up_dist(d,:),2);

    down_dist(d,:) = hist(mp_down_phase{d},phase_range);
    down_dist(d,:) = down_dist(d,:)/sum(down_dist(d,:));
    down_dist(d,:) = fgsmooth(down_dist(d,:),2);

       polar(phase_range/360*2*pi,up_dist(d,:))
       hold on
       polar(phase_range/360*2*pi,down_dist(d,:),'r')
    
    t_names = ['C:\WC_Germany\overall_calcs\phase\' num2str(cell_type(d)) '_' over_names{d}];
       print('-dpng',t_names)
       close all

    up_phase{d}(isnan(mp_up_phase{d})) = [];
    down_phase{d}(isnan(mp_down_phase{d})) = [];
    circ_up_mean(d) = circ_mean(mp_up_phase{d}/360*2*pi);
    circ_down_mean(d) = circ_mean(mp_down_phase{d}/360*2*pi);
    circ_up_var(d) = circ_var(mp_up_phase{d}/360*2*pi);
    circ_down_var(d) = circ_var(mp_down_phase{d}/360*2*pi);

end


save C:\WC_Germany\overall_calcs\phase\phase_data circ* up_dist down_dist phase_range