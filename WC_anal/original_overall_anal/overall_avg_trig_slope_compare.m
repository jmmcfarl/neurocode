 
clear all
close all
 

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_data


%% cycle through all cell types and compute the average trig averages
diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));  

    for i = 1:length(get_cells)
        
        midpt = (min(mp_utrig_mp(get_cells(i),:))+max(mp_utrig_mp(get_cells(i),:)))/2;
        mp_utrig_mp(get_cells(i),:) = mp_utrig_mp(get_cells(i),:)-midpt;
        
        midpt = (min(mp_dtrig_mp(get_cells(i),:))+max(mp_dtrig_mp(get_cells(i),:)))/2;
        mp_dtrig_mp(get_cells(i),:) = mp_dtrig_mp(get_cells(i),:)-midpt;
        
        midpt = (min(lf8_utrig_lf8(get_cells(i),:))+max(lf8_utrig_lf8(get_cells(i),:)))/2;
        lf8_utrig_lf8(get_cells(i),:) = lf8_utrig_lf8(get_cells(i),:)-midpt;
        
        midpt = (min(lf8_dtrig_lf8(get_cells(i),:))+max(lf8_dtrig_lf8(get_cells(i),:)))/2;
        lf8_dtrig_lf8(get_cells(i),:) = lf8_dtrig_lf8(get_cells(i),:)-midpt;
        
    end
    
    
    plot(lags,mp_utrig_mp(get_cells,:))
    xlim([-0.5 0.5]);grid
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\mp_up_slope_compare_' num2str(c)];
    print('-dpng',t_names);
    close 
    
    plot(lags,mp_dtrig_mp(get_cells,:))
    xlim([-0.5 0.5]);grid
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\mp_down_slope_compare_' num2str(c)];
    print('-dpng',t_names);
    close
    
    plot(lags,lf8_utrig_lf8(get_cells,:))
    xlim([-0.5 0.5]);grid
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\lf8_up_slope_compare_' num2str(c)];
    print('-dpng',t_names);
    close
    
    plot(lags,lf8_dtrig_lf8(get_cells,:))
    xlim([-0.5 0.5]);grid
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\lf8_down_slope_compare_' num2str(c)];
    print('-dpng',t_names);
 close

 
end