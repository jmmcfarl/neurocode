clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\sta\sta_data
Fsd = Fs
diff_cell_types = unique(cell_type);
cmap = colormap(jet(8));
for c = 1:length(diff_cell_types)

    get_cells = find(cell_type == diff_cell_types(c));
    
    avg_sta_wcv(c,:) = nanmean(sta_wcv(get_cells,:));
    uci_sta_wcv = avg_sta_wcv(c,:)+tinv(0.975,length(get_cells))*nanstd(sta_wcv(get_cells,:))/sqrt(length(get_cells));
    lci_sta_wcv = avg_sta_wcv(c,:)-tinv(0.975,length(get_cells))*nanstd(sta_wcv(get_cells,:))/sqrt(length(get_cells));
    
    avg_sta_lf8(c,:) = nanmean(sta_lf8(get_cells,:));
    uci_sta_lf8 = avg_sta_lf8(c,:)+tinv(0.975,length(get_cells))*nanstd(sta_lf8(get_cells,:))/sqrt(length(get_cells));
    lci_sta_lf8 = avg_sta_lf8(c,:)-tinv(0.975,length(get_cells))*nanstd(sta_lf8(get_cells,:))/sqrt(length(get_cells));
    
    avg_sta_lf3(c,:) = nanmean(sta_lf3(get_cells,:));
    uci_sta_lf3 = avg_sta_lf3(c,:)+tinv(0.975,length(get_cells))*nanstd(sta_lf3(get_cells,:))/sqrt(length(get_cells));
    lci_sta_lf3 = avg_sta_lf3(c,:)-tinv(0.975,length(get_cells))*nanstd(sta_lf3(get_cells,:))/sqrt(length(get_cells));

        avg_sta_lf2(c,:) = nanmean(sta_lf2(get_cells,:));
    uci_sta_lf2 = avg_sta_lf2(c,:)+tinv(0.975,length(get_cells))*nanstd(sta_lf2(get_cells,:))/sqrt(length(get_cells));
    lci_sta_lf2 = avg_sta_lf2(c,:)-tinv(0.975,length(get_cells))*nanstd(sta_lf2(get_cells,:))/sqrt(length(get_cells));

    plot(lags/Fsd,avg_sta_wcv(c,:),'linewidth',2)
    hold on
    plot(lags/Fsd,uci_sta_wcv,'--')
    plot(lags/Fsd,lci_sta_wcv,'--')
    xlim([-0.1 0.1]);grid on
    t_names = ['C:\WC_Germany\overall_calcs\sta_compare\sta_wcv_' num2str(c)];
    title(['mean rate ' num2str(mean(mean_rate(get_cells)))])
    print('-dpng',t_names);
    close

        plot(lags/Fsd,avg_sta_lf8(c,:),'r','linewidth',2)
        hold on
    plot(lags/Fsd,avg_sta_lf3(c,:),'g','linewidth',2)
    plot(lags/Fsd,avg_sta_lf2(c,:),'k','linewidth',2)
legend('LF8','LF3','LF2')
    plot(lags/Fsd,uci_sta_lf8,'r--')
    plot(lags/Fsd,lci_sta_lf8,'r--')
        plot(lags/Fsd,uci_sta_lf3,'g--')
    plot(lags/Fsd,lci_sta_lf3,'g--')
    plot(lags/Fsd,uci_sta_lf2,'k--')
    plot(lags/Fsd,lci_sta_lf2,'k--')
    xlim([-0.2 0.2])
    t_names = ['C:\WC_Germany\overall_calcs\sta_compare\sta_lfp_' num2str(c)];
    title(['mean rate ' num2str(mean(mean_rate(get_cells)))])
    print('-dpng',t_names);
    close


        plot(lags/Fsd,avg_sta_lf2(c,:),'color',cmap(c,:),'linewidth',2)
        hold on
        plot(lags/Fsd,uci_sta_lf2,'--','color',cmap(c,:))
        plot(lags/Fsd,lci_sta_lf2,'--','color',cmap(c,:))

    xlim([-0.2 0.2])
end


%% Plot direct comparisons
%MP-LFP
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_sta_wcv(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-0.1 0.1]);grid
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\sta_compare\overall_sta_wcv'];
    print('-dpng',t_names);
    xlim([-0.01 0.01])
        t_names = ['C:\WC_Germany\overall_calcs\sta_compare\overall_sta_wcv_zoom'];
    print('-dpng',t_names);

    close

    for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_sta_lf8(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-0.2 0.2]);grid
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\sta_compare\overall_sta_lf8'];
    print('-dpng',t_names);
    close
    
     for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_sta_lf3(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-0.2 0.2]);grid
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\sta_compare\overall_sta_lf3'];
    print('-dpng',t_names);
    close
   
        for i = [1 7]
   plot(lags/Fsd,avg_sta_lf2(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-0.2 0.2]);grid
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\sta_compare\overall_sta_lf2'];
    print('-dpng',t_names);
    close
