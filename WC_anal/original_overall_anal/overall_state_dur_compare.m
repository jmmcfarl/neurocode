clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

up_range = [0.3 20];
down_range = [0.3 20];
numBins = 100;

for d = 1:length(over_dir)
    
    [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
    [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);
    up_hist(d,:) = jmm_smooth_1d(up_hist(d,:),2);
    up_hist8(d,:) = jmm_smooth_1d(up_hist8(d,:),2);

    [down_hist(d,:),down_grid] = log_hist(synch_down_dur{d},down_range,numBins);
    [down_hist8(d,:),down_grid] = log_hist(synch_down_dur8{d},down_range,numBins);
    down_hist(d,:) = jmm_smooth_1d(down_hist(d,:),2);
    down_hist8(d,:) = jmm_smooth_1d(down_hist8(d,:),2);
    
end

%% cycle through all cell types and compute the average up and down
%% distribution of MP and LFPs
diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)

    get_cells = find(cell_type == diff_cell_types(c));
    
    avg_up(c,:) = mean(up_hist(get_cells,:));
    uci_up = avg_up(c,:)+tinv(0.975,length(get_cells))*std(up_hist(get_cells,:))/sqrt(length(get_cells));
    lci_up = avg_up(c,:)-tinv(0.975,length(get_cells))*std(up_hist(get_cells,:))/sqrt(length(get_cells));

    avg_up8(c,:) = mean(up_hist8(get_cells,:));
    uci_up8 = avg_up8(c,:)+tinv(0.975,length(get_cells))*std(up_hist8(get_cells,:))/sqrt(length(get_cells));
    lci_up8 = avg_up8(c,:)-tinv(0.975,length(get_cells))*std(up_hist8(get_cells,:))/sqrt(length(get_cells));
    
    avg_down(c,:) = mean(down_hist(get_cells,:));
    uci_down = avg_down(c,:)+tinv(0.975,length(get_cells))*std(down_hist(get_cells,:))/sqrt(length(get_cells));
    lci_down = avg_down(c,:)-tinv(0.975,length(get_cells))*std(down_hist(get_cells,:))/sqrt(length(get_cells));

    avg_down8(c,:) = mean(down_hist8(get_cells,:));
    uci_down8 = avg_down8(c,:)+tinv(0.975,length(get_cells))*std(down_hist8(get_cells,:))/sqrt(length(get_cells));
    lci_down8 = avg_down8(c,:)-tinv(0.975,length(get_cells))*std(down_hist8(get_cells,:))/sqrt(length(get_cells));

    figure
    plot(up_grid,avg_up(c,:),'linewidth',2)
    hold on
    plot(up_grid,uci_up,'--')
    plot(up_grid,lci_up,'--')
    plot(up_grid,avg_up8(c,:),'r','linewidth',2)
    plot(up_grid,uci_up8,'r--')
    plot(up_grid,lci_up8,'r--')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-3 0.07])
    xlim([0 20])
    t_names = ['C:\WC_Germany\overall_calcs\state_dur_compare\up_state_' num2str(c)];
    print('-dpng',t_names);
    close

    figure
    plot(down_grid,avg_down(c,:),'linewidth',2)
    hold on
    plot(down_grid,uci_down,'--')
    plot(down_grid,lci_down,'--')
    plot(down_grid,avg_down8(c,:),'r','linewidth',2)
    plot(down_grid,uci_down8,'r--')
    plot(down_grid,lci_down8,'r--')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
        ylim([1e-3 0.07])
    xlim([0 20])
    t_names = ['C:\WC_Germany\overall_calcs\state_dur_compare\down_state_' num2str(c)];
    print('-dpng',t_names);
    close
    
end


%% Plot direct comparisons
%MP up compare
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(up_grid,avg_up(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
        ylim([1e-3 0.07])
    xlim([0 20])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\state_dur_compare\up_compare_mp_t'];
    print('-dpng',t_names);
    close

%MP down compare
figure
for i = 1:length(diff_cell_types)
   plot(down_grid,avg_down(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
        ylim([1e-3 0.07])
    xlim([0 20])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\state_dur_compare\down_compare_mp_t'];
    print('-dpng',t_names);
    close
    
    %LFP up compare
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(up_grid,avg_up8(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
        ylim([1e-3 0.07])
    xlim([0 20])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\state_dur_compare\up_compare_lfp_t'];
    print('-dpng',t_names);
    close

%LFP down compare
figure
for i = 1:length(diff_cell_types)
   plot(down_grid,avg_down8(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
        ylim([1e-3 0.07])
    xlim([0 20])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\state_dur_compare\down_compare_lfp_t'];
    print('-dpng',t_names);
    close
