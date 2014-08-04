clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\burst\burst_data

diff_cell_types = unique(cell_type);

for d = 1:length(over_dir)
    isi_hist(d,:) = jmm_smooth_1d(isi_hist(d,:),2);
end

for c = 1:length(diff_cell_types)

    get_cells = find(cell_type == diff_cell_types(c));
    
    avg_isi_dist(c,:) = mean(isi_hist(get_cells,:));
    uci_up = avg_isi_dist(c,:)+tinv(0.975,length(get_cells))*std(isi_hist(get_cells,:))/sqrt(length(get_cells));
    lci_up = avg_isi_dist(c,:)-tinv(0.975,length(get_cells))*std(isi_hist(get_cells,:))/sqrt(length(get_cells));

    figure
    plot(isi_range,avg_isi_dist(c,:),'linewidth',2)
    hold on
    plot(isi_range,uci_up,'--')
    plot(isi_range,lci_up,'--')
    xlim([0 100])
    t_names = ['C:\WC_Germany\overall_calcs\burst_compare\' num2str(c)];
    print('-dpng',t_names);
    close
end


%% Plot direct comparisons
%MP up compare
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(isi_range,avg_isi_dist(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end

    xlim([0 100])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\burst_compare\overall'];
    print('-dpng',t_names);
    close

