clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\lag_phase_data

%% Calculate up state distribution for each cell individually
nbins = 100;
for d = 1:length(over_dir)
    mp_updur_lfp_cor{d}(isnan(mp_updur_lfp_cor{d})) = [];
    temp = abs(zscore(mp_updur_lfp_cor{d}));
   [y,x] = hist(mp_updur_lfp_cor{d}(temp < 5),nbins);
   bar(x,y)
    tname = ['C:\WC_Germany\overall_calcs\persistence_calcs\' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
    
end

%% calculate smoothed individual cell distros for averaging
hist_range = linspace(0,7,100);
clear y x
for d = 1:length(over_dir)
    
    mp_updur_lfp_cor{d}(isnan(mp_updur_lfp_cor{d})) = [];
    temp = abs(zscore(mp_updur_lfp_cor{d}));
    y(d,:) = hist(mp_updur_lfp_cor{d}(temp < 5),hist_range);
    y(d,:) = jmm_smooth_1d(y(d,:),2);
    
end

diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)

    get_cells = find(cell_type == diff_cell_types(c));
    
    avg_dist(c,:) = mean(y(get_cells,:));
    uci_dist = avg_dist(c,:)+tinv(0.975,length(get_cells))*std(y(get_cells,:))/sqrt(length(get_cells));
    lci_dist = avg_dist(c,:)-tinv(0.975,length(get_cells))*std(y(get_cells,:))/sqrt(length(get_cells));

    figure
    plot(hist_range,avg_dist(c,:),'linewidth',2)
    hold on
    plot(hist_range,uci_dist,'--')
    plot(hist_range,lci_dist,'--')
    set(gca,'yscale','log')
    xlim([0 7])
    t_names = ['C:\WC_Germany\overall_calcs\persistence_compare\avg_' num2str(c)];
    print('-dpng',t_names);
    close
    
end


%% Plot direct comparisons
%MP spec compare
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(hist_range,avg_dist(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    set(gca,'yscale','log')
    xlim([0 7])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\persistence_compare\overall_compare'];
    print('-dpng',t_names);
    close

