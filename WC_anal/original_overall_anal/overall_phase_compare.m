clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\phase\phase_data

%%
diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)

    get_cells = find(cell_type == diff_cell_types(c));
    
    avg_up_phase(c,:) = mean(up_dist(get_cells,:));
    uci_up_phase = avg_up_phase(c,:)+tinv(0.975,length(get_cells))*std(up_dist(get_cells,:))/sqrt(length(get_cells));
    lci_up_phase = avg_up_phase(c,:)-tinv(0.975,length(get_cells))*std(up_dist(get_cells,:))/sqrt(length(get_cells));

    avg_down_phase(c,:) = mean(down_dist(get_cells,:));
    uci_down_phase = avg_down_phase(c,:)+tinv(0.975,length(get_cells))*std(down_dist(get_cells,:))/sqrt(length(get_cells));
    lci_down_phase = avg_down_phase(c,:)-tinv(0.975,length(get_cells))*std(down_dist(get_cells,:))/sqrt(length(get_cells));

    h = polar(phase_range/360*2*pi,avg_up_phase(c,:))
    hold on
    set(h,'linewidth',2)
    h = polar(phase_range/360*2*pi,avg_down_phase(c,:),'r')
    set(h,'linewidth',2)
    polar(phase_range/360*2*pi,uci_up_phase,'--')
    polar(phase_range/360*2*pi,lci_up_phase,'--')
    polar(phase_range/360*2*pi,uci_down_phase,'r--')
    polar(phase_range/360*2*pi,lci_down_phase,'r--')
    legend('MP Up transition','MP Down transition')
    t_names = ['C:\WC_Germany\overall_calcs\phase_compare\' num2str(c)];
    print('-dpng',t_names);
    close
    
end


%% Plot direct comparisons
%MP spec compare
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   h= polar(phase_range/360*2*pi,avg_up_phase(i,:)) 
   hold on
   set(h,'color',cmap(i,:))
end
set(gca,'xlim',[-0.1 0.1])
    t_names = ['C:\WC_Germany\overall_calcs\phase_compare\up_transition_compare'];
    print('-dpng',t_names);

%MP spec compare
figure
for i = 1:length(diff_cell_types)
   h = polar(phase_range/360*2*pi,avg_down_phase(i,:)) 
   hold on
      set(h,'color',cmap(i,:))

end
set(gca,'xlim',[-0.1 0.1])

    t_names = ['C:\WC_Germany\overall_calcs\phase_compare\down_transition_compare'];
    print('-dpng',t_names);
