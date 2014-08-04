clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\coherence\coherence_data
load C:\WC_Germany\overall_calcs\coherence\coherence_data_hipp
load C:\WC_Germany\overall_calcs\overall_info_file_data

%% smooth all coherence functions
% for i = 1:length(over_dir)
%     
%    Cmn(i,:) = jmm_smooth_1d(Cmn(i,:),5);
%    Cmn3(i,:) = jmm_smooth_1d(Cmn3(i,:),5);
%    Cmn3s(i,:) = jmm_smooth_1d(Cmn3s(i,:),5);
%    Cmn2(i,:) = jmm_smooth_1d(Cmn2(i,:),5);
%    
% end


%% cycle through all cell types and compute the average up and down
%% distribution of MP and LFPs
diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)

%     get_cells = find(cell_type == diff_cell_types(c) & bad_lf2 == 0 & bad_lf3 == 0);
  get_cells = find(cell_type == diff_cell_types(c) & bad_lf2 == 0 & bad_lf3 == 0 & hp_lf2 == 0);

    avg_coh(c,:) = mean(Cmn(get_cells,:));
    uci_coh = avg_coh(c,:)+tinv(0.975,length(get_cells))*std(Cmn(get_cells,:))/sqrt(length(get_cells));
    lci_coh = avg_coh(c,:)-tinv(0.975,length(get_cells))*std(Cmn(get_cells,:))/sqrt(length(get_cells));
    
    avg_coh3(c,:) = mean(Cmn3(get_cells,:));
    uci_coh3 = avg_coh3(c,:)+tinv(0.975,length(get_cells))*std(Cmn3(get_cells,:))/sqrt(length(get_cells));
    lci_coh3 = avg_coh3(c,:)-tinv(0.975,length(get_cells))*std(Cmn3(get_cells,:))/sqrt(length(get_cells));
    
%     avg_coh3s(c,:) = mean(Cmn3s(get_cells,:));
%     uci_coh3s = avg_coh3s(c,:)+tinv(0.975,length(get_cells))*std(Cmn3s(get_cells,:))/sqrt(length(get_cells));
%     lci_coh3s = avg_coh3s(c,:)-tinv(0.975,length(get_cells))*std(Cmn3s(get_cells,:))/sqrt(length(get_cells));

    avg_coh2(c,:) = mean(Cmn2(get_cells,:));
    uci_coh2 = avg_coh2(c,:)+tinv(0.975,length(get_cells))*std(Cmn2(get_cells,:))/sqrt(length(get_cells));
    lci_coh2 = avg_coh2(c,:)-tinv(0.975,length(get_cells))*std(Cmn2(get_cells,:))/sqrt(length(get_cells));

    
    plot(f,avg_coh(c,:),'linewidth',2)
    hold on
    plot(f,avg_coh2(c,:),'g','linewidth',2)
    plot(f,avg_coh3(c,:),'r','linewidth',2)
%     plot(f,avg_coh3s(c,:),'k','linewidth',2)
    legend('MP-CortLFP','MP-deepHipp','MP-HippLFP')
    plot(f,uci_coh,'--')
    plot(f,lci_coh,'--')
    plot(f,uci_coh2,'g--')
    plot(f,lci_coh2,'g--')
        plot(f,uci_coh3,'r--')
    plot(f,lci_coh3,'r--')
%     plot(f,uci_coh3s,'k--')
%     plot(f,lci_coh3s,'k--')

%     xlim([0 20])
%     t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\overall_wb_' num2str(c)];
%     print('-dpng',t_names);
    xlim([0 1])
    ylim([0 1])
    t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\overall_nofilt_' num2str(c)];
    print('-dpng',t_names);
    close

    
end


%% Plot direct comparisons
% MP-LFP
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(f,avg_coh(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    xlim([0 20])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\coh_compare_wb'];
%     print('-dpng',t_names);
xlim([0 1]);    
ylim([0 1])

    t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\coh_compare_nofilt_'];
    print('-dpng',t_names);
    close

%MP-hipp LFP
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(f,avg_coh3(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
%     xlim([0 20])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\hipp_coh_compare_wb'];
%     print('-dpng',t_names);
xlim([0 1]);   
ylim([0 1])

    t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\hipp_coh_compare_nofilt'];
    print('-dpng',t_names);
    close

%     %MP-orth-hipp LFP
% cmap = colormap(jet(length(diff_cell_types)));
% figure
% for i = 1:length(diff_cell_types)
%    plot(f,avg_coh3s(i,:),'color',cmap(i,:),'linewidth',2) 
%    hold on
% end
%     xlim([0 45])
% legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\orth_hipp_coh_compare_wb'];
%     print('-dpng',t_names);
% % xlim([0 1])
% %     t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\orth_hipp_coh_compare_'];
% %     print('-dpng',t_names);
%     close
    
    cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(f,avg_coh2(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    xlim([0 20])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
% %     t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\deep_hipp_coh_compare_wb'];
%     print('-dpng',t_names);
xlim([0 1]);    
ylim([0 1])

    t_names = ['C:\WC_Germany\overall_calcs\coherence_compare\deep_hipp_coh_compare_nofilt_'];
    print('-dpng',t_names);
    close
    