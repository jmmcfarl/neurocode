clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\power_spectra\power_spec_data
load C:\WC_Germany\overall_calcs\power_spectra\power_spec_data_hipp

%% cycle through all cell types and compute the average up and down
%% distribution of MP and LFPs
diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)

    get_cells = find(cell_type == diff_cell_types(c));
    
    avg_pow(c,:) = mean(Sw(get_cells,:));
    uci_pow = avg_pow(c,:)+tinv(0.975,length(get_cells))*std(Sw(get_cells,:))/sqrt(length(get_cells));
    lci_pow = avg_pow(c,:)-tinv(0.975,length(get_cells))*std(Sw(get_cells,:))/sqrt(length(get_cells));

    avg_pow8(c,:) = mean(S8(get_cells,:));
    uci_pow8 = avg_pow8(c,:)+tinv(0.975,length(get_cells))*std(S8(get_cells,:))/sqrt(length(get_cells));
    lci_pow8 = avg_pow8(c,:)-tinv(0.975,length(get_cells))*std(S8(get_cells,:))/sqrt(length(get_cells));
    
    avg_pow3(c,:) = mean(S3(get_cells,:));
    uci_pow3 = avg_pow3(c,:)+tinv(0.975,length(get_cells))*std(S3(get_cells,:))/sqrt(length(get_cells));
    lci_pow3 = avg_pow3(c,:)-tinv(0.975,length(get_cells))*std(S3(get_cells,:))/sqrt(length(get_cells));
    
    avg_pow3s(c,:) = mean(S3s(get_cells,:));
    uci_pow3s = avg_pow3s(c,:)+tinv(0.975,length(get_cells))*std(S3s(get_cells,:))/sqrt(length(get_cells));
    lci_pow3s = avg_pow3s(c,:)-tinv(0.975,length(get_cells))*std(S3s(get_cells,:))/sqrt(length(get_cells));

    figure
    plot(f,10*log10(avg_pow(c,:)),'linewidth',2)
    hold on
    plot(f,10*log10(avg_pow8(c,:)),'r','linewidth',2)
    plot(f,10*log10(avg_pow3(c,:)),'g','linewidth',2)
    plot(f,10*log10(avg_pow3s(c,:)),'k','linewidth',2)
    plot(f,10*log10(uci_pow),'--')
    plot(f,10*log10(lci_pow),'--')
    plot(f,10*log10(uci_pow8),'r--')
    plot(f,10*log10(lci_pow8),'r--')
    plot(f,10*log10(uci_pow3),'g--')
    plot(f,10*log10(lci_pow3),'g--')
    plot(f,10*log10(uci_pow3s),'k--')
    plot(f,10*log10(lci_pow3s),'k--')
    legend('MP','Cort LFP','Hipp LFP','Orth Hipp LFP')
    set(gca,'xscale','log')
    xlim([0 45])
    t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\wb_' num2str(c)];
    print('-dpng',t_names);
    
    set(gca,'xscale','linear')
    xlim([0 2])
    t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\' num2str(c)];
    print('-dpng',t_names);
    close
    
end


%% Plot direct comparisons
%MP spec compare
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(f,10*log10(avg_pow(i,:)),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    set(gca,'xscale','log')
    xlim([0 45])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\mp_compare_wb'];
    print('-dpng',t_names);
set(gca,'xscale','linear')
xlim([0 1])
    t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\mp_compare_'];
    print('-dpng',t_names);
    close

%LFP pow compare
figure
for i = 1:length(diff_cell_types)
   plot(f,10*log10(avg_pow8(i,:)),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    set(gca,'xscale','log')
    xlim([0 45])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\lfp_compare_wb'];
    print('-dpng',t_names);
    set(gca,'xscale','linear')
    xlim([0 1])
        t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\lfp_compare_'];
    print('-dpng',t_names);
    close

%Hipp LFP pow compare
figure
for i = 1:length(diff_cell_types)
   plot(f,10*log10(avg_pow3(i,:)),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    set(gca,'xscale','log')
    xlim([0 45])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\hipp_lfp_compare_wb'];
    print('-dpng',t_names);
    set(gca,'xscale','linear')
    xlim([0 1])
        t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\hipp_lfp_compare_'];
    print('-dpng',t_names);
    close

    %Orth Hipp LFP pow compare
figure
for i = 1:length(diff_cell_types)
   plot(f,10*log10(avg_pow3s(i,:)),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    set(gca,'xscale','log')
    xlim([0 45])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\orth_hipp_lfp_compare_wb'];
    print('-dpng',t_names);
    set(gca,'xscale','linear')
    xlim([0 1])
        t_names = ['C:\WC_Germany\overall_calcs\pow_spec_compare\orth_hipp_lfp_compare_'];
    print('-dpng',t_names);
    close