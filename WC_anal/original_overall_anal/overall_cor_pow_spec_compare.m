clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\overall_info_file_data
load C:\WC_Germany\overall_calcs\cor_pow_spec\cor_pow_spec_data

diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)-2

    get_cells = find(cell_type == diff_cell_types(c) & bad_lf2 == 0 & bad_lf3 == 0);
    avg_spec(c,:) = mean(Sw(get_cells,:));
    uci_spec = avg_spec(c,:)+tinv(0.975,length(get_cells))*std(Sw(get_cells,:))/sqrt(length(get_cells));
    lci_spec = avg_spec(c,:)-tinv(0.975,length(get_cells))*std(Sw(get_cells,:))/sqrt(length(get_cells));
    
        plot(f,10*log10(avg_spec(c,:)),'linewidth',2)
    hold on
    plot(f,10*log10(uci_spec),'--')
    plot(f,10*log10(lci_spec),'--')

    xlim([0 20])
    t_names = ['C:\WC_Germany\overall_calcs\cor_pow_spec_compare\wb_' num2str(c)];
    print('-dpng',t_names);
    xlim([0 1])
    t_names = ['C:\WC_Germany\overall_calcs\cor_pow_spec_compare\' num2str(c)];
    print('-dpng',t_names);
    close

    
end

plot(f,10*log10(avg_spec),'linewidth',2)
legend('MEC L3','MEC L2','MEC L5','LEC 3','LEC 2','LEC ?')
xlim([0 20])
    t_names = ['C:\WC_Germany\overall_calcs\cor_pow_spec_compare\overall_wb'];
    print('-dpng',t_names);
    xlim([0 1])
    t_names = ['C:\WC_Germany\overall_calcs\cor_pow_spec_compare\overall'];
    print('-dpng',t_names);
    close
