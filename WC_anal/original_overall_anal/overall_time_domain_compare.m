clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\time_domain\time_domain_data

%% smooth autocors
% sm_win = 5;
% for i = 1:length(over_dir)
%     tot_wcv_acorr(i,:) = jmm_smooth_1d(tot_wcv_acorr(i,:),sm_win);
%     tot_lf8_acorr(i,:) = jmm_smooth_1d(tot_lf8_acorr(i,:),sm_win);
%     tot_lf3_acorr(i,:) = jmm_smooth_1d(tot_lf3_acorr(i,:),sm_win);
%     tot_lf2_acorr(i,:) = jmm_smooth_1d(tot_lf2_acorr(i,:),sm_win);
% %     tot_lf3s_acorr(i,:) = jmm_smooth_1d(tot_lf3s_acorr(i,:),sm_win);
% end

%% Acorrs
diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)

    get_cells = find(cell_type == diff_cell_types(c));
    
    avg_w_acorr(c,:) = mean(tot_wcv_acorr(get_cells,:));
    uci_wcor = avg_w_acorr(c,:)+tinv(0.975,length(get_cells))*std(tot_wcv_acorr(get_cells,:))/sqrt(length(get_cells));
    lci_wcor = avg_w_acorr(c,:)-tinv(0.975,length(get_cells))*std(tot_wcv_acorr(get_cells,:))/sqrt(length(get_cells));
    
    avg_8_acorr(c,:) = mean(tot_lf8_acorr(get_cells,:));
    uci_8cor = avg_8_acorr(c,:)+tinv(0.975,length(get_cells))*std(tot_lf8_acorr(get_cells,:))/sqrt(length(get_cells));
    lci_8cor = avg_8_acorr(c,:)-tinv(0.975,length(get_cells))*std(tot_lf8_acorr(get_cells,:))/sqrt(length(get_cells));

    avg_3_acorr(c,:) = mean(tot_lf3_acorr(get_cells,:));
    uci_3cor = avg_3_acorr(c,:)+tinv(0.975,length(get_cells))*std(tot_lf3_acorr(get_cells,:))/sqrt(length(get_cells));
    lci_3cor = avg_3_acorr(c,:)-tinv(0.975,length(get_cells))*std(tot_lf3_acorr(get_cells,:))/sqrt(length(get_cells));

    avg_2_acorr(c,:) = mean(tot_lf2_acorr(get_cells,:));
    uci_2cor = avg_2_acorr(c,:)+tinv(0.975,length(get_cells))*std(tot_lf2_acorr(get_cells,:))/sqrt(length(get_cells));
    lci_2cor = avg_2_acorr(c,:)-tinv(0.975,length(get_cells))*std(tot_lf2_acorr(get_cells,:))/sqrt(length(get_cells));

    %         avg_3s_acorr(c,:) = mean(tot_lf3s_acorr(get_cells,:));
%     uci_3scor = avg_3s_acorr(c,:)+tinv(0.975,length(get_cells))*std(tot_lf3s_acorr(get_cells,:))/sqrt(length(get_cells));
%     lci_3scor = avg_3s_acorr(c,:)-tinv(0.975,length(get_cells))*std(tot_lf3s_acorr(get_cells,:))/sqrt(length(get_cells));

    plot(lags/Fsd,avg_w_acorr(c,:),'linewidth',2)
    hold on
    plot(lags/Fsd,avg_8_acorr(c,:),'r','linewidth',2)
    plot(lags/Fsd,avg_3_acorr(c,:),'g','linewidth',2)
    plot(lags/Fsd,avg_2_acorr(c,:),'k','linewidth',2)
%     plot(lags/Fsd,avg_3s_acorr(c,:),'k','linewidth',2)
    legend('MP-acorr','LF8-acorr','LF3-acorr','LF2-acorr')
    plot(lags/Fsd,uci_wcor,'--')
    plot(lags/Fsd,lci_wcor,'--')
        plot(lags/Fsd,uci_8cor,'r--')
    plot(lags/Fsd,lci_8cor,'r--')
    plot(lags/Fsd,uci_3cor,'g--')
    plot(lags/Fsd,lci_3cor,'g--')
    plot(lags/Fsd,uci_2cor,'k--')
    plot(lags/Fsd,lci_2cor,'k--')

    xlim([0 10])
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\overall_acorr_' num2str(c)];
    print('-dpng',t_names);
    close

    
    
end

%% xcors
for c = 1:length(diff_cell_types)

    get_cells = find(cell_type == diff_cell_types(c));
    
    avg_w8_x(c,:) = mean(tot_w8_x(get_cells,:));
    uci_w8_x = avg_w8_x(c,:)+tinv(0.975,length(get_cells))*std(tot_w8_x(get_cells,:))/sqrt(length(get_cells));
    lci_w8_x = avg_w8_x(c,:)-tinv(0.975,length(get_cells))*std(tot_w8_x(get_cells,:))/sqrt(length(get_cells));
    
    avg_w3_x(c,:) = mean(tot_w3_x(get_cells,:));
    uci_w3_x = avg_w3_x(c,:)+tinv(0.975,length(get_cells))*std(tot_w3_x(get_cells,:))/sqrt(length(get_cells));
    lci_w3_x = avg_w3_x(c,:)-tinv(0.975,length(get_cells))*std(tot_w3_x(get_cells,:))/sqrt(length(get_cells));

%     avg_w3s_x(c,:) = mean(tot_w3s_x(get_cells,:));
%     uci_w3s_x = avg_w3s_x(c,:)+tinv(0.975,length(get_cells))*std(tot_w3s_x(get_cells,:))/sqrt(length(get_cells));
%     lci_w3s_x = avg_w3s_x(c,:)-tinv(0.975,length(get_cells))*std(tot_w3s_x(get_cells,:))/sqrt(length(get_cells));
    avg_w2_x(c,:) = mean(tot_w2_x(get_cells,:));
    uci_w2_x = avg_w2_x(c,:)+tinv(0.975,length(get_cells))*std(tot_w2_x(get_cells,:))/sqrt(length(get_cells));
    lci_w2_x = avg_w2_x(c,:)-tinv(0.975,length(get_cells))*std(tot_w2_x(get_cells,:))/sqrt(length(get_cells));

    plot(lags/Fsd,avg_w8_x(c,:),'linewidth',2)
    hold on
    plot(lags/Fsd,avg_w3_x(c,:),'r','linewidth',2)
    plot(lags/Fsd,avg_w2_x(c,:),'k','linewidth',2)
%     plot(lags/Fsd,avg_w3s_x(c,:),'k','linewidth',2)
    legend('MP-Lf8','MP-LF3','MP-LF2')
    plot(lags/Fsd,uci_w8_x,'--')
    plot(lags/Fsd,lci_w8_x,'--')
        plot(lags/Fsd,uci_w3_x,'r--')
    plot(lags/Fsd,lci_w3_x,'r--')
    plot(lags/Fsd,uci_w2_x,'k--')
    plot(lags/Fsd,lci_w2_x,'k--')

    xlim([-10 10])
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\overall_xcorr_wb' num2str(c)];
    print('-dpng',t_names);
    xlim([-2 2])
    grid on
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\overall_xcorr_' num2str(c)];
    print('-dpng',t_names);
    close

    
    
end

%% Plot direct comparisons Acorrs
%MP-acorr
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_w_acorr(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    xlim([0 10])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\mp_acorr_compare'];
    print('-dpng',t_names);
    close

%LFP-acorr
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_8_acorr(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    xlim([0 10])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\lf8_acorr_compare'];
    print('-dpng',t_names);
    close

    %LF3-acorr
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_3_acorr(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    xlim([0 10])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\lf3_acorr_compare'];
    print('-dpng',t_names);
    close

    %LF3s-acorr
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_2_acorr(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    xlim([0 10])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\lf2_acorr_compare'];
    print('-dpng',t_names);
    close

    %MP-LF8_xcor
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_w8_x(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    xlim([-10 10])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\mp_lf8_xcorr_compare_wb'];
    print('-dpng',t_names);
    xlim([-2 2])
    grid on
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\mp_lf8_xcorr_compare_'];
    print('-dpng',t_names);
    close
    
        %MP-LF3_xcor
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_w3_x(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    xlim([-10 10])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\mp_lf3_xcorr_compare_wb'];
    print('-dpng',t_names);
    xlim([-2 2])
    grid on
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\mp_lf3_xcorr_compare_'];
    print('-dpng',t_names);
    close

            %MP-LF3s_xcor
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags/Fsd,avg_w2_x(i,:),'color',cmap(i,:),'linewidth',2) 
   hold on
end
    xlim([-10 10])
legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\mp_lf2_xcorr_compare_wb'];
    print('-dpng',t_names);
    xlim([-2 2])
    grid on
    t_names = ['C:\WC_Germany\overall_calcs\time_domain_compare\mp_lf2_xcorr_compare_'];
    print('-dpng',t_names);
    close