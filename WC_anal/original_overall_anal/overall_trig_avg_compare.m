%%
 
clear all
close all
 

load('C:\WC_Germany\overall_calcs\overall_dir.mat')

%% cycle through orth lf3 and convert
load C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_data

% z_time = find(lags > 0,1,'first');
% flip_sign = ones(size(over_dir));

for i = 1:(length(over_dir))

% 	if mp_utrig_lf3(i,z_time) < 0
% 		flip_sign(i) = -1;
% 	end
% 	mp_utrig_lf3s(i,:) = flip_sign(i)*mp_utrig_lf3(i,:);
% 	mp_dtrig_lf3s(i,:) = flip_sign(i)*mp_dtrig_lf3(i,:);
% 	lf8_utrig_lf3s(i,:) = flip_sign(i)*lf8_utrig_lf3(i,:);
% 	lf8_dtrig_lf3s(i,:) = flip_sign(i)*lf8_dtrig_lf3(i,:);	

end


load C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_data
load C:\WC_Germany\overall_calcs\overall_info_file_data

% beg_lag = find(lags>-2,1,'first');
% end_lag = find(lags>2,1,'first');
% 


%% cycle through all cell types and compute the average trig averages
diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));
 
%% mp utrig
    avg_mp_utrig_mp(c,:) = mean(abs(mp_utrig_mp(get_cells,:)));
    uci_mp_utrig_mp = avg_mp_utrig_mp(c,:)+tinv(0.975,length(get_cells))*std(abs(mp_utrig_mp(get_cells,:)))/sqrt(length(get_cells));
    lci_mp_utrig_mp = avg_mp_utrig_mp(c,:)-tinv(0.975,length(get_cells))*std(abs(mp_utrig_mp(get_cells,:)))/sqrt(length(get_cells));
    
    avg_mp_utrig_lf8(c,:) = mean(abs(mp_utrig_lf8(get_cells,:)));
    uci_mp_utrig_lf8 = avg_mp_utrig_lf8(c,:)+tinv(0.975,length(get_cells))*std(abs(mp_utrig_lf8(get_cells,:)))/sqrt(length(get_cells));
    lci_mp_utrig_lf8 = avg_mp_utrig_lf8(c,:)-tinv(0.975,length(get_cells))*std(abs(mp_utrig_lf8(get_cells,:)))/sqrt(length(get_cells));
  
        get_cells = find(cell_type == diff_cell_types(c) & bad_lf3 == 0);

    avg_mp_utrig_lf3(c,:) = mean(abs(mp_utrig_lf3(get_cells,:)));
    uci_mp_utrig_lf3 = avg_mp_utrig_lf3(c,:)+tinv(0.975,length(get_cells))*std(abs(mp_utrig_lf3(get_cells,:)))/sqrt(length(get_cells));
    lci_mp_utrig_lf3 = avg_mp_utrig_lf3(c,:)-tinv(0.975,length(get_cells))*std(abs(mp_utrig_lf3(get_cells,:)))/sqrt(length(get_cells));
    
%     avg_mp_utrig_lf3s(c,:) = mean(mp_utrig_lf3s(get_cells,:));
%     uci_mp_utrig_lf3s = avg_mp_utrig_lf3s(c,:)+tinv(0.975,length(get_cells))*std(mp_utrig_lf3s(get_cells,:))/sqrt(length(get_cells));
%     lci_mp_utrig_lf3s = avg_mp_utrig_lf3s(c,:)-tinv(0.975,length(get_cells))*std(mp_utrig_lf3s(get_cells,:))/sqrt(length(get_cells));
    get_cells = find(cell_type == diff_cell_types(c) & bad_lf2 == 0);

    avg_mp_utrig_lf2(c,:) = mean(abs(mp_utrig_lf2(get_cells,:)));
    uci_mp_utrig_lf2 = avg_mp_utrig_lf2(c,:)+tinv(0.975,length(get_cells))*std(abs(mp_utrig_lf2(get_cells,:)))/sqrt(length(get_cells));
    lci_mp_utrig_lf2 = avg_mp_utrig_lf2(c,:)-tinv(0.975,length(get_cells))*std(abs(mp_utrig_lf2(get_cells,:)))/sqrt(length(get_cells));
    
    plot(lags,avg_mp_utrig_mp(c,:),'linewidth',1)
    hold on
    plot(lags,avg_mp_utrig_lf8(c,:),'r','linewidth',1)
    plot(lags,avg_mp_utrig_lf3(c,:),'g','linewidth',1)
% plot(lags,avg_mp_utrig_lf3s(c,:),'m','linewidth',1)
    plot(lags,avg_mp_utrig_lf2(c,:),'k','linewidth',1)
    legend('MP','LF8','LF3','LF2')
    plot(lags,uci_mp_utrig_mp,'--')
    plot(lags,lci_mp_utrig_mp,'--')
        plot(lags,uci_mp_utrig_lf8,'r--')
    plot(lags,lci_mp_utrig_lf8,'r--')
    plot(lags,uci_mp_utrig_lf3,'g--')
    plot(lags,lci_mp_utrig_lf3,'g--')    
%     plot(lags,uci_mp_utrig_lf3s,'m--')
%     plot(lags,lci_mp_utrig_lf3s,'m--')    

    plot(lags,uci_mp_utrig_lf2,'k--')
    plot(lags,lci_mp_utrig_lf2,'k--')
    xlim([-5 10]);
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_mp_utrig_' num2str(c)];
    print('-dpng',t_names);
    xlim([-2 2]); grid
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_mp_utrig_zoom_' num2str(c)];
    print('-dpng',t_names);
    close
 
%% mp dtrig
    get_cells = find(cell_type == diff_cell_types(c));

    avg_mp_dtrig_mp(c,:) = mean(abs(mp_dtrig_mp(get_cells,:)));
    uci_mp_dtrig_mp = avg_mp_dtrig_mp(c,:)+tinv(0.975,length(get_cells))*std(abs(mp_dtrig_mp(get_cells,:)))/sqrt(length(get_cells));
    lci_mp_dtrig_mp = avg_mp_dtrig_mp(c,:)-tinv(0.975,length(get_cells))*std(abs(mp_dtrig_mp(get_cells,:)))/sqrt(length(get_cells));
    
    avg_mp_dtrig_lf8(c,:) = mean(abs(mp_dtrig_lf8(get_cells,:)));
    uci_mp_dtrig_lf8 = avg_mp_dtrig_lf8(c,:)+tinv(0.975,length(get_cells))*std(abs(mp_dtrig_lf8(get_cells,:)))/sqrt(length(get_cells));
    lci_mp_dtrig_lf8 = avg_mp_dtrig_lf8(c,:)-tinv(0.975,length(get_cells))*std(abs(mp_dtrig_lf8(get_cells,:)))/sqrt(length(get_cells));
           
    get_cells = find(cell_type == diff_cell_types(c) & bad_lf3 == 0);

    avg_mp_dtrig_lf3(c,:) = mean(abs(mp_dtrig_lf3(get_cells,:)));
    uci_mp_dtrig_lf3 = avg_mp_dtrig_lf3(c,:)+tinv(0.975,length(get_cells))*std(abs(mp_dtrig_lf3(get_cells,:)))/sqrt(length(get_cells));
    lci_mp_dtrig_lf3 = avg_mp_dtrig_lf3(c,:)-tinv(0.975,length(get_cells))*std(abs(mp_dtrig_lf3(get_cells,:)))/sqrt(length(get_cells));
        
% 	avg_mp_dtrig_lf3s(c,:) = mean(mp_dtrig_lf3s(get_cells,:));
%     uci_mp_dtrig_lf3s = avg_mp_dtrig_lf3s(c,:)+tinv(0.975,length(get_cells))*std(mp_dtrig_lf3s(get_cells,:))/sqrt(length(get_cells));
%     lci_mp_dtrig_lf3s = avg_mp_dtrig_lf3s(c,:)-tinv(0.975,length(get_cells))*std(mp_dtrig_lf3s(get_cells,:))/sqrt(length(get_cells));
    
get_cells = find(cell_type == diff_cell_types(c) & bad_lf2 == 0);

    avg_mp_dtrig_lf2(c,:) = mean(abs(mp_dtrig_lf2(get_cells,:)));
    uci_mp_dtrig_lf2 = avg_mp_dtrig_lf2(c,:)+tinv(0.975,length(get_cells))*std(abs(mp_dtrig_lf2(get_cells,:)))/sqrt(length(get_cells));
    lci_mp_dtrig_lf2 = avg_mp_dtrig_lf2(c,:)-tinv(0.975,length(get_cells))*std(abs(mp_dtrig_lf2(get_cells,:)))/sqrt(length(get_cells));
 
     
 
    plot(lags,avg_mp_dtrig_mp(c,:),'linewidth',1)
    hold on
    plot(lags,avg_mp_dtrig_lf8(c,:),'r','linewidth',1)
    plot(lags,avg_mp_dtrig_lf3(c,:),'g','linewidth',1)
% 	plot(lags,avg_mp_dtrig_lf3s(c,:),'m','linewidth',1)
    plot(lags,avg_mp_dtrig_lf2(c,:),'k','linewidth',1)
    legend('MP','LF8','LF3','LF2')
    plot(lags,uci_mp_dtrig_mp,'--')
    plot(lags,lci_mp_dtrig_mp,'--')
        plot(lags,uci_mp_dtrig_lf8,'r--')
    plot(lags,lci_mp_dtrig_lf8,'r--')
    plot(lags,uci_mp_dtrig_lf3,'g--')
    plot(lags,lci_mp_dtrig_lf3,'g--')  
%     plot(lags,uci_mp_dtrig_lf3s,'m--')
%     plot(lags,lci_mp_dtrig_lf3s,'m--')     
    plot(lags,uci_mp_dtrig_lf2,'k--')
    plot(lags,lci_mp_dtrig_lf2,'k--')
    xlim([-5 10]);
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_mp_dtrig_' num2str(c)];
    print('-dpng',t_names);
    xlim([-2 2]); grid
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_mp_dtrig_zoom_' num2str(c)];
    print('-dpng',t_names);
    close
   
%% lf8 utrig
%     avg_lf8_utrig_mp(c,:) = mean(lf8_utrig_mp(get_cells,:));
%     uci_lf8_utrig_mp = avg_lf8_utrig_mp(c,:)+tinv(0.975,length(get_cells))*std(lf8_utrig_mp(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_utrig_mp = avg_lf8_utrig_mp(c,:)-tinv(0.975,length(get_cells))*std(lf8_utrig_mp(get_cells,:))/sqrt(length(get_cells));
%     
%     avg_lf8_utrig_lf8(c,:) = mean(lf8_utrig_lf8(get_cells,:));
%     uci_lf8_utrig_lf8 = avg_lf8_utrig_lf8(c,:)+tinv(0.975,length(get_cells))*std(lf8_utrig_lf8(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_utrig_lf8 = avg_lf8_utrig_lf8(c,:)-tinv(0.975,length(get_cells))*std(lf8_utrig_lf8(get_cells,:))/sqrt(length(get_cells));
%     
%     avg_lf8_utrig_lf3(c,:) = mean(lf8_utrig_lf3(get_cells,:));
%     uci_lf8_utrig_lf3 = avg_lf8_utrig_lf3(c,:)+tinv(0.975,length(get_cells))*std(lf8_utrig_lf3(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_utrig_lf3 = avg_lf8_utrig_lf3(c,:)-tinv(0.975,length(get_cells))*std(lf8_utrig_lf3(get_cells,:))/sqrt(length(get_cells));
%     
%     avg_lf8_utrig_lf3s(c,:) = mean(lf8_utrig_lf3s(get_cells,:));
%     uci_lf8_utrig_lf3s = avg_lf8_utrig_lf3s(c,:)+tinv(0.975,length(get_cells))*std(lf8_utrig_lf3s(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_utrig_lf3s = avg_lf8_utrig_lf3s(c,:)-tinv(0.975,length(get_cells))*std(lf8_utrig_lf3s(get_cells,:))/sqrt(length(get_cells));
% 
%     avg_lf8_utrig_lf2(c,:) = mean(lf8_utrig_lf2(get_cells,:));
%     uci_lf8_utrig_lf2 = avg_lf8_utrig_lf2(c,:)+tinv(0.975,length(get_cells))*std(lf8_utrig_lf2(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_utrig_lf2 = avg_lf8_utrig_lf2(c,:)-tinv(0.975,length(get_cells))*std(lf8_utrig_lf2(get_cells,:))/sqrt(length(get_cells));
%  
%     plot(lags,mp_utrig_mp(get_cells,:))
%     xlim([-1 1]);grid
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\mp_up_slope_compare_' num2str(c)];
%     print('-dpng',t_names);
%     close 
%     
%     plot(lags,mp_dtrig_mp(get_cells,:))
%     xlim([-1 1]);grid
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\mp_down_slope_compare_' num2str(c)];
%     print('-dpng',t_names);
%     close
%     
%     plot(lags,lf8_utrig_lf8(get_cells,:))
%     xlim([-1 1]);grid
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\lf8_up_slope_compare_' num2str(c)];
%     print('-dpng',t_names);
%     close
%     
%     plot(lags,lf8_dtrig_lf8(get_cells,:))
%     xlim([-1 1]);grid
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\lf8_down_slope_compare_' num2str(c)];
%     print('-dpng',t_names);
%  close
%  
%     plot(lags,avg_lf8_utrig_mp(c,:),'linewidth',1)
%     hold on
%     plot(lags,avg_lf8_utrig_lf8(c,:),'r','linewidth',1)
% %     plot(lags,avg_lf8_utrig_lf3(c,:),'g','linewidth',1)
% % 	plot(lags,avg_lf8_utrig_lf3s(c,:),'m','linewidth',1)
%     plot(lags,avg_lf8_utrig_lf2(c,:),'k','linewidth',1)
%     legend('MP','LF8','LF2')
%     plot(lags,uci_lf8_utrig_mp,'--')
%     plot(lags,lci_lf8_utrig_mp,'--')
%         plot(lags,uci_lf8_utrig_lf8,'r--')
%     plot(lags,lci_lf8_utrig_lf8,'r--')
% %     plot(lags,uci_lf8_utrig_lf3,'g--')
% %     plot(lags,lci_lf8_utrig_lf3,'g--')    
% %     plot(lags,uci_lf8_utrig_lf3s,'m--')
% %     plot(lags,lci_lf8_utrig_lf3s,'m--')    
%     plot(lags,uci_lf8_utrig_lf2,'k--')
%     plot(lags,lci_lf8_utrig_lf2,'k--')
%     xlim([-5 10]);
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\lf8_utrig_' num2str(c)];
%     print('-dpng',t_names);
%     xlim([-2 2]); grid
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\lf8_utrig_zoom_' num2str(c)];
%     print('-dpng',t_names);
%     close
%  
% %% lf8 dtrig
%     avg_lf8_dtrig_mp(c,:) = mean(lf8_dtrig_mp(get_cells,:));
%     uci_lf8_dtrig_mp = avg_lf8_dtrig_mp(c,:)+tinv(0.975,length(get_cells))*std(lf8_dtrig_mp(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_dtrig_mp = avg_lf8_dtrig_mp(c,:)-tinv(0.975,length(get_cells))*std(lf8_dtrig_mp(get_cells,:))/sqrt(length(get_cells));
%     
%     avg_lf8_dtrig_lf8(c,:) = mean(lf8_dtrig_lf8(get_cells,:));
%     uci_lf8_dtrig_lf8 = avg_lf8_dtrig_lf8(c,:)+tinv(0.975,length(get_cells))*std(lf8_dtrig_lf8(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_dtrig_lf8 = avg_lf8_dtrig_lf8(c,:)-tinv(0.975,length(get_cells))*std(lf8_dtrig_lf8(get_cells,:))/sqrt(length(get_cells));
%     
%     avg_lf8_dtrig_lf3(c,:) = mean(lf8_dtrig_lf3(get_cells,:));
%     uci_lf8_dtrig_lf3 = avg_lf8_dtrig_lf3(c,:)+tinv(0.975,length(get_cells))*std(lf8_dtrig_lf3(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_dtrig_lf3 = avg_lf8_dtrig_lf3(c,:)-tinv(0.975,length(get_cells))*std(lf8_dtrig_lf3(get_cells,:))/sqrt(length(get_cells));
%     
%     avg_lf8_dtrig_lf3s(c,:) = mean(lf8_dtrig_lf3s(get_cells,:));
%     uci_lf8_dtrig_lf3s = avg_lf8_dtrig_lf3s(c,:)+tinv(0.975,length(get_cells))*std(lf8_dtrig_lf3s(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_dtrig_lf3s = avg_lf8_dtrig_lf3s(c,:)-tinv(0.975,length(get_cells))*std(lf8_dtrig_lf3s(get_cells,:))/sqrt(length(get_cells));
% 
%     avg_lf8_dtrig_lf2(c,:) = mean(lf8_dtrig_lf2(get_cells,:));
%     uci_lf8_dtrig_lf2 = avg_lf8_dtrig_lf2(c,:)+tinv(0.975,length(get_cells))*std(lf8_dtrig_lf2(get_cells,:))/sqrt(length(get_cells));
%     lci_lf8_dtrig_lf2 = avg_lf8_dtrig_lf2(c,:)-tinv(0.975,length(get_cells))*std(lf8_dtrig_lf2(get_cells,:))/sqrt(length(get_cells));
%   
%     
%     plot(lags,avg_lf8_dtrig_mp(c,:),'linewidth',1)
%     hold on
%     plot(lags,avg_lf8_dtrig_lf8(c,:),'r','linewidth',1)
% %     plot(lags,avg_lf8_dtrig_lf3(c,:),'g','linewidth',1)
% % 	plot(lags,avg_lf8_dtrig_lf3s(c,:),'m','linewidth',1)
%     plot(lags,avg_lf8_dtrig_lf2(c,:),'k','linewidth',1)
%     legend('MP','LF8','LF2')
%     plot(lags,uci_lf8_dtrig_mp,'--')
%     plot(lags,lci_lf8_dtrig_mp,'--')
%         plot(lags,uci_lf8_dtrig_lf8,'r--')
%     plot(lags,lci_lf8_dtrig_lf8,'r--')
% %     plot(lags,uci_lf8_dtrig_lf3,'g--')
% %     plot(lags,lci_lf8_dtrig_lf3,'g--')
% %        plot(lags,uci_lf8_dtrig_lf3s,'m--')
% %     plot(lags,lci_lf8_dtrig_lf3s,'m--')
%  
%     plot(lags,uci_lf8_dtrig_lf2,'k--')
%     plot(lags,lci_lf8_dtrig_lf2,'k--')
%     xlim([-5 10]);
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\lf8_dtrig_' num2str(c)];
%     print('-dpng',t_names);
%     xlim([-2 2]); grid
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\lf8_dtrig_zoom_' num2str(c)];
%     print('-dpng',t_names);
%     close
 
end
 
 
%% Plot direct comparisons
%MP utrig MP
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags,avg_mp_utrig_mp(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-5 10])
    legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_utrig_mp'];
    print('-dpng',t_names);
    xlim([-2 2]);grid on
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_utrig_mp_zoom'];
    print('-dpng',t_names);
    close
 
%MP utrig LF8
for i = 1:length(diff_cell_types)
   plot(lags,avg_mp_utrig_lf8(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-5 10])
    legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_utrig_lf8'];
    print('-dpng',t_names);
    xlim([-2 2]);grid on
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_utrig_lf8_zoom'];
    print('-dpng',t_names);
    close
 
%MP utrig LF3
for i = 1:length(diff_cell_types)
   plot(lags,avg_mp_utrig_lf3(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-5 10])
    legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_utrig_lf3'];
    print('-dpng',t_names);
    xlim([-2 2]);grid on
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_utrig_lf3_zoom'];
    print('-dpng',t_names);
    close
  
% %MP utrig LF3s
% for i = 1:length(diff_cell_types)
%    plot(lags,avg_mp_utrig_lf3s(i,:),'color',cmap(i,:),'linewidth',1) 
%    hold on
% end
%     xlim([-5 10])
%     legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_mp_utrig_lf3s'];
%     print('-dpng',t_names);
%     xlim([-2 2]);grid on
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_mp_utrig_lf3s_zoom'];
%     print('-dpng',t_names);
%     close
 
%MP utrig LF2
for i = 1:length(diff_cell_types)
   plot(lags,avg_mp_utrig_lf2(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-5 10])
    legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_utrig_lf2'];
    print('-dpng',t_names);
    xlim([-2 2]);grid on
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_utrig_lf2_zoom'];
    print('-dpng',t_names);
    close
    
%MP dtrig MP
cmap = colormap(jet(length(diff_cell_types)));
figure
for i = 1:length(diff_cell_types)
   plot(lags,avg_mp_dtrig_mp(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-5 10])
    legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_dtrig_mp'];
    print('-dpng',t_names);
    xlim([-2 2]);grid on
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_dtrig_mp_zoom'];
    print('-dpng',t_names);
    close
 
%MP dtrig LF8
for i = 1:length(diff_cell_types)
   plot(lags,avg_mp_dtrig_lf8(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-5 10])
    legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_dtrig_lf8'];
    print('-dpng',t_names);
    xlim([-2 2]);grid on
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_dtrig_lf8_zoom'];
    print('-dpng',t_names);
    close
 
%MP dtrig LF3
for i = 1:length(diff_cell_types)
   plot(lags,avg_mp_dtrig_lf3(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-5 10])
    legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_dtrig_lf3'];
    print('-dpng',t_names);
    xlim([-2 2]);grid on
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_dtrig_lf3_zoom'];
    print('-dpng',t_names);
    close
   
%MP dtrig LF2
for i = 1:length(diff_cell_types)
   plot(lags,avg_mp_utrig_lf2(i,:),'color',cmap(i,:),'linewidth',1) 
   hold on
end
    xlim([-5 10])
    legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_dtrig_lf2'];
    print('-dpng',t_names);
    xlim([-2 2]);grid on
    t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\abs_overall_mp_dtrig_lf2_zoom'];
    print('-dpng',t_names);
    close
    
 
    
%     
% %LF8 utrig MP
% cmap = colormap(jet(length(diff_cell_types)));
% figure
% for i = 1:length(diff_cell_types)
%    plot(lags,avg_lf8_utrig_mp(i,:),'color',cmap(i,:),'linewidth',1) 
%    hold on
% end
%     xlim([-5 10])
%     legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_utrig_mp'];
%     print('-dpng',t_names);
%     xlim([-2 2]);grid on
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_utrig_mp_zoom'];
%     print('-dpng',t_names);
%     close
%  
% %LF8 utrig LF8
% for i = 1:length(diff_cell_types)
%    plot(lags,avg_lf8_utrig_lf8(i,:),'color',cmap(i,:),'linewidth',1) 
%    hold on
% end
%     xlim([-5 10])
%     legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_utrig_lf8'];
%     print('-dpng',t_names);
%     xlim([-2 2]);grid on
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_utrig_lf8_zoom'];
%     print('-dpng',t_names);
%     close
%  
% %LF8 utrig LF3
% for i = 1:length(diff_cell_types)
%    plot(lags,avg_lf8_utrig_lf3(i,:),'color',cmap(i,:),'linewidth',1) 
%    hold on
% end
%     xlim([-5 10])
%     legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_utrig_lf3'];
%     print('-dpng',t_names);
%     xlim([-2 2]);grid on
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_utrig_lf3_zoom'];
%     print('-dpng',t_names);
%     close
%    
% %LF8 utrig LF2
% for i = 1:length(diff_cell_types)
%    plot(lags,avg_lf8_utrig_lf2(i,:),'color',cmap(i,:),'linewidth',1) 
%    hold on
% end
%     xlim([-5 10])
%     legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_utrig_lf2'];
%     print('-dpng',t_names);
%     xlim([-2 2]);grid on
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_utrig_lf2_zoom'];
%     print('-dpng',t_names);
%     close
%     
% %LF8 dtrig MP
% cmap = colormap(jet(length(diff_cell_types)));
% figure
% for i = 1:length(diff_cell_types)
%    plot(lags,avg_lf8_dtrig_mp(i,:),'color',cmap(i,:),'linewidth',1) 
%    hold on
% end
%     xlim([-5 10])
%     legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_dtrig_mp'];
%     print('-dpng',t_names);
%     xlim([-2 2]);grid on
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_dtrig_mp_zoom'];
%     print('-dpng',t_names);
%     close
%  
% %LF8 dtrig LF8
% for i = 1:length(diff_cell_types)
%    plot(lags,avg_lf8_dtrig_lf8(i,:),'color',cmap(i,:),'linewidth',1) 
%    hold on
% end
%     xlim([-5 10])
%     legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_dtrig_lf8'];
%     print('-dpng',t_names);
%     xlim([-2 2]);grid on
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_dtrig_lf8_zoom'];
%     print('-dpng',t_names);
%     close
%  
% %LF8 dtrig LF3
% for i = 1:length(diff_cell_types)
%    plot(lags,avg_lf8_dtrig_lf3(i,:),'color',cmap(i,:),'linewidth',1) 
%    hold on
% end
%     xlim([-5 10])
%     legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_dtrig_lf3'];
%     print('-dpng',t_names);
%     xlim([-2 2]);grid on
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_dtrig_lf3_zoom'];
%     print('-dpng',t_names);
%     close
%    
% %LF8 dtrig LF2
% for i = 1:length(diff_cell_types)
%    plot(lags,avg_lf8_utrig_lf2(i,:),'color',cmap(i,:),'linewidth',1) 
%    hold on
% end
%     xlim([-5 10])
%     legend('MEC L3 Pyr','MEC L2 stell','MEC L5','LEC L3','LEC L2','LEC ?','MEC L3 pyr atr','MEC L2 ste atr')
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_dtrig_lf2'];
%     print('-dpng',t_names);
%     xlim([-2 2]);grid on
%     t_names = ['C:\WC_Germany\overall_calcs\avg_trig_compare\overall_lf8_dtrig_lf2_zoom'];
%     print('-dpng',t_names);
%     close
%  
%     save C:\WC_Germany\overall_calcs\avg_trig_compare\avg_trig_data_v2 mp* lf8*
