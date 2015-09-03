load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


cd 'E:\WC_Germany\JMM_Analysis\multLFP_all\depthANDDVAxis_plots'


pairIDs = [15,17,23,25,35,38,41,43];

depth(pairIDs) = [];
DVAxis(pairIDs) = [];
med_corr2_off(14) = mean([med_corr2_off(14) med_corr2_off(15)]);
med_corr2_off(16) = mean([med_corr2_off(16) med_corr2_off(17)]);
med_corr2_off(22) = mean([med_corr2_off(22) med_corr2_off(23)]);
med_corr2_off(24) = mean([med_corr2_off(24) med_corr2_off(25)]);
med_corr2_off(34) = mean([med_corr2_off(34) med_corr2_off(34)]);
med_corr2_off(37) = mean([med_corr2_off(37) med_corr2_off(38)]);
med_corr2_off(40) = mean([med_corr2_off(40) med_corr2_off(41)]);
med_corr2_off(42) = mean([med_corr2_off(42) med_corr2_off(43)]);

med_corr3_off(14) = mean([med_corr3_off(14) med_corr3_off(15)]);
med_corr3_off(16) = mean([med_corr3_off(16) med_corr3_off(17)]);
med_corr3_off(22) = mean([med_corr3_off(22) med_corr3_off(23)]);
med_corr3_off(24) = mean([med_corr3_off(24) med_corr3_off(25)]);
med_corr3_off(34) = mean([med_corr3_off(34) med_corr3_off(34)]);
med_corr3_off(37) = mean([med_corr3_off(37) med_corr3_off(38)]);
med_corr3_off(40) = mean([med_corr3_off(40) med_corr3_off(41)]);
med_corr3_off(42) = mean([med_corr3_off(42) med_corr3_off(43)]);

med_corr5_off(14) = mean([med_corr5_off(14) med_corr5_off(15)]);
med_corr5_off(16) = mean([med_corr5_off(16) med_corr5_off(17)]);
med_corr5_off(22) = mean([med_corr5_off(22) med_corr5_off(23)]);
med_corr5_off(24) = mean([med_corr5_off(24) med_corr5_off(25)]);
med_corr5_off(34) = mean([med_corr5_off(34) med_corr5_off(34)]);
med_corr5_off(37) = mean([med_corr5_off(37) med_corr5_off(38)]);
med_corr5_off(40) = mean([med_corr5_off(40) med_corr5_off(41)]);
med_corr5_off(42) = mean([med_corr5_off(42) med_corr5_off(43)]);

med_corr7_off(14) = mean([med_corr7_off(14) med_corr7_off(15)]);
med_corr7_off(16) = mean([med_corr7_off(16) med_corr7_off(17)]);
med_corr7_off(22) = mean([med_corr7_off(22) med_corr7_off(23)]);
med_corr7_off(24) = mean([med_corr7_off(24) med_corr7_off(25)]);
med_corr7_off(34) = mean([med_corr7_off(34) med_corr7_off(34)]);
med_corr7_off(37) = mean([med_corr7_off(37) med_corr7_off(38)]);
med_corr7_off(40) = mean([med_corr7_off(40) med_corr7_off(41)]);
med_corr7_off(42) = mean([med_corr7_off(42) med_corr7_off(43)]);

med_corr8_off(14) = mean([med_corr8_off(14) med_corr8_off(15)]);
med_corr8_off(16) = mean([med_corr8_off(16) med_corr8_off(17)]);
med_corr8_off(22) = mean([med_corr8_off(22) med_corr8_off(23)]);
med_corr8_off(24) = mean([med_corr8_off(24) med_corr8_off(25)]);
med_corr8_off(34) = mean([med_corr8_off(34) med_corr8_off(34)]);
med_corr8_off(37) = mean([med_corr8_off(37) med_corr8_off(38)]);
med_corr8_off(40) = mean([med_corr8_off(40) med_corr8_off(41)]);
med_corr8_off(42) = mean([med_corr8_off(42) med_corr8_off(43)]);


med_corr2_off(pairIDs) = [];
med_corr3_off(pairIDs) = [];
med_corr5_off(pairIDs) = [];
med_corr7_off(pairIDs) = [];
med_corr8_off(pairIDs) = [];


%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr2_amp(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr2_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr2_amp'];
%     print('-dpng',tname)
%     close all
% 
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr3_amp(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr3_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr3_amp'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr5_amp(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr5_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr5_amp'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr7_amp(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr7_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr7_amp'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr8_amp(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr8_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr8_amp'];
%     print('-dpng',tname)
%     close all
%     
%     %%%%% 
%     
%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr2_amp(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr2_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr2_amp'];
%     print('-dpng',tname)
%     close all
% 
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr3_amp(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr3_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr3_amp'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr5_amp(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr5_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr5_amp'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr7_amp(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr7_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr7_amp'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr8_amp(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr8_amp(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr8_amp'];
%     print('-dpng',tname)
%     close all
%     
%     %%%%%
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr2_coef(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr2_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr2_coef'];
%     print('-dpng',tname)
%     close all
% 
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr3_coef(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr3_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr3_coef'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr5_coef(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr5_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr5_coef'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr7_coef(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr7_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr7_coef'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(depth(1:17),med_corr8_coef(1:17),'o')
%     hold on
%     plot(depth(18:end),med_corr8_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr8_coef'];
%     print('-dpng',tname)
%     close all
%     
%     
%     %%%%%%
%     
%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr2_coef(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr2_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr2_coef'];
%     print('-dpng',tname)
%     close all
% 
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr3_coef(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr3_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr3_coef'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr5_coef(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr5_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr5_coef'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr7_coef(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr7_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr7_coef'];
%     print('-dpng',tname)
%     close all
%     
%         figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(DVAxis(1:17),med_corr8_coef(1:17),'o')
%     hold on
%     plot(DVAxis(18:end),med_corr8_coef(18:end),'ro')
%     legend('Stellate','Pyramidal')
%     xlabel('Depth (um)','FontSize',14);
%     ylabel('Median Corr Amp','FontSize',14)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr8_coef'];
%     print('-dpng',tname)
%     close all
%         
    
    %%%%%
    
    
     
        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(depth(1:15),med_corr2_off(1:15),'o')
    hold on
    plot(depth(16:end),med_corr2_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr2_off'];
    print('-dpng',tname)
    close all

        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(depth(1:15),med_corr3_off(1:15),'o')
    hold on
    plot(depth(16:end),med_corr3_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr3_off'];
    print('-dpng',tname)
    close all
    
        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(depth(1:15),med_corr5_off(1:15),'o')
    hold on
    plot(depth(16:end),med_corr5_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr5_off'];
    print('-dpng',tname)
    close all
    
        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(depth(1:15),med_corr7_off(1:15),'o')
    hold on
    plot(depth(16:end),med_corr7_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr7_off'];
    print('-dpng',tname)
    close all
    
        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(depth(1:15),med_corr8_off(1:15),'o')
    hold on
    plot(depth(16:end),med_corr8_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\depth_corr8_off'];
    print('-dpng',tname)
    close all
    
    
    %%%%%
    
        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(DVAxis(1:15),med_corr2_off(1:15),'o')
    hold on
    plot(DVAxis(16:end),med_corr2_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr2_off'];
    print('-dpng',tname)
    close all

        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(DVAxis(1:15),med_corr3_off(1:15),'o')
    hold on
    plot(DVAxis(16:end),med_corr3_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr3_off'];
    print('-dpng',tname)
    close all
    
        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(DVAxis(1:15),med_corr5_off(1:15),'o')
    hold on
    plot(DVAxis(16:end),med_corr5_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr5_off'];
    print('-dpng',tname)
    close all
    
        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(DVAxis(1:15),med_corr7_off(1:15),'o')
    hold on
    plot(DVAxis(16:end),med_corr7_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr7_off'];
    print('-dpng',tname)
    close all
    
        figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(DVAxis(1:15),med_corr8_off(1:15),'o')
    hold on
    plot(DVAxis(16:end),med_corr8_off(16:end),'ro')
    legend('Stellate','Pyramidal')
    xlabel('Depth (um)','FontSize',14);
    ylabel('Median Corr Amp','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_all\depth_plots\dV_corr8_off'];
    print('-dpng',tname)
    close all