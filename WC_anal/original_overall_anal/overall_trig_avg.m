clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);
dsf = 8;
Fsd = Fs/dsf;

backlag = 5*Fsd;
forwardlag = 10*Fsd;
lags = (-backlag:forwardlag)/Fsd;

for d = 1:length(over_dir)
    
    
    cd(over_dir{d})
    disp(num2str(d))
    
    load used_data lf8 lf2 lf3 wcv_minus_spike
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    lf2_f = filtfilt(b,a,lf2);
    lf3_f = filtfilt(b,a,lf3);
    
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    down_2 = downsample(lf2_f,dsf);
    down_3 = downsample(lf3_f,dsf);
    
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
    down_2 = zscore(down_2);
    
%     temp_cor = corrcoef(down_3,down_8);
%     hc_cor = temp_cor(2,1);
%     down3_sub = down_3-hc_cor*down_8;
    
    %initialize
%     mp_utrig_mp_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf8_mat = zeros(length(synch_ups{d}),length(lags));
%     mp_utrig_lf3_mat = zeros(length(synch_ups{d}),length(lags));
%     mp_utrig_lf2_mat = zeros(length(synch_ups{d}),length(lags));
    
%     lf8_utrig_mp_mat = zeros(length(synch_ups8{d}),length(lags));
%     lf8_utrig_lf8_mat = zeros(length(synch_ups8{d}),length(lags));
%     lf8_utrig_lf3_mat = zeros(length(synch_ups8{d}),length(lags));
%     lf8_utrig_lf2_mat = zeros(length(synch_ups8{d}),length(lags));
%     
%     mp_dtrig_mp_mat = zeros(length(synch_downs{d}),length(lags));
%     mp_dtrig_lf8_mat = zeros(length(synch_downs{d}),length(lags));
%     mp_dtrig_lf3_mat = zeros(length(synch_downs{d}),length(lags));
%     mp_dtrig_lf2_mat = zeros(length(synch_downs{d}),length(lags));
%     
%     lf8_dtrig_mp_mat = zeros(length(synch_downs8{d}),length(lags));
%     lf8_dtrig_lf8_mat = zeros(length(synch_downs8{d}),length(lags));
%     lf8_dtrig_lf3_mat =zeros(length(synch_downs8{d}),length(lags));
%     lf8_dtrig_lf2_mat = zeros(length(synch_downs8{d}),length(lags));
    
    
    %calculate mp utrigs
    for i = 1:length(synch_ups{d})
        
        if up_trans{d}(synch_ups{d}(i)) > backlag && ...
                length(down_w) - up_trans{d}(synch_ups{d}(i)) > forwardlag
           
%             mp_utrig_mp_mat(i,:) = down_w(up_trans{d}(synch_ups{d}(i))-backlag:...
%                 up_trans{d}(synch_ups{d}(i))+forwardlag);
            mp_utrig_lf8_mat(i,:) = down_8(up_trans{d}(synch_ups{d}(i))-backlag:...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
%             mp_utrig_lf3_mat(i,:) = down_3(up_trans{d}(synch_ups{d}(i))-backlag:...
%                 up_trans{d}(synch_ups{d}(i))+forwardlag);
%             mp_utrig_lf2_mat(i,:) = down_2(up_trans{d}(synch_ups{d}(i))-backlag:...
%                 up_trans{d}(synch_ups{d}(i))+forwardlag);

        else
            
%             mp_utrig_mp_mat(i,:) = nan;
            mp_utrig_lf8_mat(i,:) = nan;
%             mp_utrig_lf3_mat(i,:) = nan;
%             mp_utrig_lf2_mat(i,:) = nan;
            
        end
        
    end
    
       %calculate mp dtrigs
%     for i = 1:length(synch_downs{d})
%         
%         if down_trans{d}(synch_downs{d}(i)) > backlag && ...
%                 length(down_w) - down_trans{d}(synch_downs{d}(i)) > forwardlag
%            
%             mp_dtrig_mp_mat(i,:) = down_w(down_trans{d}(synch_downs{d}(i))-backlag:...
%                 down_trans{d}(synch_downs{d}(i))+forwardlag);
%             mp_dtrig_lf8_mat(i,:) = down_8(down_trans{d}(synch_downs{d}(i))-backlag:...
%                 down_trans{d}(synch_downs{d}(i))+forwardlag);
%             mp_dtrig_lf3_mat(i,:) = down_3(down_trans{d}(synch_downs{d}(i))-backlag:...
%                 down_trans{d}(synch_downs{d}(i))+forwardlag);
%             mp_dtrig_lf2_mat(i,:) = down_2(down_trans{d}(synch_downs{d}(i))-backlag:...
%                 down_trans{d}(synch_downs{d}(i))+forwardlag);
% 
%         else
%             
%             mp_dtrig_mp_mat(i,:) = nan;
%             mp_dtrig_lf8_mat(i,:) = nan;
%             mp_dtrig_lf3_mat(i,:) = nan;
%             mp_dtrig_lf2_mat(i,:) = nan;
%             
%         end
%         
%     end
 
      %calculate lfp utrigs
%     for i = 1:length(synch_ups8{d})
%         
%         if up_trans8{d}(synch_ups8{d}(i)) > backlag && ...
%                 length(down_w) - up_trans8{d}(synch_ups8{d}(i)) > forwardlag
% %            
%             lf8_utrig_mp_mat(i,:) = down_w(up_trans8{d}(synch_ups8{d}(i))-backlag:...
%                 up_trans8{d}(synch_ups8{d}(i))+forwardlag);
%             lf8_utrig_lf8_mat(i,:) = down_8(up_trans8{d}(synch_ups8{d}(i))-backlag:...
%                 up_trans8{d}(synch_ups8{d}(i))+forwardlag);
%             lf8_utrig_lf3_mat(i,:) = down_3(up_trans8{d}(synch_ups8{d}(i))-backlag:...
%                 up_trans8{d}(synch_ups8{d}(i))+forwardlag);
%             lf8_utrig_lf2_mat(i,:) = down_2(up_trans8{d}(synch_ups8{d}(i))-backlag:...
%                 up_trans8{d}(synch_ups8{d}(i))+forwardlag);
% 
%         else
%             
%             lf8_utrig_mp_mat(i,:) = nan;
%             lf8_utrig_lf8_mat(i,:) = nan;
%             lf8_utrig_lf3_mat(i,:) = nan;
%             lf8_utrig_lf2_mat(i,:) = nan;
%             
%         end
%         
%     end
%     
%        %calculate mp dtrigs
%     for i = 1:length(synch_downs8{d})
%         
%         if down_trans8{d}(synch_downs8{d}(i)) > backlag && ...
%                 length(down_w) - down_trans8{d}(synch_downs8{d}(i)) > forwardlag
%            
%             lf8_dtrig_mp_mat(i,:) = down_w(down_trans8{d}(synch_downs8{d}(i))-backlag:...
%                 down_trans8{d}(synch_downs8{d}(i))+forwardlag);
%             lf8_dtrig_lf8_mat(i,:) = down_8(down_trans8{d}(synch_downs8{d}(i))-backlag:...
%                 down_trans8{d}(synch_downs8{d}(i))+forwardlag);
%             lf8_dtrig_lf3_mat(i,:) = down_3(down_trans8{d}(synch_downs8{d}(i))-backlag:...
%                 down_trans8{d}(synch_downs8{d}(i))+forwardlag);
%             lf8_dtrig_lf2_mat(i,:) = down_2(down_trans8{d}(synch_downs8{d}(i))-backlag:...
%                 down_trans8{d}(synch_downs8{d}(i))+forwardlag);
% 
%         else
%             
%             lf8_dtrig_mp_mat(i,:) = nan;
%             lf8_dtrig_lf8_mat(i,:) = nan;
%             lf8_dtrig_lf3_mat(i,:) = nan;
%             lf8_dtrig_lf2_mat(i,:) = nan;
%             
%         end
%         
%     end
  
%     mp_utrig_mp(d,:) = nanmean(mp_utrig_mp_mat);
    mp_utrig_lf8(d,:) = nanmean(mp_utrig_lf8_mat);
%     mp_utrig_lf3(d,:) = nanmean(mp_utrig_lf3_mat);
%     mp_utrig_lf2(d,:) = nanmean(mp_utrig_lf2_mat);
%     mp_dtrig_mp(d,:) = nanmean(mp_dtrig_mp_mat);
%     mp_dtrig_lf8(d,:) = nanmean(mp_dtrig_lf8_mat);
%     mp_dtrig_lf3(d,:) = nanmean(mp_dtrig_lf3_mat);
%     mp_dtrig_lf2(d,:) = nanmean(mp_dtrig_lf2_mat);
%     lf8_utrig_mp(d,:) = nanmean(lf8_utrig_mp_mat);
%     lf8_utrig_lf8(d,:) = nanmean(lf8_utrig_lf8_mat);
%     lf8_utrig_lf3(d,:) = nanmean(lf8_utrig_lf3_mat);
%     lf8_utrig_lf2(d,:) = nanmean(lf8_utrig_lf2_mat);
%     lf8_dtrig_mp(d,:) = nanmean(lf8_dtrig_mp_mat);
%     lf8_dtrig_lf8(d,:) = nanmean(lf8_dtrig_lf8_mat);
%     lf8_dtrig_lf3(d,:) = nanmean(lf8_dtrig_lf3_mat);
%     lf8_dtrig_lf2(d,:) = nanmean(lf8_dtrig_lf2_mat);
    
%     [dummy,up_order] = sort(up_state_dur{d}(synch_ups{d}));
%     [dummy,down_order] = sort(down_state_dur{d}(synch_downs{d}));
%     [dummy,up_order8] = sort(up_state_dur8{d}(synch_ups8{d}));
%     [dummy,down_order8] = sort(down_state_dur8{d}(synch_downs8{d}));
    
% %% plot mp trig averages    
%    plot(lags,mp_utrig_mp(d,:),'linewidth',2)
%    hold on
%    plot(lags,mp_utrig_lf8(d,:),'r','linewidth',2)
%    plot(lags,mp_utrig_lf3(d,:),'g','linewidth',2)
%    plot(lags,mp_utrig_lf2(d,:),'k','linewidth',2)
%    legend('MP','LF8','LF3','LF2')
%    title('MP Up Triggered Avg')
%    t_names = ['C:\WC_Germany\overall_calcs\trig_avgs\mp_utrig_lf3orth' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
%    xlim([-2 2]); grid on
%       t_names = ['C:\WC_Germany\overall_calcs\trig_avgs\mp_utrig_lf3orth_zoom_' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
% 
% close
%     
%    plot(lags,mp_dtrig_mp(d,:),'linewidth',2)
%    hold on
%    plot(lags,mp_dtrig_lf8(d,:),'r','linewidth',2)
%    plot(lags,mp_dtrig_lf3(d,:),'g','linewidth',2)
%    plot(lags,mp_dtrig_lf2(d,:),'k','linewidth',2)
%    legend('MP','LF8','LF3','LF2')
%    title('MP Down Triggered Avg')
%    t_names = ['C:\WC_Germany\overall_calcs\trig_avgs\mp_dtrig_lf3orth_' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
%       xlim([-2 2]); grid on
%       t_names = ['C:\WC_Germany\overall_calcs\trig_avgs\mp_dtrig_lf3orth_zoom_' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
% 
% close

% %% plot lf8 trig averages    
%    plot(lags,lf8_utrig_mp(d,:),'linewidth',2)
%    hold on
%    plot(lags,lf8_utrig_lf8(d,:),'r','linewidth',2)
%    plot(lags,lf8_utrig_lf3(d,:),'g','linewidth',2)
%    plot(lags,lf8_utrig_lf2(d,:),'k','linewidth',2)
%    legend('MP','LF8','LF3','LF2')
%    title('LF8 Up Triggered Avg')
%    t_names = ['C:\WC_Germany\overall_calcs\trig_avgs\lf8_utrig_lf3orth_' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
%       xlim([-2 2]); grid on
%       t_names = ['C:\WC_Germany\overall_calcs\trig_avgs\lf8_utrig_lf3orth_zoom_' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
% 
% close
%     
%    plot(lags,lf8_dtrig_mp(d,:),'linewidth',2)
%    hold on
%    plot(lags,lf8_dtrig_lf8(d,:),'r','linewidth',2)
%    plot(lags,lf8_dtrig_lf3(d,:),'g','linewidth',2)
%    plot(lags,lf8_dtrig_lf2(d,:),'k','linewidth',2)
%    legend('MP','LF8','LF3','LF2')
%    title('LF8 Down Triggered Avg')
%    t_names = ['C:\WC_Germany\overall_calcs\trig_avgs\lf8_dtrig_lf3orth_' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
%       xlim([-2 2]); grid on
%       t_names = ['C:\WC_Germany\overall_calcs\trig_avgs\lf8_dtrig_lf3orth_zoom_' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
% 
% close

%% plot mp up trig matrices
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,mp_utrig_mp_mat(up_order,:));shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%     caxis([-3 3]);colorbar
%     xlim([-2 3])
%     line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\mp_utrig_mp_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
   
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    pcolor(lags,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
        ,mp_utrig_lf8_mat);shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
        xlim([-2 3])
caxis([-3 3]);colorbar
line([0 0],[0 1],'Color','k')
    t_names = ['C:\WC_Germany\overall_calcs\trig_mats\no_ord_mp_utrig_lf8_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

%         Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,mp_utrig_lf3_mat(up_order,:));shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%         xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\mp_utrig_lf3orth_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
%         Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,mp_utrig_lf2_mat(up_order,:));shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%         xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\mp_utrig_lf2_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
% %% plot mp down trig matrices
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_dtrig_mp_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%         xlim([-2 3])
% caxis([-3 3]);colorbar
% t_names = ['C:\WC_Germany\overall_calcs\trig_mats\mp_dtrig_mp_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
%    
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_dtrig_lf8_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%         xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\mp_dtrig_lf8_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
%         Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_dtrig_lf3_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%         xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\mp_dtrig_lf3orth_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
%         Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_dtrig_lf2_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%         xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\mp_dtrig_lf2_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
% %% plot lfp up trig matrices
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_mp_mat(up_order8,:));shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%         xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\lf8_utrig_mp_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
%    
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf8_mat(up_order8,:));shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%             xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\lf8_utrig_lf8_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
%         Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf3_mat(up_order8,:));shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%             xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\lf8_utrig_lf3orth_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
%         Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf2_mat(up_order8,:));shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%             xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\lf8_utrig_lf2_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
% %% plot mp down trig matrices
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_downs8{d}))/length(synch_downs8{d}) ...
%         ,lf8_dtrig_mp_mat(down_order8,:));shading flat; colorbar; hold on
%     plot(down_state_dur8{d}(synch_downs8{d}(down_order8)),(1:length(synch_downs8{d}))/length(synch_downs8{d}),'w','linewidth',2)
%             xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\lf8_dtrig_mp_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
%    
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_downs8{d}))/length(synch_downs8{d}) ...
%         ,lf8_dtrig_lf8_mat(down_order8,:));shading flat; colorbar; hold on
%     plot(down_state_dur8{d}(synch_downs8{d}(down_order8)),(1:length(synch_downs8{d}))/length(synch_downs8{d}),'w','linewidth',2)
%             xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\lf8_dtrig_lf8_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
%         Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_downs8{d}))/length(synch_downs8{d}) ...
%         ,lf8_dtrig_lf3_mat(down_order8,:));shading flat; colorbar; hold on
%     plot(down_state_dur8{d}(synch_downs8{d}(down_order8)),(1:length(synch_downs8{d}))/length(synch_downs8{d}),'w','linewidth',2)
%             xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\lf8_dtrig_lf3orth_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
%         Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(lags,(1:length(synch_downs8{d}))/length(synch_downs8{d}) ...
%         ,lf8_dtrig_lf2_mat(down_order8,:));shading flat; colorbar; hold on
%     plot(down_state_dur8{d}(synch_downs8{d}(down_order8)),(1:length(synch_downs8{d}))/length(synch_downs8{d}),'w','linewidth',2)
%             xlim([-2 3])
% caxis([-3 3]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats\lf8_dtrig_lf2_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close

end

% clear *mat
% save C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_data lags mp* lf8*