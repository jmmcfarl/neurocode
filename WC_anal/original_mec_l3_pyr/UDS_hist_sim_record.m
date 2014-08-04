cd 'C:\WC_Germany\2007-11-02_Par_EC_LFP_B\2007-11-2_14-36-23'

g = [0:0.1:20];

segDur = 30;
winSlide = 1;
dsf = 8;
Fsd = 2016/dsf;


load C:\WC_Germany\JMM_Analysis_pyr\sim_record\up_per_data2

    load used_data lf8 lf15

    sm_threshold_15 = jmm_smooth_1d(threshold_15,10);
    sm_threshold_8 = jmm_smooth_1d(threshold_8,10);
    
    disp('WCV')
[up_state_dur,down_state_dur,up_trans,down_trans,sig_15,int_thresh_15] = UDS_extract_run_hist(lf15,sm_threshold_15);
disp('LFP')
[up_state_dur8,down_state_dur8,up_trans8,down_trans8,sig_8,int_thresh_8] = UDS_extract_run_hist(lf8,sm_threshold_8);

% dataViewer_run_hist(sig_w{d},sig_8{d},int_thresh_w{d},int_thresh_8{d},up_trans{d},up_trans8{d},down_trans{d},down_trans8{d})


%calculate average up state duration in segmented data
nsamples = length(sig_15);
numWins = floor((nsamples - segDur*Fsd)/(winSlide*Fsd));
for w = 1:numWins
    
    begSamp = round(1+(w-1)*winSlide*Fsd);
    endSamp = begSamp + round(segDur*Fsd);
   
    %find current up states
    curUps_15 = find(up_trans>=begSamp&up_trans<endSamp);
    curUps_8 = find(up_trans8>=begSamp&up_trans8<endSamp);

    curDowns_15 = find(down_trans>=begSamp&down_trans<endSamp);
    curDowns_8 = find(down_trans8>=begSamp&down_trans8<endSamp);
    
    avg_up_dur15(w) = mean(up_state_dur(curUps_15));
    avg_up_dur8(w) = mean(up_state_dur8(curUps_8));
    
    if ~isempty(curUps_15)
        max_up_dur(w) = max(up_state_dur(curUps_15));
    else
        max_up_dur(w) = nan;
    end
    
    if ~isempty(curUps_8)
        max_up_dur8(w) = max(up_state_dur8(curUps_8));
    else
        max_up_dur8(w) = nan;
    end
    
    curDowns_15(curDowns_15==length(down_trans)) = [];
    curDowns_8(curDowns_8==length(down_trans8)) = [];
    
    avg_down_dur(w) = mean(down_state_dur(curDowns_15));
    avg_down_dur8(w) = mean(down_state_dur8(curDowns_8));
    
end

%for each membrane up transition find nearest preceeding lfp up 
for t = 1:length(up_trans)
    
    prev_up = up_trans8(find(up_trans8<up_trans(t),1,'last'));
    if ~isempty(prev_up)
        near_lf8_up(t) = prev_up;
    else
        near_lf8_up(t) = nan;
    end
    
    
end

% up_dist{d} = hist(up_state_dur{d},g);
% up_dist{d} = up_dist{d}/sum(up_dist{d});
% down_dist{d} = hist(down_state_dur{d},g);
% down_dist{d} = down_dist{d}/sum(down_dist{d});
% 
% up_dist8{d} = hist(up_state_dur8{d},g);
% up_dist8{d} = up_dist8{d}/sum(up_dist8{d});
% down_dist8{d} = hist(down_state_dur8{d},g);
% down_dist8{d} = down_dist8{d}/sum(down_dist8{d});
  
[up_dist,up_r] = log_hist(up_state_dur,[0.3 50],100);
[down_dist,down_r] = log_hist(down_state_dur,[0.3 50],100);
[up_dist8,up_r] = log_hist(up_state_dur8,[0.3 50],100);
[down_dist8,down_r]  = log_hist(down_state_dur8,[0.3 50],100);

% Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
% plot(up_state_dur,(up_trans-near_lf8_up)/Fsd,'.')
% xlabel('Up State Duration','FontSize',14)
% ylabel('Up Transition Lag','FontSize',14)
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% [a,b] = corrcoef(up_state_dur(~isnan(near_lf8_up)),(up_trans(~isnan(near_lf8_up))-near_lf8_up(~isnan(near_lf8_up))));
% title(sprintf('Corr = %d   P = %d',a(2,1),b(2,1)))
% tname = 'C:\WC_Germany\JMM_Analysis_pyr\sim_record\sess1_dur_lag_corr';
%     print('-dpng',tname)
% 
% Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
% subplot(2,1,1)
% plot(up_trans{d}/Fsd,(up_trans{d}-near_lf8_up{d})/Fsd,'o')
% xlim([0 length(avg_up_dur{d})])
% set(gca,'yscale','log')
% xlabel('Time (s)','FontSize',14)
% ylabel('Up Transition Lag','FontSize',14)
% subplot(2,1,2)
% plot(max_up_dur{d},'linewidth',2)
% hold on
% plot(max_up_dur8{d},'r','linewidth',2)
% set(gca,'yscale','log')
% xlabel('Time','FontSize',14)
% ylabel('Max Up State Duration','FontSize',14)
% xlim([0 length(max_up_dur{d})])
%     tname = ['C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\dur_lag' f_names{d}];
%     print('-dpng',tname)
    
    
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
% %     max_up = max([prctile(up_state_dur{d},90) prctile(up_state_dur8{d},90)]);
%    stairs(up_r,up_dist)
%     hold on
%     stairs(up_r,up_dist8,'r')
%     xlim([0.2 20])
%     set(gca,'xscale','log')
%     legend('LF15','LF8')
%     title('Up State Duration','FontSize',14)
%     grid on
%     subplot(2,1,2)
%     stairs(down_r,down_dist)
%     hold on
%     stairs(down_r,down_dist8,'r')
% %     max_down = max([prctile(down_state_dur{d},90) prctile(down_state_dur8{d},90)]);
%     xlim([0.2 20])
%         set(gca,'xscale','log')
%         grid on
%     legend('LF15','LF8')
%     title('Down State Duration','FontSize',14)
%     grid on
%     tname = 'C:\WC_Germany\JMM_Analysis_pyr\sim_record\sess2_UDS_hist';
%     print('-dpng',tname)
%     close

    save C:\WC_Germany\JMM_Analysis_pyr\sim_record_UDS_dur_2 avg* up_state_dur* down_state_dur* up_trans* down_trans* 
    


% load C:\WC_Germany\JMM_Analysis_update\meta_data_use_2
% reps = find(diff(cell_ID_map) == 0);
% 
% % merge data for split sessions
% 
% for i = 2:41
% if cell_ID_map(i) == cell_ID_map(i-1)
% up_state_dur{i-1} = [up_state_dur{i-1} up_state_dur{i}];
% up_state_dur8{i-1} = [up_state_dur8{i-1} up_state_dur8{i}];
% down_state_dur{i-1} = [down_state_dur{i-1} down_state_dur{i}];
% down_state_dur8{i-1} = [down_state_dur8{i-1} down_state_dur8{i}];
% end
% end
% 
% up_state_dur(reps) = [];
% up_state_dur8(reps) = [];
% down_state_dur(reps) = [];
% down_state_dur8(reps) = [];
% 
% nbins = 40;
% binrange = [0.05 30];
% for i = 1:33
%     [nup(i,:),x] = log_hist(up_state_dur{i},binrange,nbins);
%     [nup8(i,:),x] = log_hist(up_state_dur8{i},binrange,nbins);
%     [ndown(i,:),x] = log_hist(down_state_dur{i},binrange,nbins);
%     [ndown8(i,:),x] = log_hist(down_state_dur8{i},binrange,nbins);
% end
% 
% % plot overall up state dists
% figure
% cmap = colormap(jet(16));
% subplot(2,1,1)
% for i = 1:16
%     plot(x,nup(i,:),'Color',cmap(i,:))
%     hold on
% end
% set(gca,'xscale','log')
% xlim([0.1 20])
% xlabel('Time (s)')
% title('Stellate Up State Dur','FontSize',14)
% cmap = colormap(jet(17));
% subplot(2,1,2)
% for i = 1:17
%     plot(x,nup(i+16,:),'Color',cmap(i,:))
%     hold on
% end
% set(gca,'xscale','log')
% xlim([0.1 20])
% xlabel('Time (s)')
% title('Pyramidal Up State Dur','FontSize',14)
% 
% % plot overall down state dists
% figure
% cmap = colormap(jet(16));
% subplot(2,1,1)
% for i = 1:16
%     plot(x,ndown(i,:),'Color',cmap(i,:))
%     hold on
% end
% set(gca,'xscale','log')
% xlim([0.1 20])
% xlabel('Time (s)')
% title('Stellate Down State Dur','FontSize',14)
% cmap = colormap(jet(17));
% subplot(2,1,2)
% for i = 1:17
%     plot(x,ndown(i+16,:),'Color',cmap(i,:))
%     hold on
% end
% set(gca,'xscale','log')
% xlim([0.1 20])
% xlabel('Time (s)')
% title('Pyramidal Down State Dur','FontSize',14)
% 
% 
% % get averages and plot
% nup_stells = mean(nup(1:16,:),1);
% nup_pyrs = mean(nup(17:end,:),1);
% ndown_stells = mean(ndown(1:16,:),1);
% ndown_pyrs = mean(ndown(17:end,:),1);
% nup8av = mean(nup8,1);
% ndown8av = mean(ndown8,1);
% 
% figure
% stairs(x,nup_stells,'linewidth',2)
% hold on
% stairs(x,nup_pyrs,'r','linewidth',2)
% stairs(x,nup8av,'k','linewidth',2)
% set(gca,'xscale','log')
% xlabel('Time (s)')
% title('Avg Up State Duration','FontSize',14)
% xlim([0.2 25])
% 
% figure
% stairs(x,ndown_stells,'linewidth',2)
% hold on
% stairs(x,ndown_pyrs,'r','linewidth',2)
% stairs(x,ndown8av,'k','linewidth',2)
% set(gca,'xscale','log')
% xlabel('Time (s)')
% title('Avg Down State Duration','FontSize',14)
% 
