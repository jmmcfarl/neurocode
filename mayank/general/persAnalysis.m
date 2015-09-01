%% ANALYZE PERSISTENT STATES BATCH

clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
load E:\WC_Germany\JMM_Analysis\persPoints.mat

%initialize overall variables

d = 1;

pers_ratio = zeros(1,43);

while d <= length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data wcv_minus_spike lf8
%     load used_data_multLFP lf2 lf5
%     load used_data_multLFP_2 lf3

    Fs = 2016;


%     %filter data
%     nyqf = Fs/2;
%     hif = 2/nyqf;
%     lif = 0.05/nyqf;
%     [b,a] = butter(2, [lif hif]);
%     lofwcv = filtfilt(b,a,wcv_minus_spike);
%     loflf8 = filtfilt(b,a,lf8);
%     loflf2 = filtfilt(b,a,lf2);
%     loflf3 = filtfilt(b,a,lf3);
%     loflf5 = filtfilt(b,a,lf5);
%     loflf7 = filtfilt(b,a,lf7);

    %zscore data
%     lofwcv = lofwcv - mean(lofwcv);
%     lofwcv = lofwcv/std(lofwcv);
%     loflf8 = loflf8 - mean(loflf8);
%     loflf8 = loflf8/std(loflf8);
    wcv_z = wcv_minus_spike - mean(wcv_minus_spike);
    wcv_z = wcv_z/std(wcv_z);
    lf8_z = lf8 - mean(lf8);
    lf8_z = lf8_z/std(lf8_z);
%     loflf2 = loflf2 - mean(loflf2);
%     loflf2 = loflf2/std(loflf2);
%     loflf3 = loflf3 - mean(loflf3);
%     loflf3 = loflof3/std(loflf3);
%     loflf5 = loflf5 - mean(loflf5);
%     loflf5 = loflf5/std(loflf5);
%     loflf7 = loflf7 - mean(loflf7);
%     loflf7 = loflf7/std(loflf7);
    
    
    
    z_grid = [-5:.01:5];

%     %separate data into pers and nonpers
%     persOn = persPoints{d};
%     persOff = setdiff([1:length(wcv_z)],persOn);
%     pers_lofwcv = lofwcv(persOn);
%     noPers_lofwcv = lofwcv(persOff);
%     pers_wcv = wcv_z(persOn);
%     noPers_wcv = wcv_z(persOff);
%     pers_lf8 = lf8_z(persOn);
%     noPers_lf8 = lf8_z(persOff);
%     pers_loflf2 = loflf2(persOn);
%     noPers_loflf2 = loflf2(persOff);
%     pers_loflf3 = loflf3(persOn);
%     noPers_loflf3 = loflf3(persOff);
%     pers_loflf5 = loflf5(persOn);
%     noPers_loflf5 = loflf5(persOff);
%     pers_loflf7 = loflf7(persOn);
%     noPers_loflf7 = loflf7(persOff);
%     pers_loflf8 = loflf8(persOn);
%     noPers_loflf8 = loflf8(persOff);

%     %calculate histogram of pers and non pers
%     if ~isempty(persOn)
%         persHist{d} = hist(pers_wcv,z_grid);
%         persHist{d} = persHist{d}/sum(persHist{d});
%         persHist_LF8{d} = hist(pers_lf8,z_grid);
%         persHist_LF8{d} = persHist_LF8{d}/sum(persHist_LF8{d});
%     end
%     if ~isempty(persOff)
%         noPersHist{d} = hist(noPers_wcv,z_grid);
%         noPersHist{d} = noPersHist{d}/sum(noPersHist{d});
%         noPersHist_LF8{d} = hist(noPers_lf8,z_grid);
%         noPersHist_LF8{d} = noPersHist_LF8{d}/sum(noPersHist_LF8{d});
%     end
% 
%     %calculate ratio of time spent in pers state
%     pers_ratio(d) = round(100*length(persOn)/length(wcv_z));


wcv_hist = hist(wcv_z,z_grid);
lf8_hist = hist(lf8_z,z_grid);

  figure
  set(gcf,'PaperUnits','centimeters');
  set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
  set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
  hold on
  
% if ~isempty(persOn)
%     plot(z_grid,fgsmooth(persHist{d},3),'linewidth',2)
%     plot(z_grid,fgsmooth(persHist_LF8{d},3),'b-')
% else
%     plot(z_grid,zeros(size(z_grid)),'linewidth',2)
%     plot(z_grid,zeros(size(z_grid)),'b-')
% end
% if ~isempty(persOff)
%     plot(z_grid,fgsmooth(noPersHist{d},3),'r','linewidth',2)
%     plot(z_grid,fgsmooth(noPersHist_LF8{d},3),'r-')
% else
%     plot(z_grid,zeros(size(z_grid)),'r','linewidth',2)
%     plot(z_grid,zeros(size(z_grid)),'r-')
% end
% legend('Persistent Act State WCV','Persistent Act State LF8', ... 
%     'Non-Persistent Act State WCV','Non-Persistent Act State')
% xlabel('Z Score','FontSize',14)
% ylabel('Probability','FontSize',14)
% title(sprintf('%d PerCent Persistent Act',pers_ratio(d)),'FontSize',16)
% tname = ['E:\WC_Germany\JMM_Analysis\persistentAct\state_dep_hist_' f_names{d}];
% print('-dpng',tname)
% close all

plot(z_grid,fgsmooth(lf8_hist,3),'linewidth',2)
plot(z_grid,fgsmooth(wcv_hist,3),'r','linewidth',2)
legend('LF8','WCV')
 xlabel('Z Score','FontSize',14)
ylabel('Probability','FontSize',14)
   tname = ['E:\WC_Germany\JMM_Analysis\persistentAct\raw_hist_' f_names{d}];
    print('-dpng',tname)
    close all
    
    
%     %calculate xcorr in segments for pers and nonpers separately
%     windowSize = 50; %window size in s
%     windowSizeBins = round(windowSize*Fs);
%     timeSlide = 10; %amount to slide time window in s
%     binSlide = round(timeSlide*Fs);
%     numWinsPers = round((length(pers_wcv)-windowSizeBins)/binSlide); %number of windows
%     numWinsNoPers = round((length(noPers_wcv) - windowSizeBins)/binSlide);
%     maxLag = 5*Fs; %max lag in sec

    %if there is some pers activity
%     if ~isempty(persOn) & numWinsPers > 0
%         
%         for i = 1:numWinsPers
% 
%             windowCent = (i-1)*binSlide+round(windowSizeBins/2);
%             windowBeg = windowCent - round(windowSizeBins/2)+1;
%             windowEnd = min([(windowCent + round(windowSizeBins/2)) length(persOn)]);
% 
%             wcvSeg = pers_lofwcv(windowBeg:windowEnd);
% %             lf2Seg = pers_loflf2(windowBeg:windowEnd);
% %             lf3Seg = pers_loflf3(windowBeg:windowEnd);
% %             lf5Seg = pers_loflf5(windowBeg:windowEnd);
%             %         lf7Seg = loflf7(windowBeg:windowEnd);
%             lf8Seg = pers_loflf8(windowBeg:windowEnd);
% 
% %             [temp_pxcorr_LF2(i,:),lagVec] = xcov(wcvSeg,lf2Seg,maxLag,'coef');
% %             [temp_pxcorr_LF3(i,:),lagVec] = xcov(wcvSeg,lf3Seg,maxLag,'coef');
% %             [temp_pxcorr_LF5(i,:),lagVec] = xcov(wcvSeg,lf5Seg,maxLag,'coef');
%             %         [temp_xcorr_LF7(i,:),lagVec] = xcov(wcvSeg,lf7Seg,maxLag,'coef');
%             [temp_pxcorr_LF8(i,:),lagVec] = xcov(wcvSeg,lf8Seg,maxLag,'coef');
% 
%         end
% 
% %         pxcorr_LF2(d,:) = mean(temp_pxcorr_LF2);
% %         pxcorr_LF3(d,:) = mean(temp_pxcorr_LF3);
% %         pxcorr_LF5(d,:) = mean(temp_pxcorr_LF5);
%         pxcorr_LF8{d} = mean(temp_pxcorr_LF8);
%         
%     else
%         pxcorr_LF8{d} = [];
%     end
% 
%     %if there is some nonpersistent activity
%        if ~isempty(persOff) & numWinsNoPers > 0
%         
%         for i = 1:numWinsNoPers
% 
%             windowCent = (i-1)*binSlide+round(windowSizeBins/2);
%             windowBeg = windowCent - round(windowSizeBins/2)+1;
%             windowEnd = min([(windowCent + round(windowSizeBins/2)) length(persOff)]);
% 
%             wcvSeg = noPers_lofwcv(windowBeg:windowEnd);
% %             lf2Seg = noPers_loflf2(windowBeg:windowEnd);
% %             lf3Seg = noPers_loflf3(windowBeg:windowEnd);
% %             lf5Seg = noPers_loflf5(windowBeg:windowEnd);
%             %         lf7Seg = loflf7(windowBeg:windowEnd);
%             lf8Seg = noPers_loflf8(windowBeg:windowEnd);
% 
% %             [temp_nxcorr_LF2(i,:),lagVec] = xcov(wcvSeg,lf2Seg,maxLag,'coef');
% %             [temp_nxcorr_LF3(i,:),lagVec] = xcov(wcvSeg,lf3Seg,maxLag,'coef');
% %             [temp_nxcorr_LF5(i,:),lagVec] = xcov(wcvSeg,lf5Seg,maxLag,'coef');
% %             %         [temp_xcorr_LF7(i,:),lagVec] = xcov(wcvSeg,lf7Seg,maxLag,'coef');
%             [temp_nxcorr_LF8(i,:),lagVec] = xcov(wcvSeg,lf8Seg,maxLag,'coef');
% 
%         end
% % 
% %         nxcorr_LF2(d,:) = mean(temp_nxcorr_LF2);
% %         nxcorr_LF3(d,:) = mean(temp_nxcorr_LF3);
% %         nxcorr_LF5(d,:) = mean(temp_nxcorr_LF5);
%         nxcorr_LF8{d} = mean(temp_nxcorr_LF8);
%         
%        else
%            nxcorr_LF8{d} = [];
%        end
    
       
       
       
%          figure
%   set(gcf,'PaperUnits','centimeters');
%   set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%   set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%   hold on
% if ~isempty(persOn) & numWinsPers > 0
%     plot(lagVec/Fs,pxcorr_LF8{d},'linewidth',2)
% else
%     plot(lagVec/Fs,zeros(size(lagVec)),'linewidth',2)
% end
% if ~isempty(persOff) & numWinsNoPers > 0
%     plot(lagVec/Fs,nxcorr_LF8{d},'r','linewidth',2)
% else
%     plot(lagVec/Fs,zeros(size(lagVec)),'r','linewidth',2)
% end
% legend('Persistent Activity','Non-Persistent Activity')
% xlabel('Time Lag (s)','FontSize',14)
% ylabel('XCorr','FontSize',14)
% title(sprintf('%d PerCent Persistent Act',pers_ratio(d)),'FontSize',16)
% tname = ['E:\WC_Germany\JMM_Analysis\persistentAct\lf8_state_dep_xcorr_' f_names{d}];
% print('-dpng',tname)
% close all
% 
%           figure
%   set(gcf,'PaperUnits','centimeters');
%   set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%   set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%   hold on
% if ~isempty(persOn) & numWinsPers > 0
%     plot(lagVec/Fs,pxcorr_LF8{d},'linewidth',2)
% else
%     plot(lagVec/Fs,zeros(size(lagVec)),'linewidth',2)
% end
% if ~isempty(persOff) & numWinsNoPers > 0
%     plot(lagVec/Fs,nxcorr_LF8{d},'r','linewidth',2)
% else
%     plot(lagVec/Fs,zeros(size(lagVec)),'r','linewidth',2)
% end
% legend('Persistent Activity','Non-Persistent Activity')
% xlabel('Time Lag (s)','FontSize',14)
% ylabel('XCorr','FontSize',14)
% xlim([-0.5 1.5])
% grid on
% title(sprintf('%d PerCent Persistent Act',pers_ratio(d)),'FontSize',16)
% tname = ['E:\WC_Germany\JMM_Analysis\persistentAct\lf8_state_dep_xcorr_zoom' f_names{d}];
% print('-dpng',tname)
% close all
      
       
%   clear lf8 lf8Seg lf8_z loflf8 lofwcv noPers_lf8 noPers_loflf8 noPers_lofwcv noPers_wcv ... 
%       numWins* persOff persOn pers_lf8 pers_loflf8 pers_lofwcv pers_wcv ...
%       temp* wcv_minus_spike wcv_z window* wcvSeg
  
clear lf8 lf8_z loflf8 lofwcv wcv_z wcv_minus_spike

  pack       

    d = d +1


end


%     save E:\WC_Germany\JMM_Analysis\overall_data_jmm_pers_analysis lagVec noPersHist* ...
%         nxcorr* persHist* pers_ratio pxcorr*
