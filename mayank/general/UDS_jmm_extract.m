%extract UDS properties

clear all
load hip_jmm_calc


if ~exist('./jmm')
    !mkdir jmm
end

%% extract sliding correlation between flp and wcv
% windowSize = 20; %window size in s
% windowSizeBins = round(windowSize/dt);
% timeSlide = 0.5; %amount to slide time window in s
% binSlide = round(timeSlide/dt);
% numCorrCalcs = round((length(loflf8)-windowSizeBins)/binSlide);
% corrFun = zeros(1,numCorrCalcs);
% 
% for i = 1:numCorrCalcs
%     windowCent = (i-1)*binSlide+round(windowSizeBins/2);
%     windowBeg = windowCent - round(windowSizeBins/2)+1;
%     windowEnd = windowCent + round(windowSizeBins/2);
%     temp = corrcoef(loflf8(windowBeg:windowEnd),lofwcv(windowBeg:windowEnd));
%     corrFun(i) = temp(2,1);
% end
% 
% corrTimeAxis = [1:numCorrCalcs]*timeSlide;
% 
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(corrTimeAxis,corrFun,'k','linewidth',2)
% xlabel('Time (s)','FontSize',14);
% ylabel('Correlation Coefficient','FontSize',14)
% title(sprintf('LFP WCV Correlation Vs Time %d s windows',windowSize),'FontSize',16)
% print('-dpng','./jmm/lfp_wcv_corr_v_time')
% close all


%% Calculate thresholded duty cycle, cross covariance peak amp and offset
clear *duty*
windowSize = 20; %window size in s
windowSizeBins = round(windowSize/dt);
timeSlide = 1; %amount to slide time window in s
binSlide = round(timeSlide/dt);
numDutyCalcs = round((length(loflf8)-windowSizeBins)/binSlide);
dutyFun_lf8 = zeros(1,numDutyCalcs);
dutyFun_wcv = dutyFun_lf8;
corrAmpFun = dutyFun_lf8;
corrOffFun = dutyFun_lf8;

maxLag = round(1.5/dt); %max lag in bins

for i = 1:numDutyCalcs
    windowCent = (i-1)*binSlide+round(windowSizeBins/2);
    windowBeg = windowCent - round(windowSizeBins/2)+1;
    windowEnd = windowCent + round(windowSizeBins/2);
    seg_lf8 = loflf8z(windowBeg:windowEnd);
    seg_wcv = lofwcvz(windowBeg:windowEnd);
    %detrend data segments
    seg_lf8 = detrend(seg_lf8);
    seg_wcv = detrend(seg_wcv);
    
    %find ratio of time spent above mean to time spent below mean
%     dutyFun_lf8(i) = length(find(seg_lf8>0))/length(find(seg_lf8<0));
%     dutyFun_wcv(i) = length(find(seg_wcv>0))/length(find(seg_wcv<0));
%     dutyFun_ratio(i) = dutyFun_wcv(i)/dutyFun_lf8(i);
    
    %find cross corr peak amp and offset
    [temp,crossCorWin] = xcov(seg_lf8,seg_wcv,maxLag,'coef');
    [peakAmp,peakLoc] = max(temp);
    corrAmpFun(i) = peakAmp;
    corrOffFun(i) = (peakLoc-length(crossCorWin)/2)*dt;
    i
end

dutyTimeAxis = [1:numDutyCalcs]*timeSlide;

% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(dutyTimeAxis,fgsmooth(dutyFun_lf8,3),'linewidth',2)
% hold on
% plot(dutyTimeAxis,fgsmooth(dutyFun_wcv,3),'r','linewidth',2)
% plot(dutyTimeAxis,fgsmooth(dutyFun_ratio,3),'k','linewidth',2)
% xlabel('Time (s)','FontSize',14)
% ylabel('Duty Cycle ratio','FontSize',14)
% title('Duty Cycle vs Time','FontSize',16)
% legend('LFP','WCV','WCV/LFP')
% print('-dpng','./jmm/thresh_duty_cycle')
close all


%% calculate duration of each up state
% 
% lf8_up_dur = zeros(1,length(lf8_up_start_points));
% lf8_down_dur = lf8_up_dur;
% for i = 1:length(lf8_up_start_points)-1
%     lf8_up_dur(i) = (lf8_up_end_points(i)-lf8_up_start_points(i))*dt;
%     lf8_down_dur(i) = (lf8_up_start_points(i+1) - lf8_up_end_points(i))*dt;
% end
% smooth_lf8_up_dur = fgsmooth(lf8_up_dur,2);
% smooth_lf8_down_dur = fgsmooth(lf8_down_dur,2);
% lf8_smooth_duty = smooth_lf8_up_dur./smooth_lf8_down_dur;
% 
% wcv_up_dur = zeros(1,length(wcv_up_start_points));
% wcv_down_dur = wcv_up_dur;
% for i = 1:length(wcv_up_start_points)-1
%     wcv_up_dur(i) = (wcv_up_end_points(i)-wcv_up_start_points(i))*dt;
%     wcv_down_dur(i) = (wcv_up_start_points(i+1) - wcv_up_end_points(i))*dt;
% end
% 
% if length(lf8_up_dur) > length(lf8_up_start_time)
%     lf8_up_dur(length(lf8_up_start_time)+1:end) = [];
% end
% if length(lf8_down_dur) > length(lf8_up_end_time)
%     lf8_down_dur(length(lf8_up_end_time)+1:end) = [];
% end
% if length(wcv_up_dur) > length(wcv_up_start_time)
%     wcv_up_dur(length(wcv_up_start_time)+1:end) = [];
% end
% if length(wcv_down_dur) > length(wcv_up_end_time)
%     wcv_down_dur(length(wcv_up_end_time)+1:end) = [];
% end
% 
% 
% smooth_wcv_up_dur = fgsmooth(wcv_up_dur,2);
% smooth_wcv_down_dur = fgsmooth(wcv_down_dur,2);
% wcv_smooth_duty = smooth_wcv_up_dur./smooth_wcv_down_dur;
% 
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(lf8_up_start_time,smooth_lf8_up_dur,'linewidth',2)
% hold on
% plot(lf8_up_end_time,smooth_lf8_down_dur,'r','linewidth',2)
% legend('Up State','Down State')
% xlabel('Time (s)','FontSize',14)
% ylabel('State Duration (s)','FontSize',14)
% title('LFP Up and Down State Duration','FontSize',16)
% print('-dpng','./jmm/lf8_UDS_duration')
% close all
% 
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(wcv_up_start_time,smooth_wcv_up_dur,'linewidth',2)
% hold on
% plot(wcv_up_end_time,smooth_wcv_down_dur,'r','linewidth',2)
% legend('Up State','Down State')
% xlabel('Time (s)','FontSize',14)
% ylabel('State Duration (s)','FontSize',14)
% title('WCV Up and Down State Duration','FontSize',16)
% print('-dpng','./jmm/wcv_UDS_duration')
% close all
% 
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(lf8_up_start_time,smooth_lf8_up_dur,'linewidth',2)
% hold on
% plot(wcv_up_start_time,smooth_wcv_up_dur,'r','linewidth',2)
% legend('LF8 Up State','WCV Up State')
% xlabel('Time (s)','FontSize',14)
% ylabel('State Duration (s)','FontSize',14)
% title('LF8 and WCV Up State Duration','FontSize',16)
% print('-dpng','./jmm/lf8_wcv_up_dur')
% close all
% 
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(lf8_up_end_time,smooth_lf8_down_dur,'linewidth',2)
% hold on
% plot(wcv_up_end_time,smooth_wcv_down_dur,'r','linewidth',2)
% legend('LF8 Down State','WCV Down State')
% xlabel('Time (s)','FontSize',14)
% ylabel('State Duration (s)','FontSize',14)
% title('LF8 and WCV Down State Duration','FontSize',16)
% print('-dpng','./jmm/lf8_wcv_down_dur')
% close all
% 
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(lf8_up_start_time,lf8_smooth_duty,'linewidth',2)
% hold on
% plot(wcv_up_start_time,wcv_smooth_duty,'r','linewidth',2)
% legend('LF8','WCV')
% xlabel('Time (s)','FontSize',14)
% ylabel('Duty Cycle','FontSize',14)
% title('LF8 and WCV Duty Cycle vs Time','FontSize',16)
% print('-dpng','./jmm/lf8_wcv_dutycycle')
% close all
% 
%% find delay to wcv up trans from lf8 up trans
% wcv_up_delay = zeros(size(lf8_up_start_time));
% lf8_trans_cnt = 1;
% while lf8_trans_cnt < length(lf8_up_start_time)
% 
%     %find closest wcv up transition
%     [temp_delay,closest_wcv_up] = min(abs(wcv_up_start_time - lf8_up_start_time(lf8_trans_cnt)));
% 
%     %if delay is more than a second ignore and move to next transition
%     if temp_delay > 1
%         wcv_up_delay(lf8_trans_cnt) = NaN;
%         lf8_trans_cnt = lf8_trans_cnt + 1;
%     else
% 
%         wcv_up_delay(lf8_trans_cnt) = wcv_up_start_time(closest_wcv_up) - lf8_up_start_time(lf8_trans_cnt);
%         
%         %any lf8 up transitions that occur before the wcv comes back down should
%         %be ignored
%         temp = find(lf8_up_start_time > wcv_up_end_time(closest_wcv_up),1,'first');
%         if temp == lf8_trans_cnt
%             lf8_trans_cnt = lf8_trans_cnt + 1
%         else
%             lf8_trans_cnt = temp
%         end
%         
%     end
% 
% end
% wcv_up_delay(wcv_up_delay == 0) = NaN;
% 
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(lf8_up_start_time,wcv_up_delay,'o')
% xlabel('Time (s)','FontSize',14)
% ylabel('WCV up transition Delay (s)','FontSize',14)
% print('-dpng','./jmm/wcv_up_delay')
% close all
% % smooth_wcv_up_delay = fgsmooth(wcv_up_delay,2);
% % plot(lf8_up_start_time,smooth_wcv_up_delay,'linewidth',2)

%% fit slope to up state

% %first need to filter out high frequencies
f=CSC8_SampleFrequencies(1);
lof = .1; hif = 5; nyqf = f/2;
lof = lof/nyqf; hif = hif/nyqf;
[b,a] = butter(2, [lof hif]);

loflf8_1 = filtfilt(b,a,lf8);
lofwcv_1 = filtfilt(b,a,wcv_minus_spike);

%zscore
loflf8_1 = loflf8_1 - mean(loflf8_1);
loflf8_1 = loflf8_1/std(loflf8_1);
lofwcv_1 = lofwcv_1 - mean(lofwcv_1);
lofwcv_1 = lofwcv_1/std(lofwcv_1);

lf8_up_start_sh_pts = zeros(size(lf8_up_start_time));
lf8_up_end_sh_pts = lf8_up_start_sh_pts;
%for lf8
lf8_diff = [0;diff(loflf8_1)];
for i = 1:length(lf8_up_start_time)
    
   %find first point after start of up with negative slope and where signal
   %is above mean (z = 0)
   [tempx] = find((lf8_diff(lf8_up_start_points(i):lf8_up_end_points(i)) < 0) & (loflf8_1(lf8_up_start_points(i):lf8_up_end_points(i)) > 0),1,'first');
   if ~isempty(tempx)    
       lf8_up_start_sh_pts(i) = tempx+lf8_up_start_points(i);
%        lf8_up_start_sh_pts(i,2) = loflf8_1(lf8_up_start_sh_pts(i,1));
   else
       lf8_up_start_sh_pts(i) = NaN;
   end
   
   %find last point before end of up with positive slope
   [tempx] = find((lf8_diff(lf8_up_start_points(i):lf8_up_end_points(i)) > 0) & (loflf8_1(lf8_up_start_points(i):lf8_up_end_points(i)) > 0),1,'last');
   if ~isempty(tempx)
       lf8_up_end_sh_pts(i) = tempx+lf8_up_start_points(i);
%        lf8_up_end_sh_pts(i,2) = loflf8_1(lf8_up_end_sh_pts(i,1));
   else
       lf8_up_end_sh_pts(i) = NaN;
   end
   
%    %fit first order poly to up state
%    upIndLength = length(lf8_up_start_sh_pts(i):lf8_up_end_sh_pts(i));
%    p = polyfit([1:upIndLength]'*dt,loflf8_1(lf8_up_start_sh_pts(i):lf8_up_end_sh_pts(i)),1);
%    lf8_up_slope(i) = p(1);
%    lf8_up_offset(i) = p(2);
   
end

%for wcv
wcv_up_start_sh_pts = zeros(size(wcv_up_start_time));
wcv_up_end_sh_pts = wcv_up_start_sh_pts;
%for lf8
wcv_diff = [0;diff(lofwcv_1)];
for i = 1:length(wcv_up_start_time)
    
   %find first point after start of up with negative slope and where signal
   %is above mean (z = 0)
   [tempx] = find((wcv_diff(wcv_up_start_points(i):wcv_up_end_points(i)) < 0) & (lofwcv_1(wcv_up_start_points(i):wcv_up_end_points(i)) > 0),1,'first');
   if ~isempty(tempx)    
       wcv_up_start_sh_pts(i) = tempx+wcv_up_start_points(i);
   else
       wcv_up_start_sh_pts(i) = NaN;
   end
   
   %find last point before end of up with positive slope
   [tempx] = find((wcv_diff(wcv_up_start_points(i):wcv_up_end_points(i)) > 0) & (lofwcv_1(wcv_up_start_points(i):wcv_up_end_points(i)) > 0),1,'last');
   if ~isempty(tempx)
       wcv_up_end_sh_pts(i) = tempx+wcv_up_start_points(i);
   else
       wcv_up_end_sh_pts(i) = NaN;
   end
   
   %% fit first order poly to up state
%    upIndLength = length(wcv_up_start_sh_pts(i):wcv_up_end_sh_pts(i));
%    p = polyfit([1:upIndLength]'*dt,lofwcv_1(wcv_up_start_sh_pts(i):wcv_up_end_sh_pts(i)),1);
%    wcv_up_slope(i) = p(1);
%    wcv_up_offset(i) = p(2);

end

% %plot running average of slope for lfp and wcv
% 
% smooth_lf8_up_slope = fgsmooth(lf8_up_slope,2);
% smooth_wcv_up_slope = fgsmooth(wcv_up_slope,2);
% 
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(lf8_up_start_time,smooth_lf8_up_slope,'.-','linewidth',2)
% hold on
% plot(wcv_up_start_time,smooth_wcv_up_slope,'r.-','linewidth',2)
% legend('LF8','WCV')
% xlabel('Time (s)','FontSize',14)
% ylabel('Smoothed Slope (s)','FontSize',14)
% title('Smoothed Up State Slope vs Time')
% print('-dpng','./jmm/sm_up_slope_v_time')
% close all

%% visualize lfp up slope fit results

% segLength = 40 %segment length in seconds
% lnwdth = 2;
% 
% firstTime = synct(1);
% nsynct = synct-firstTime;
% times1 = synct(lf8_up_start_sh_pts(~isnan(lf8_up_start_sh_pts)));
% times2 = synct(lf8_up_end_sh_pts(~isnan(lf8_up_end_sh_pts)));
% times1 = times1 - firstTime;
% times2 = times2 - firstTime;
% sig1 = loflf8_1;
% pointsPerSeg = round(segLength/dt);
% 
% numSegs = floor(length(sig1)/pointsPerSeg);
% 
% upStateCnt = 1;
% for i = 1:numSegs
% 
%     begID = pointsPerSeg*(i-1)+1;
%     endID = pointsPerSeg*i;
% 
%     %zscore signals within segment window
% %     sig1(begID:endID) = sig1(begID:endID) - mean(sig1(begID:endID));
% %     sig1(begID:endID) = sig1(begID:endID)/std(sig1(begID:endID));
% % 
%     begTime = (i-1)*segLength;
%     endTime = i*segLength;
%     plot(nsynct(begID:endID),sig1(begID:endID),'linewidth',lnwdth)
%     hold on
%     plot(times1(times1>begTime&times1<endTime),sig1(round(times1(times1>begTime&times1<endTime)/dt)),'go','linewidth',4,'MarkerSize',4)
%     plot(times2(times2>begTime&times2<endTime),sig1(round(times2(times2>begTime&times2<endTime)/dt)),'ro','linewidth',4,'MarkerSize',4)
% 
% 
%     curUpStates = find(times1>begTime&times1<endTime);
%     for j = 1:length(curUpStates)
%     
%         upBegInd = round(times1(curUpStates(j))/dt);
%         upEndInd = round(times2(curUpStates(j))/dt);
%         upStateDur = upEndInd-upBegInd + 1;
%         plot(nsynct(upBegInd:upEndInd),polyval([lf8_up_slope(j-1+upStateCnt) lf8_up_offset(j-1+upStateCnt)],dt*[1:upStateDur]),'r','linewidth',2)
%         
%     end
%     upStateCnt = upStateCnt + length(curUpStates);
%     clear curUpStates
%     pause
%     clf
% 
% 
% end

%% visualize wcv up slope fit results

% segLength = 40 %segment length in seconds
% lnwdth = 2;
% 
% firstTime = synct(1);
% nsynct = synct-firstTime;
% times1 = synct(wcv_up_start_sh_pts(~isnan(wcv_up_start_sh_pts)));
% times2 = synct(wcv_up_end_sh_pts(~isnan(wcv_up_end_sh_pts)));
% times1 = times1 - firstTime;
% times2 = times2 - firstTime;
% sig1 = lofwcv_1;
% pointsPerSeg = round(segLength/dt);
% 
% numSegs = floor(length(sig1)/pointsPerSeg);
% 
% upStateCnt = 1;
% for i = 1:numSegs
% 
%     begID = pointsPerSeg*(i-1)+1;
%     endID = pointsPerSeg*i;
% 
%     %zscore signals within segment window
% %     sig1(begID:endID) = sig1(begID:endID) - mean(sig1(begID:endID));
% %     sig1(begID:endID) = sig1(begID:endID)/std(sig1(begID:endID));
% % 
%     begTime = (i-1)*segLength;
%     endTime = i*segLength;
%     plot(nsynct(begID:endID),sig1(begID:endID),'linewidth',lnwdth)
%     hold on
%     plot(times1(times1>begTime&times1<endTime),sig1(round(times1(times1>begTime&times1<endTime)/dt)),'go','linewidth',4,'MarkerSize',4)
%     plot(times2(times2>begTime&times2<endTime),sig1(round(times2(times2>begTime&times2<endTime)/dt)),'ro','linewidth',4,'MarkerSize',4)
% 
% 
%     curUpStates = find(times1>begTime&times1<endTime);
%     for j = 1:length(curUpStates)
%     
%         upBegInd = round(times1(curUpStates(j))/dt);
%         upEndInd = round(times2(curUpStates(j))/dt);
%         upStateDur = upEndInd-upBegInd + 1;
%         plot(nsynct(upBegInd:upEndInd),polyval([wcv_up_slope(j-1+upStateCnt) wcv_up_offset(j-1+upStateCnt)],dt*[1:upStateDur]),'r','linewidth',2)
%         
%     end
%     upStateCnt = upStateCnt + length(curUpStates);
%     clear curUpStates
%     pause
%     clf
% 
% 
% end
% 

%% Get variance during up state
% 
% lf8_up_state_var = zeros(size(lf8_up_start_sh_pts));
% wcv_up_state_var = zeros(size(wcv_up_start_sh_pts));
% 
% %for lf8
% for i = 1:length(lf8_up_start_sh_pts)
%    
%     dataSeg = lf8(lf8_up_start_sh_pts(i):lf8_up_end_sh_pts(i));
%     lf8_up_state_var(i) = var(dataSeg);
%     clear dataSeg
%     
% end
% 
% %for wcv
% for i = 1:length(wcv_up_start_sh_pts)
%     
%    dataSeg = wcv(wcv_up_start_sh_pts(i):wcv_up_end_sh_pts(i));
%    wcv_up_state_var(i) = var(dataSeg);
%    clear dataSeg
%     
% end
% 
% % figure
% % plot(lf8_up_start_time,lf8_up_state_var,'.-')
% % hold on
% % plot(wcv_up_start_time,wcv_up_state_var,'r.-')
% % legend('LF8','WCV')
% % xlabel('Time (s)','FontSize',14)
% % ylabel('Up State Variance','FontSize',14)
% % title('Up State Variance vs Time','FontSize',16)
% 
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% plot(lf8_up_start_time,fgsmooth(lf8_up_state_var,2),'.-')
% hold on
% plot(wcv_up_start_time,fgsmooth(wcv_up_state_var,2),'r.-')
% legend('LF8','WCV')
% xlabel('Time (s)','FontSize',14)
% ylabel('Up State Variance','FontSize',14)
% title('Up State Variance vs Time','FontSize',16)
% print('-dpng','./jmm/up_state_var_v_time')
% close all

%% spectrum (2-50) during up state
%for wcv
minupdur = 0.3; %min up state duration in s
NFFT = 2^12;
minFreq = 2;
maxFreq = 50;
minFreqID = round(minFreq/(2*nyqf/NFFT))+1;
maxFreqID = round(maxFreq/(2*nyqf/NFFT))+1;

wcv_up_power = zeros(length(wcv_up_end_sh_pts),maxFreqID-minFreqID+1);
wcv_trashStates = [];
for i = 1:length(wcv_up_start_sh_pts)
    
    upStateDur = (wcv_up_end_sh_pts(i) - wcv_up_start_sh_pts(i))*dt;
    if upStateDur > minupdur
        
        dataSeg = wcv_minus_spike(wcv_up_start_sh_pts(i):wcv_up_end_sh_pts(i));
        [tempPow,F] = pwelch(dataSeg,[],[],2^12,2*nyqf);
%         tempPow = tempPow/sum(tempPow);  %if you want relative power
        wcv_up_power(i,:) = tempPow(minFreqID:maxFreqID);
        
    else
       wcv_trashStates = [wcv_trashStates i]; 
    end
    clear dataSeg
end

F = F(minFreqID:maxFreqID);



lf8_up_power = zeros(length(lf8_up_end_sh_pts),maxFreqID-minFreqID+1);
%for lf8
lf8_trashStates = [];
for i = 1:length(lf8_up_start_sh_pts)
    
    upStateDur = (lf8_up_end_sh_pts(i) - lf8_up_start_sh_pts(i))*dt;
    if upStateDur > minupdur
        
        dataSeg = lf8(lf8_up_start_sh_pts(i):lf8_up_end_sh_pts(i));
        [tempPow,F] = pwelch(dataSeg,[],[],2^12,2*nyqf);
%         tempPow = tempPow/sum(tempPow); %if you want relative power
        lf8_up_power(i,:) = tempPow(minFreqID:maxFreqID);
        
    else
        lf8_trashStates = [lf8_trashStates i];
    end
    clear dataSeg
end

F = F(minFreqID:maxFreqID);

good_lf8_states = [1:length(lf8_up_start_sh_pts)];
good_lf8_states(lf8_trashStates) = [];
good_wcv_states = [1:length(wcv_up_start_sh_pts)];
good_wcv_states(wcv_trashStates) = [];

num_good_lf8_ups = length(good_lf8_states);
num_good_wcv_ups = length(good_wcv_states);

%if there are more lf8 up states find the closest lf8 up state to every wcv
%up state
if num_good_lf8_ups > num_good_wcv_ups
    
   corr_lf8_up = zeros(size(good_wcv_states));
    for i= 1:num_good_wcv_ups
        
        [a,corr_lf8_up(i)] = min(abs(lf8_up_start_time(good_lf8_states) - wcv_up_start_time(good_wcv_states(i))));
        
    end
    
else
    %if there are more wcv up states find the closest wcv up state to every
    %lf8 up state
    corr_wcv_up = zeros(size(good_lf8_states));
    for i= 1:num_good_lf8_ups
        
        [a,corr_wcv_up(i)] = min(abs(wcv_up_start_time(good_wcv_states) - lf8_up_start_time(good_lf8_states(i))));
        
    end
    
end


% close all
% figure
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
% subplot(3,1,1)
% pcolor(wcv_up_start_time(good_wcv_states),F,log(wcv_up_power(good_wcv_states,:)'));shading flat;colorbar
% set(gca,'Layer','Top')
% xlabel('Time (s)','FontSize',12)
% ylabel('Frequency (Hz)','FontSize',12)
% title('WCV Up State Spectra','FontSize',14)
% subplot(3,1,2)
% pcolor(lf8_up_start_time(good_lf8_states),F,log(lf8_up_power(good_lf8_states,:)'));shading flat;colorbar
% set(gca,'Layer','Top')
% xlabel('Time (s)','FontSize',12)
% ylabel('Frequency (Hz)','FontSize',12)
% title('LF8 Up State Spectra','FontSize',14)
% subplot(3,1,3)
% if num_good_lf8_ups > num_good_wcv_ups
%     pcolor(wcv_up_start_time(good_wcv_states),F,log(wcv_up_power(good_wcv_states,:)'./lf8_up_power(good_lf8_states(corr_lf8_up),:)'));shading flat;colorbar
% else
%     pcolor(lf8_up_start_time(good_lf8_states),F,log(wcv_up_power(good_wcv_states(corr_wcv_up),:)'./lf8_up_power(good_lf8_states,:)'));shading flat;colorbar
% end
% set(gca,'Layer','Top')
% xlabel('Time (s)','FontSize',12)
% ylabel('Frequency (Hz)','FontSize',12)
% title('Relative Power in WCV/LF8 Up States','FontSize',14)
% print('-dpng','./jmm/up_state_spectra')

close all
figure
set(gcf,'PaperUnits','centimeters');
set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
subplot(2,1,1)
pcolor(wcv_up_start_time(good_wcv_states),F,log(wcv_up_power(good_wcv_states,:)'));shading flat;colorbar
caxis([0 12]);colorbar
set(gca,'Layer','Top')
xlabel('Time (s)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)
title('WCV Up State Spectra','FontSize',14)
subplot(2,1,2)
pcolor(lf8_up_start_time(good_lf8_states),F,log(lf8_up_power(good_lf8_states,:)'));shading flat;colorbar
caxis([0 12]);colorbar
set(gca,'Layer','Top')
xlabel('Time (s)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)
title('LF8 Up State Spectra','FontSize',14)
print('-dpng','./jmm/up_state_spectra_norel_br')
