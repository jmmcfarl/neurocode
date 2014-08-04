clear all

load C:\WC_Germany\Layer5\layer_5_dir
load C:\WC_Germany\Layer5\UDS_dur_run_hist\data
Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

lags = [-1*Fsd:1*Fsd];
maxlag = 1*Fsd;

for d = [1:6 8:length(dir_array)]
    
   cd(dir_array{d});
   pwd
   
   if ~exist('spike_time.mat')
   load spike_time_br
   else
       load spike_time
   end
   load used_data wcv_minus_spike lf8
   
   wcv = filtfilt(b,a,wcv_minus_spike);
   wcv = downsample(wcv,dsf);
   wcv = zscore(wcv);
   wcv = wcv';
   
   lf8 = filtfilt(b,a,lf8);
   lf8 = downsample(lf8,dsf);
   lf8 = zscore(lf8);
   lf8 = lf8';
      
   %find lfp up transitions that occur during a wcv up state (wcv must have
   %already been in up state for at least 0.5 s)
   good_ups = [];
   good_wcv_ups = [];
   for i = 1:length(up_trans8{d})
      
       %find previous wcv up transition
       prev_up = find(up_trans{d} < up_trans8{d}(i) - Fsd,1,'last');
       
       if up_trans8{d}(i) < down_trans{d}(prev_up)
           good_ups = [good_ups i];
           good_wcv_ups = [good_wcv_ups prev_up];
       end
       
   end
   
    %find lfp down transitions that occur during wcv up state
    good_downs = [];
    for i = 1:length(down_trans8{d})
       
        prev_up = find(up_trans{d} < down_trans8{d}(i) - 0.3*Fsd,1,'last');
        
        if down_trans8{d}(i) <= down_trans{d}(prev_up)
            good_downs = [good_downs i];
        end
        
    end
  
    
   %find lfp up transitions that occur during a wcv down state
   lup_wdown = [];
%    wcv_lup_wdown = [];
    next_wcv_delay{d} = [];
   for i = 1:length(up_trans8{d})
      
       %find previous wcv down transition
       prev_down = find(down_trans{d} < up_trans8{d}(i),1,'last');
       
       if prev_down < length(up_trans{d})
           if up_trans8{d}(i) < up_trans{d}(prev_down+1)
               lup_wdown = [lup_wdown i];
%                wcv_lup_wdown = [wcv_lup_wdown prev_down];
                next_wcv_delay{d} = [next_wcv_delay{d} (up_trans{d}(prev_down+1)-up_trans8{d}(i))];
           end
       end
       
   end

   %find lfp down transitions that occur during a wcv down state
   ldown_wdown = [];
   for i = 1:length(down_trans8{d})
      
       %find previous wcv up transition
       prev_down = find(down_trans{d} < down_trans8{d}(i),1,'last');
       
       if prev_down < length(up_trans{d})
           if down_trans8{d}(i) < up_trans{d}(prev_down+1)
               ldown_wdown = [ldown_wdown i];
           end
       end
       
   end
 
   
   
    cor_up_trans = up_trans8{d}(good_ups);        
    cor_lup_wdown_trans = up_trans8{d}(lup_wdown);
    cor_down_trans = down_trans8{d}(good_downs);
    cor_ldown_wdown = down_trans8{d}(ldown_wdown);
    
    %for wcv down conditional lfp ups
    if length(cor_lup_wdown_trans) > 0
        wcv_lup_wdown = zeros(length(cor_lup_wdown_trans),length(lags));     
        lf8_lup_wdown = wcv_lup_wdown;
        for i = 1:length(cor_lup_wdown_trans)
            if cor_lup_wdown_trans(i) > maxlag & cor_lup_wdown_trans(i) < length(wcv) - maxlag
                wcv_lup_wdown(i,:) = wcv(cor_lup_wdown_trans(i)-maxlag:cor_lup_wdown_trans(i)+maxlag);
                lf8_lup_wdown(i,:) = lf8(cor_lup_wdown_trans(i)-maxlag:cor_lup_wdown_trans(i)+maxlag);
            else
                wcv_lup_wdown(i,:) = nan;
                lf8_lup_wdown(i,:) = nan;
            end
        end
    end

 %for wup_lup
    wcv_lup_wup = zeros(length(cor_up_trans),length(lags));
    lf8_lup_wup = wcv_lup_wup;
    for i = 1:length(cor_up_trans)     
       
        if cor_up_trans(i) > maxlag & cor_up_trans(i) < length(wcv)-maxlag

        wcv_lup_wup(i,:) = wcv(cor_up_trans(i)-maxlag:cor_up_trans(i)+maxlag);
        lf8_lup_wup(i,:) = lf8(cor_up_trans(i)-maxlag:cor_up_trans(i)+maxlag);
        
        else
            wcv_lup_wup(i,:) = nan;
            lf8_lup_wup(i,:) = nan;
        end

    end
 
    %for wup_ldown
        wcv_ldown_wup = zeros(length(cor_down_trans),length(lags));
        lf8_ldown_wup = wcv_ldown_wup;
    for i = 1:length(cor_down_trans)     
       
        if cor_down_trans(i) > maxlag & cor_down_trans(i) < length(wcv)-maxlag

        wcv_ldown_wup(i,:) = wcv(cor_down_trans(i)-maxlag:cor_down_trans(i)+maxlag);
        lf8_ldown_wup(i,:) = lf8(cor_down_trans(i)-maxlag:cor_down_trans(i)+maxlag);
        else
            wcv_ldown_wup(i,:) = nan;
            lf8_ldown_wup(i,:) = nan;
        end

    end

    %for wdown_ldown
        wcv_ldown_wdown = zeros(length(cor_ldown_wdown),length(lags));
        lf8_ldown_wdown = wcv_ldown_wdown;
    for i = 1:length(cor_ldown_wdown)     
       
        if cor_ldown_wdown(i) > maxlag & cor_ldown_wdown(i) < length(wcv)-maxlag

        wcv_ldown_wdown(i,:) = wcv(cor_ldown_wdown(i)-maxlag:cor_ldown_wdown(i)+maxlag);
        lf8_ldown_wdown(i,:) = lf8(cor_ldown_wdown(i)-maxlag:cor_ldown_wdown(i)+maxlag);
        else
            wcv_ldown_wdown(i,:) = nan;
            lf8_ldown_wdown(i,:) = nan;
        end

    end
    
    
    m_wcv_lup_wup(d,:) = nanmean(wcv_lup_wup);
    m_wcv_lup_wdown(d,:) = nanmean(wcv_lup_wdown);
    m_lf8_lup_wup(d,:) = nanmean(lf8_lup_wup);
    m_lf8_lup_wdown(d,:) = nanmean(lf8_lup_wdown);
    m_wcv_ldown_wup(d,:) = nanmean(wcv_ldown_wup);
    m_wcv_ldown_wdown(d,:) = nanmean(wcv_ldown_wdown);
    m_lf8_ldown_wup(d,:) = nanmean(lf8_ldown_wup);
    m_lf8_ldown_wdown(d,:) = nanmean(lf8_ldown_wdown);
%% calculate up transition delays

%estimate density of up delays
[y,x] = gpkde(next_wcv_delay{d}'/Fsd);

[dummy,most_likely_delay(d)] = max(y);
most_likely_delay(d) = x(most_likely_delay(d));

max_wcv_val = max(m_wcv_lup_wdown(d,:));
max_lf8_val = max(m_lf8_lup_wdown(d,:));
min_wcv_val = min(m_wcv_lup_wdown(d,:));
min_lf8_val = min(m_lf8_lup_wdown(d,:));

lf8_tpoint = find(m_lf8_lup_wdown(d,:) > (max_lf8_val+min_lf8_val)/2,1,'first');
wcv_tpoint = lf8_tpoint+find(m_wcv_lup_wdown(d,lf8_tpoint:end) > (max_wcv_val+min_wcv_val)/2,1,'first');

wcv_up_delay(d) = (wcv_tpoint - lf8_tpoint)/Fsd;


plot(lags/Fsd,m_wcv_lup_wdown(d,:))
hold on
plot(lags/Fsd,m_lf8_lup_wdown(d,:),'r')
plot(lags(wcv_tpoint)/Fsd,m_wcv_lup_wdown(d,wcv_tpoint),'o')
plot(lags(lf8_tpoint)/Fsd,m_lf8_lup_wdown(d,lf8_tpoint),'ro')
title(sprintf('midpoint: %0.2g  mode: %0.2g',wcv_up_delay(d),most_likely_delay(d)))
t_names = ['C:\WC_Germany\Layer5\State_Trig\delay\' f_names{d}];
print('-dpng',t_names)
close all
    
    
    %smooth average
%     m_wcv_lup_wup = jmm_smooth_1d(m_wcv_lup_wup,5);
%     m_wcv_lup_wdown = jmm_smooth_1d(m_wcv_lup_wdown,5);
%     m_wcv_ldown_wup = jmm_smooth_1d(m_wcv_ldown_wup,5);
%     m_wcv_ldown_wdown = jmm_smooth_1d(m_wcv_ldown_wdown,5);
    
%% use jitter method to estimate significance of wcv modulation
%     max_jit = [0.1 0.2 0.3 0.4 0.5];
%     minBins = 20;
%     sig_t = 0;
%     
%     first_check = find(lags > -0.5*Fsd,1,'first');
%     last_check = find(lags > 0.5*Fsd,1,'first');
 
    %for wcvuplup
%     for j = 1:length(max_jit)
%         
%         [jit_avg,jit_ci] = jitter_average(wcv,cor_up_trans,round(max_jit(j)*Fsd),maxlag);
%         %check if average is outside ci 
% %         above_ci = find(m_wcv_lup_wup(first_check:last_check) > jit_ci(2,first_check:last_check));
%         below_ci = find(m_wcv_lup_wup(first_check:last_check) < jit_ci(1,first_check:last_check));
%         if  length(below_ci) > minBins
%             sig_t = max_jit(j);
%             break
%         else
%             clear jit_avg jit_ci above_ci below_ci
%         end
%     end
%     
%     if sig_t > 0
%         plot(lags/Fsd,jit_ci(1,:),'b')
%         hold on
%         plot(lags/Fsd,jit_ci(2,:),'r')
%         plot(lags/Fsd,m_wcv_lup_wup,'k','linewidth',2)
%         title(sprintf('significance at %0.2g',sig_t))
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wup_lup\sig_test_' f_names{d}];
%         print('-dpng',t_names)
%         close all
%     else
%         disp('no significant modulation')
%     end
    
   
    
%   for j = 1:length(max_jit)
%         
%         [jit_avg,jit_ci] = jitter_average(wcv,cor_ldown_wdown,round(max_jit(j)*Fsd),maxlag);
%         %check if average is outside ci 
%         above_ci = find(m_wcv_ldown_wdown(first_check:last_check) > jit_ci(2,first_check:last_check));
%         below_ci = find(m_wcv_ldown_wdown(first_check:last_check) < jit_ci(1,first_check:last_check));
%         if length(above_ci) > minBins | length(below_ci) > minBins
%             sig_t = max_jit(j);
%             break
%         else
%             clear jit_avg jit_ci above_ci below_ci
%         end
%     end
%     
%     if sig_t > 0
%         plot(lags/Fsd,jit_ci(1,:),'b')
%         hold on
%         plot(lags/Fsd,jit_ci(2,:),'r')
%         plot(lags/Fsd,m_wcv_ldown_wdown,'k','linewidth',2)
%         title(sprintf('significance at %0.2g',sig_t))
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wdown_ldown\sig_test_' f_names{d}];
%         print('-dpng',t_names)
%         close all
%     else
%         disp('no significant modulation')
%     end

    
end

%% analyze average responses

% mean_w_lup_wup = mean(wcv_lup_wup);
% u_w_lup_wup = mean_w_lup_wup+2*std(wcv_lup_wup)/sqrt(17);
% l_w_lup_wup = mean_w_lup_wup-2*std(wcv_lup_wup)/sqrt(17);
% 
% mean_l_lup_wup = mean(lf8_lup_wup);
% u_l_lup_wup = mean_l_lup_wup+2*std(lf8_lup_wup)/sqrt(17);
% l_l_lup_wup = mean_l_lup_wup-2*std(lf8_lup_wup)/sqrt(17);
% 
% plot(lags/Fsd,mean_w_lup_wup,'linewidth',2)
% hold on
% plot(lags/Fsd,mean_l_lup_wup,'r','linewidth',2)
% legend('WCV','LFP')
% plot(lags/Fsd,u_w_lup_wup,'--')
% plot(lags/Fsd,l_w_lup_wup,'--')
% plot(lags/Fsd,u_l_lup_wup,'r--')
% plot(lags/Fsd,l_l_lup_wup,'r--')
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)
% 
% 
% mean_w_lup_wdown = mean(wcv_lup_wdown);
% u_w_lup_wdown = mean_w_lup_wdown+2*std(wcv_lup_wdown)/sqrt(17);
% l_w_lup_wdown = mean_w_lup_wdown-2*std(wcv_lup_wdown)/sqrt(17);
% 
% mean_l_lup_wdown = mean(lf8_lup_wdown);
% u_l_lup_wdown = mean_l_lup_wdown+2*std(lf8_lup_wdown)/sqrt(17);
% l_l_lup_wdown = mean_l_lup_wdown-2*std(lf8_lup_wdown)/sqrt(17);
% 
% plot(lags/Fsd,mean_w_lup_wdown,'linewidth',2)
% hold on
% plot(lags/Fsd,mean_l_lup_wdown,'r','linewidth',2)
% legend('WCV','LFP')
% plot(lags/Fsd,u_w_lup_wdown,'--')
% plot(lags/Fsd,l_w_lup_wdown,'--')
% plot(lags/Fsd,u_l_lup_wdown,'r--')
% plot(lags/Fsd,l_l_lup_wdown,'r--')
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)

% mean_w_ldown_wup = mean(wcv_ldown_wup);
% u_w_ldown_wup = mean_w_ldown_wup+2*std(wcv_ldown_wup)/sqrt(17);
% l_w_ldown_wup = mean_w_ldown_wup-2*std(wcv_ldown_wup)/sqrt(17);
% 
% mean_l_ldown_wup = mean(lf8_ldown_wup);
% u_l_ldown_wup = mean_l_ldown_wup+2*std(lf8_ldown_wup)/sqrt(17);
% l_l_ldown_wup = mean_l_ldown_wup-2*std(lf8_ldown_wup)/sqrt(17);
% 
% plot(lags/Fsd,mean_w_ldown_wup,'linewidth',2)
% hold on
% plot(lags/Fsd,mean_l_ldown_wup,'r','linewidth',2)
% legend('WCV','LFP')
% plot(lags/Fsd,u_w_ldown_wup,'--')
% plot(lags/Fsd,l_w_ldown_wup,'--')
% plot(lags/Fsd,u_l_ldown_wup,'r--')
% plot(lags/Fsd,l_l_ldown_wup,'r--')
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)


% mean_w_ldown_wdown = mean(wcv_ldown_wdown);
% u_w_ldown_wdown = mean_w_ldown_wdown+2*std(wcv_ldown_wdown)/sqrt(17);
% l_w_ldown_wdown = mean_w_ldown_wdown-2*std(wcv_ldown_wdown)/sqrt(17);
% 
% mean_l_ldown_wdown = mean(lf8_ldown_wdown);
% u_l_ldown_wdown = mean_l_ldown_wdown+2*std(lf8_ldown_wdown)/sqrt(17);
% l_l_ldown_wdown = mean_l_ldown_wdown-2*std(lf8_ldown_wdown)/sqrt(17);
% 
% plot(lags/Fsd,mean_w_ldown_wdown,'linewidth',2)
% hold on
% plot(lags/Fsd,mean_l_ldown_wdown,'r','linewidth',2)
% legend('WCV','LFP')
% plot(lags/Fsd,u_w_ldown_wdown,'--')
% plot(lags/Fsd,l_w_ldown_wdown,'--')
% plot(lags/Fsd,u_l_ldown_wdown,'r--')
% plot(lags/Fsd,l_l_ldown_wdown,'r--')
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)