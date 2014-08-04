clear all

load C:\WC_Germany\JMM_Analysis_ste\dir_tree_ste
load C:\WC_Germany\JMM_Analysis_ste\run_hist_thresh\up_per_data_ste
g = [0:0.1:20];

segDur = 30;
winSlide = 1;
dsf = 8;
Fsd = 2016/dsf;


for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike
    if exist('spike_time.mat')
        load spike_time
    else
        load spike_time_br
    end
    
    sm_threshold_w{d} = jmm_smooth_1d(threshold_w{d},10);
    sm_threshold_8{d} = jmm_smooth_1d(threshold_8{d},10);
    
    disp('WCV')
[up_state_dur{d},down_state_dur{d},up_trans{d},down_trans{d},sig_w{d},int_thresh_w{d}] = UDS_extract_run_hist(wcv_minus_spike,sm_threshold_w{d});
disp('LFP')
[up_state_dur8{d},down_state_dur8{d},up_trans8{d},down_trans8{d},sig_8{d},int_thresh_8{d}] = UDS_extract_run_hist(lf8,sm_threshold_8{d});

% dataViewer_run_hist(sig_w{d},sig_8{d},int_thresh_w{d},int_thresh_8{d},up_trans{d},up_trans8{d},down_trans{d},down_trans8{d})

%calculate average rate per up state
up_spk_rate{d} = zeros(1,length(up_trans{d}));
for i = 1:length(up_trans{d})
   
    num_spikes = length(find(spkid>up_trans{d}(i)*dsf & spkid<down_trans{d}(i)*dsf));
    up_spk_rate{d}(i) = num_spikes/up_state_dur{d}(i);
    
end

%create state indicator vectors
wcv_state = zeros(1,round(length(lf8)/dsf));
for i = 1:length(up_trans{d})

    wcv_state(up_trans{d}(i):down_trans{d}(i)) = 1;
    
end

lf8_state = zeros(1,round(length(lf8)/dsf));
for i = 1:length(up_trans8{d})
   
    lf8_state(up_trans8{d}(i):down_trans8{d}(i)) = 1;
    
end


%count all spikes that occured when both lfp and wcv were in up state
%count all spikes that occured when only wcv was in up state
double_up_spikes = [];
single_up_spikes = [];
down_spikes = [];
for i = 1:length(spkid)
   
    first_up8 = find(up_trans8{d}*dsf < spkid(i),1,'last');
    first_up = find(up_trans{d}*dsf < spkid(i),1,'last');
    
    if ~isempty(first_up) & ~isempty(first_up8)
        if spkid(i) < down_trans{d}(first_up)*dsf
            if spkid(i) < down_trans8{d}(first_up8)*dsf
                double_up_spikes = [double_up_spikes i];
            else
                single_up_spikes = [single_up_spikes i];
            end
        elseif spkid(i) < down_trans{d}(end)*dsf
            down_spikes = [down_spikes i];
        end
    end
    
end

combined_up_time = sum(wcv_state == 1 & lf8_state == 1)/Fsd;
up_down8_time = sum(wcv_state == 1 & lf8_state == 0)/Fsd;

total_time = length(lf8)/2016;

up8_cond_spike(d) = length(double_up_spikes)/combined_up_time;
down8_cond_spike(d) = length(single_up_spikes)/up_down8_time;
up_rate(d) = (length(double_up_spikes)+length(single_up_spikes))/sum(wcv_state)*Fsd;

%calculate average up state duration in segmented data
nsamples = length(sig_w{d});
numWins = floor((nsamples - segDur*Fsd)/(winSlide*Fsd));
for w = 1:numWins
    
    begSamp = round(1+(w-1)*winSlide*Fsd);
    endSamp = begSamp + round(segDur*Fsd);
   
    %find current up states
    curUps_w = find(up_trans{d}>=begSamp&up_trans{d}<endSamp);
    curUps_8 = find(up_trans8{d}>=begSamp&up_trans8{d}<endSamp);

    curDowns_w = find(down_trans{d}>=begSamp&down_trans{d}<endSamp);
    curDowns_8 = find(down_trans8{d}>=begSamp&down_trans8{d}<endSamp);
    
    avg_up_dur{d}(w) = mean(up_state_dur{d}(curUps_w));
    avg_up_dur8{d}(w) = mean(up_state_dur8{d}(curUps_8));
    
    if ~isempty(curUps_w)
        max_up_dur{d}(w) = max(up_state_dur{d}(curUps_w));
    else
        max_up_dur{d}(w) = nan;
    end
    
    if ~isempty(curUps_8)
        max_up_dur8{d}(w) = max(up_state_dur8{d}(curUps_8));
    else
        max_up_dur8{d}(w) = nan;
    end
    curDowns_w(curDowns_w==length(down_trans{d})) = [];
    curDowns_8(curDowns_8==length(down_trans8{d})) = [];
    if ~isempty(curDowns_w)
        max_down_dur{d}(w) = max(down_state_dur{d}(curDowns_w));
    else
        max_down_dur{d}(w) = nan;
    end
    
    if ~isempty(curDowns_8)
        max_down_dur8{d}(w) = max(down_state_dur8{d}(curDowns_8));
    else
        max_down_dur8{d}(w) = nan;
    end
    
    
    avg_down_dur{d}(w) = mean(down_state_dur{d}(curDowns_w));
    avg_down_dur8{d}(w) = mean(down_state_dur8{d}(curDowns_8));
    
end
% 

end_time = length(max_up_dur{d});

% Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
% plot(max_up_dur{d},'linewidth',2)
% hold on
% plot(max_up_dur8{d},'r','linewidth',2)
% xlim([0 end_time])
% subplot(2,1,2)
% plot(up_trans{d}/Fsd,up_spk_rate{d},'o-')
% xlim([0 end_time])
% tname = ['C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\up_spike_rate_oversmooth_' f_names{d}];
%     print('-dpng',tname)




% %for each membrane up transition find nearest preceeding lfp up 
for t = 1:length(up_trans{d})
    
    prev_up = up_trans8{d}(find(up_trans8{d}<up_trans{d}(t),1,'last'));
    if ~isempty(prev_up)
        near_lf8_up{d}(t) = prev_up;
    else
        near_lf8_up{d}(t) = nan;
    end
    
    
end
% 
% % up_dist{d} = hist(up_state_dur{d},g);
% % up_dist{d} = up_dist{d}/sum(up_dist{d});
% % down_dist{d} = hist(down_state_dur{d},g);
% % down_dist{d} = down_dist{d}/sum(down_dist{d});
% % 
% % up_dist8{d} = hist(up_state_dur8{d},g);
% % up_dist8{d} = up_dist8{d}/sum(up_dist8{d});
% % down_dist8{d} = hist(down_state_dur8{d},g);
% % down_dist8{d} = down_dist8{d}/sum(down_dist8{d});
%   
[up_dist{d},up_r] = log_hist(up_state_dur{d},[0.3 20],100);
[down_dist{d},down_r] = log_hist(down_state_dur{d},[0.3 20],100);
[up_dist8{d},up_r] = log_hist(up_state_dur8{d},[0.3 20],100);
[down_dist8{d},down_r]  = log_hist(down_state_dur8{d},[0.3 20],100);
% 
% Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
% plot(up_state_dur{d},(up_trans{d}-near_lf8_up{d})/Fsd,'.')
% xlabel('Up State Duration','FontSize',14)
% ylabel('Up Transition Lag','FontSize',14)
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% [a,b] = corrcoef(up_state_dur{d}(~isnan(near_lf8_up{d})),(up_trans{d}(~isnan(near_lf8_up{d}))-near_lf8_up{d}(~isnan(near_lf8_up{d}))));
% title(sprintf('Corr = %d   P = %d',a(2,1),b(2,1)))
% tname = ['C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\dur_lag_corr' f_names{d}];
%     print('-dpng',tname)
% 
Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
subplot(3,1,1)
plot(up_trans{d}/Fsd,(up_trans{d}-near_lf8_up{d})/Fsd,'o')
xlim([0 length(avg_up_dur{d})])
set(gca,'yscale','log')
xlabel('Time (s)','FontSize',14)
ylabel('Up Transition Lag','FontSize',14)
subplot(3,1,2)
plot(max_up_dur{d},'linewidth',2)
hold on
plot(max_up_dur8{d},'r','linewidth',2)
my = max(max_up_dur{d});
xlabel('Time','FontSize',14)
ylabel('Max Up State Duration','FontSize',14)
xlim([0 length(max_up_dur{d})])
subplot(3,1,3)
plot(avg_up_dur{d},'linewidth',2)
hold on
plot(avg_up_dur8{d},'r','linewidth',2)
xlim([0 length(max_up_dur{d})])

xlabel('Time','FontSize',14)
ylabel('Avg Up State Duration','FontSize',14)

    tname = ['C:\WC_Germany\JMM_Analysis_ste\UDS_dur_run_hist\dur_lag_10_28' f_names{d}];
    print('-dpng',tname)
%     
%     
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
%     max_up = max([prctile(up_state_dur{d},90) prctile(up_state_dur8{d},90)]);
   stairs(up_r,up_dist{d})
    hold on
    stairs(up_r,up_dist8{d},'r')
    xlim([0.2 20])
    set(gca,'xscale','log')
    legend('WCV','LF8')
    title('Up State Duration','FontSize',14)
    grid on
    subplot(2,1,2)
    stairs(down_r,down_dist{d})
    hold on
    stairs(down_r,down_dist8{d},'r')
%     max_down = max([prctile(down_state_dur{d},90) prctile(down_state_dur8{d},90)]);
    xlim([0.2 20])
        set(gca,'xscale','log')
        grid on
    legend('WCV','LF8')
    title('Down State Duration','FontSize',14)
    grid on
    tname = ['C:\WC_Germany\JMM_Analysis_ste\UDS_dur_run_hist\' f_names{d}];
    print('-dpng',tname)
    close

    save C:\WC_Germany\JMM_Analysis_ste\UDS_dur_run_hist\data avg* max* up_state_dur* down_state_dur* up_trans* down_trans* g up_dist* down_dist* 
%     
    clear lf8 wcv_minus_spike spkid
   d 
end

for d = 1:16
    up_dist_s(d,:) = up_dist{d};
    up_dist_s8(d,:) = up_dist8{d};
    down_dist_s(d,:) = down_dist{d};
    down_dist_s8(d,:) = down_dist8{d};
end

m_u_s = mean(up_dist_s);
u_u_s = m_u_s+2*std(up_dist_s)/sqrt(16);
l_u_s = m_u_s-2*std(up_dist_s)/sqrt(16);

m_d_s = mean(down_dist_s);
u_d_s = m_d_s + 2*std(down_dist_s)/sqrt(16);
l_d_s = m_d_s - 2*std(down_dist_s)/sqrt(16);

m_u_s8 = mean(up_dist_s8);
u_u_s8 = m_u_s8+2*std(up_dist_s8)/sqrt(16);
l_u_s8 = m_u_s8-2*std(up_dist_s8)/sqrt(16);

m_d_s8 = mean(down_dist_s8);
u_d_s8 = m_d_s8 + 2*std(down_dist_s8)/sqrt(16);
l_d_s8 = m_d_s8 - 2*std(down_dist_s8)/sqrt(16);


stairs(up_r,m_u_s,'linewidth',2)
hold on
stairs(up_r,m_u_s8,'r','linewidth',2)
legend('WCV','LFP')
stairs(up_r,u_u_s,'--')
stairs(up_r,l_u_s,'--')
stairs(up_r,u_u_s8,'r--')
stairs(up_r,l_u_s8,'r--')
set(gca,'xscale','log')


stairs(down_r,m_d_s,'linewidth',2)
hold on
stairs(down_r,m_d_s8,'r','linewidth',2)
legend('WCV','LFP')
stairs(down_r,u_d_s,'--')
stairs(down_r,l_d_s,'--')
stairs(down_r,u_d_s8,'r--')
stairs(down_r,l_d_s8,'r--')
set(gca,'xscale','log')
grid