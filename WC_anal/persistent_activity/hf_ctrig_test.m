clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
% load C:\WC_Germany\Persistent_activity\sigmoid_fit\sig_fit_data
load C:\WC_Germany\Persistent_activity\lf8_period_f_data2
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 4;
Fsd = 2016/dsf;
backlag = 5*Fsd;
forwardlag = 10*Fsd
lags = -backlag:forwardlag;
% dlags = -backlag:backlag;

%hf parameters
sm_length = round(0.4*Fsd);
sm_win = ones(1,sm_length)/sm_length;

lcf = 40/niqf;
hcf = 100/niqf;
[b2,a2] = butter(2,[lcf hcf]);

lcf = 60/niqf;
hcf = 250/niqf;
[b3,a3] = butter(2,[lcf hcf]);



% for d = 1:length(dir_array)
d = 14
    cd(dir_array{d})
    pwd
up_trans{d} = up_trans{d}*2;
down_trans{d} = down_trans{d}*2;
up_trans8{d} = up_trans8{d}*2;
down_trans8{d} = down_trans8{d}*2;

    load used_data lf8 wcv_minus_spike
    
    if ~exist('spike_time.mat')
        load spike_time_br
    else
        load spike_time
    end

    spike_id = round(spkid/dsf);
    
    
wcv = wcv_minus_spike;
    lf8_f = filtfilt(b,a,lf8);
    wcv_f = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8_f,dsf);
    wcv_d = downsample(wcv_f,dsf);

    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);
    
    spike_times = zeros(size(lf8_d));
    spike_times(spike_id) = 1;
    
    spk_rate = jmm_smooth_1d(spike_times,sm_length);
    spk_rate = zscore(spk_rate);
    
wcv_f_med = filtfilt(b2,a2,wcv);
wcv_med = downsample(wcv_f_med,dsf);
wcv_med = sqrt(wcv_med.^2);
wcv_med = conv(wcv_med,sm_win);
wcv_med(1:floor(sm_length/2)) = [];
wcv_med(end-floor(sm_length/2)+1:end) = [];
wcv_med = zscore(wcv_med);

% wcv_f_high = filtfilt(b3,a3,wcv);
% wcv_high = downsample(wcv_f_high,dsf);
% wcv_high = sqrt(wcv_high.^2);
% wcv_high = conv(wcv_high,sm_win);
% wcv_high(1:floor(sm_length/2)) = [];
% wcv_high(end-floor(sm_length/2)+1:end) = [];
% wcv_high = zscore(wcv_high);

    wcv_up_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    lf8_up_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    med_up_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    spk_up_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    high_up_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    wcv_down_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    lf8_down_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    spk_down_ctrig_mat = zeros(length(synch_ups{d}),length(lags));

    max_sec_lag = 2*Fsd;
    sec_lag = -max_sec_lag:max_sec_lag;
    
    %initialize
    mp_down_state = zeros(length(synch_ups{d}),1);
    bad_rows = [];
    cnt = 1;
    for i = 1:length(synch_ups{d})

                %find first lfp down state after MP up transition
        first_in_lfp_down = find(down_trans8{d} > up_trans{d}(synch_ups{d}(i)),1,'first');
        %find nearest lfp up transition
        [dummy,near_lfp_up] = min(abs(up_trans8{d} - up_trans{d}(synch_ups{d}(i))));
        
        if ~isempty(first_in_lfp_down)
            %find all lfp up transitions that occur before mp comes back down
            secondary_lfp_states{i} = find(down_trans8{d} < down_trans{d}(synch_ups{d}(i)) ...
                & up_trans8{d} > down_trans8{d}(first_in_lfp_down));
            for su = 1:length(secondary_lfp_states{i})
                if up_trans8{d}(secondary_lfp_states{i}(su)) < length(lf8_d)-max_sec_lag
               mp_med_sec_up_trig(cnt,:) = wcv_med(up_trans8{d}(secondary_lfp_states{i}(su))-max_sec_lag:...
                   up_trans8{d}(secondary_lfp_states{i}(su))+max_sec_lag);
               mp_sec_up_trig(cnt,:) = wcv_d(up_trans8{d}(secondary_lfp_states{i}(su))-max_sec_lag:...
                   up_trans8{d}(secondary_lfp_states{i}(su))+max_sec_lag);
               lf8_sec_up_trig(cnt,:) = lf8_d(up_trans8{d}(secondary_lfp_states{i}(su))-max_sec_lag:...
                   up_trans8{d}(secondary_lfp_states{i}(su))+max_sec_lag);
               spk_sec_up_trig(cnt,:) = spk_rate(up_trans8{d}(secondary_lfp_states{i}(su))-max_sec_lag:...
                   up_trans8{d}(secondary_lfp_states{i}(su))+max_sec_lag);
               
               cnt = cnt+1;
                end
            end
            
            
        end
        
        if up_trans8{d}(near_lfp_up) > max_sec_lag & up_trans8{d}(near_lfp_up) < length(lf8_d) - max_sec_lag
            mp_med_first_up_trig(i,:) = wcv_med(up_trans8{d}(near_lfp_up)-max_sec_lag:...
                up_trans8{d}(near_lfp_up)+max_sec_lag);
            mp_first_up_trig(i,:) = wcv_d(up_trans8{d}(near_lfp_up)-max_sec_lag:...
                up_trans8{d}(near_lfp_up)+max_sec_lag);
            lf8_first_up_trig(i,:) = lf8_d(up_trans8{d}(near_lfp_up)-max_sec_lag:...
                up_trans8{d}(near_lfp_up)+max_sec_lag);
            spk_first_up_trig(i,:) = spk_rate(up_trans8{d}(near_lfp_up)-max_sec_lag:...
                up_trans8{d}(near_lfp_up)+max_sec_lag);

        else
            mp_med_first_up_trig(i,:) = nan(1,length(sec_lag));
            mp_first_up_trig(i,:) = nan(1,length(sec_lag));
            lf8_first_up_trig(i,:) = nan(1,length(sec_lag));
            spk_first_up_trig(i,:) = nan(1,length(sec_lag));
        end
        
        %calculate mup triggered averages
        if up_trans{d}(synch_ups{d}(i)) > backlag & up_trans{d}(synch_ups{d}(i)) ...
                < length(wcv_d)-forwardlag
            wcv_up_ctrig_mat(i,:) = wcv_d(up_trans{d}(synch_ups{d}(i))-backlag: ...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
            lf8_up_ctrig_mat(i,:) = lf8_d(up_trans{d}(synch_ups{d}(i))-backlag: ...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
            med_up_ctrig_mat(i,:) = wcv_med(up_trans{d}(synch_ups{d}(i))-backlag: ...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
%             high_up_ctrig_mat(i,:) = wcv_high(up_trans{d}(synch_ups{d}(i))-backlag: ...
%                 up_trans{d}(synch_ups{d}(i))+forwardlag);

        else
            wcv_up_ctrig_mat(i,:) = nan;
            lf8_up_ctrig_mat(i,:) = nan;
            med_up_ctrig_mat(i,:) = nan;
%             high_up_ctrig_mat(i,:) = nan;

            bad_rows = [bad_rows i];
        end

                %calculate mdown triggered averages
        if down_trans{d}(synch_ups{d}(i)) > backlag & down_trans{d}(synch_ups{d}(i)) ...
                < length(wcv_d)-forwardlag
            wcv_down_ctrig_mat(i,:) = wcv_d(down_trans{d}(synch_ups{d}(i))-backlag: ...
                down_trans{d}(synch_ups{d}(i))+forwardlag);
            lf8_down_ctrig_mat(i,:) = lf8_d(down_trans{d}(synch_ups{d}(i))-backlag: ...
                down_trans{d}(synch_ups{d}(i))+forwardlag);
        else
            wcv_down_ctrig_mat(i,:) = nan;
            lf8_down_ctrig_mat(i,:) = nan;
        end
        
    end


%     m_wcv_up_ctrig_mat(d,:) = nanmean(wcv_up_ctrig_mat);
%     u_wcv_up_ctrig_mat(d,:) = m_wcv_up_ctrig_mat(d,:)+2*nanstd(wcv_up_ctrig_mat)/sqrt(length(synch_ups{d}));
%         l_wcv_up_ctrig_mat(d,:) = m_wcv_up_ctrig_mat(d,:)-2*nanstd(wcv_up_ctrig_mat)/sqrt(length(synch_ups{d}));
% 
%     m_lf8_up_ctrig_mat(d,:) = nanmean(lf8_up_ctrig_mat);
%         u_lf8_up_ctrig_mat(d,:) = m_lf8_up_ctrig_mat(d,:)+2*nanstd(lf8_up_ctrig_mat)/sqrt(length(synch_ups{d}));
%         l_lf8_up_ctrig_mat(d,:) = m_lf8_up_ctrig_mat(d,:)-2*nanstd(lf8_up_ctrig_mat)/sqrt(length(synch_ups{d}));
% 
%     m_wcv_down_ctrig_mat(d,:) = nanmean(wcv_down_ctrig_mat);
%         u_wcv_down_ctrig_mat(d,:) = m_wcv_down_ctrig_mat(d,:)+2*nanstd(wcv_down_ctrig_mat)/sqrt(length(synch_ups{d}));
%         l_wcv_down_ctrig_mat(d,:) = m_wcv_down_ctrig_mat(d,:)-2*nanstd(wcv_down_ctrig_mat)/sqrt(length(synch_ups{d}));
% 
%     m_lf8_down_ctrig_mat(d,:) = nanmean(lf8_down_ctrig_mat);
%             u_lf8_down_ctrig_mat(d,:) = m_lf8_down_ctrig_mat(d,:)+2*nanstd(lf8_down_ctrig_mat)/sqrt(length(synch_ups{d}));
%         l_lf8_down_ctrig_mat(d,:) = m_lf8_down_ctrig_mat(d,:)-2*nanstd(lf8_down_ctrig_mat)/sqrt(length(synch_ups{d}));

    [dummy,up_order] = sort(up_state_dur{d}(synch_ups{d}));
med_up_ctrig_mat(:,1:2520) = nan;
for i = 1:length(synch_ups{d})
    down_ind = find(lags/Fsd > up_state_dur{d}(synch_ups{d}(i)),1,'first');
    med_up_ctrig_mat(i,down_ind:end) = nan;
end
    
    
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    
%     subplot(3,1,1)
%     pcolor(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,wcv_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))...
%         /length(synch_ups{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     xlim([-2 10])
        subplot(2,1,1)
    pcolor(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d})...
        ,lf8_up_ctrig_mat(up_order,:));shading flat;
    colorbar
    hold on
    plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
    line([0 0],[0 1],'Color','k')
    xlim([-2 10])

    subplot(2,1,2)
    pcolor(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d})...
        ,med_up_ctrig_mat(up_order,:));shading flat;
    colorbar
    hold on
    plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
    line([0 0],[0 1],'Color','k')
    xlim([-2 10])

m_mp_first(d,:) = nanmean(mp_first_up_trig);
m_lf8_first(d,:) = nanmean(lf8_first_up_trig);
m_med_first(d,:) = nanmean(mp_med_first_up_trig);
m_spk_first(d,:) = nanmean(spk_first_up_trig);

m_mp_sec(d,:) = nanmean(mp_sec_up_trig);
m_lf8_sec(d,:) = nanmean(lf8_sec_up_trig);
m_med_sec(d,:) = nanmean(mp_med_sec_up_trig);
m_spk_sec(d,:) = nanmean(spk_sec_up_trig);


% plot_with_ci(sec_lag/Fsd,mp_first_up_trig,'b')
% hold on
% plot_with_ci(sec_lag/Fsd,lf8_first_up_trig,'r')
% plot_with_ci(sec_lag/Fsd,mp_med_first_up_trig,'k')
% plot_with_ci(sec_lag/Fsd,spk_first_up_trig,'g')
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (z)','FontSize',14)
% title('MP Down / LFP UP trig Avg','FontSize',16)
%     tname = ['C:\WC_Germany\Persistent_activity\high_freq_lup\mp_down_' f_names{d}];
%     print('-dpng',tname);
%     close
% 
% plot_with_ci(sec_lag/Fsd,mp_sec_up_trig,'b')
% hold on
% plot_with_ci(sec_lag/Fsd,lf8_sec_up_trig,'r')
% plot_with_ci(sec_lag/Fsd,mp_med_sec_up_trig,'k')
% plot_with_ci(sec_lag/Fsd,spk_sec_up_trig,'g')
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (z)','FontSize',14)
% title('MP UP / LFP UP trig Avg','FontSize',16)
%     tname = ['C:\WC_Germany\Persistent_activity\high_freq_lup\mp_up_' f_names{d}];
%     print('-dpng',tname);
%     close

% plot_with_ci(sec_lag/Fsd,mp_med_first_up_trig,'b')
% hold on
% plot_with_ci(sec_lag/Fsd,mp_med_sec_up_trig,'k')
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (z)','FontSize',14)
% title('MP UP / LFP UP trig Avg','FontSize',16)
%     tname = ['C:\WC_Germany\Persistent_activity\high_freq_lup\mp_up_v_mp_down' f_names{d}];
%     print('-dpng',tname);
%     close

%     clear near_lfp_ind
%     clear prev_lfp_down
%     clear *up_trig
% 
% end

