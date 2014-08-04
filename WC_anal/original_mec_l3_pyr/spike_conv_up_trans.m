clear all

load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28
load C:\WC_Germany\JMM_analysis_pyr\t_90_data
load C:\WC_Germany\JMM_analysis_pyr\wcv_up_tau
load C:\WC_Germany\JMM_analysis_pyr\sig_fit_data
load C:\WC_Germany\JMM_analysis_pyr\sig_down_fit_data

Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

lags = -2*Fsd:2*Fsd;
maxlag = 2*Fsd;
% kern_dur = 10;
% half_point = kern_dur/2*Fsd;
% act_hist_kern = zeros(kern_dur*Fsd+1,1);
% rel_t = (1:half_point+1)/Fsd;
% tau = 0.1;
% act_hist_kern(half_point+1:end) = exp(-rel_t/tau);


for d = 1:length(dir_array)

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

    t_axis = 0:length(wcv);
    
    spike_id = round(spkid/8);
    t_90_wup{d} = round(t_90_wup{d}/8);
    spike_vec = hist(spike_id,t_axis);
    t_10_wdown{d} = round(t_10_wdown{d}/8);
    
%     hist_vec = conv(spike_vec,act_hist_kern);
%     hist_vec(1:floor(length(act_hist_kern)/2)) = [];
%     hist_vec(end-floor(length(act_hist_kern)/2)+1:end) = [];
    
    %% calculate wcv up triggered history average
%     num_wcv_ups = length(t_90_wup{d});
%     wup_spike_avg = zeros(num_wcv_ups,length(lags));
%     wup_hist_avg = zeros(num_wcv_ups,length(lags));
%     wup_wcv_avg = zeros(num_wcv_ups,length(lags));
%     wup_lf8_avg = zeros(num_wcv_ups,length(lags));
%     for i = 1:num_wcv_ups
%         if t_90_wup{d}(i) > maxlag && t_90_wup{d}(i) < length(wcv)-maxlag
%             wcv_down = down_trans{d}(find(down_trans{d} > t_90_wup{d}(i),1,'first'));
%             wcv_up_dur = wcv_down-t_90_wup{d}(i);
%             wup_spike_avg(i,:) = jmm_smooth_1d(spike_vec(t_90_wup{d}(i)-maxlag:t_90_wup{d}(i)+maxlag),10);
%             wup_hist_avg(i,:) = hist_vec(t_90_wup{d}(i)-maxlag:t_90_wup{d}(i)+maxlag);
%             if wcv_up_dur < maxlag
%                 wup_spike_avg(i,maxlag+wcv_up_dur:end) = nan;
%                 wup_hist_avg(i,maxlag+wcv_up_dur:end) = nan;
%             end
%             wup_wcv_avg(i,:) = wcv(t_90_wup{d}(i)-maxlag:t_90_wup{d}(i)+maxlag);
% %             wup_lf8_avg(i,:) = lf8(t_90_wup{d}(i)-maxlag:t_90_wup{d}(i)+maxlag);
%         else
%             wup_spike_avg(i,:) = nan;
%             wup_wcv_avg(i,:) = nan;
%             wup_hist_avg(i,:) = nan;
% %             wup_lf8_avg(i,:) = nan;
%         end
% 
%     end
% 
%     wup_spike_avg = 250*wup_spike_avg;

    
        %% calculate wcv down triggered spike average
    num_wcv_ups = length(t_10_wdown{d});
    wup_spike_avg = zeros(num_wcv_ups,length(lags));
%     wup_hist_avg = zeros(num_wcv_ups,length(lags));
    wup_wcv_avg = zeros(num_wcv_ups,length(lags));
    wup_lf8_avg = zeros(num_wcv_ups,length(lags));
    for i = 1:num_wcv_ups
        if t_10_wdown{d}(i) > maxlag && t_10_wdown{d}(i) < length(wcv)-maxlag
            wcv_down = down_trans{d}(find(down_trans{d} < t_10_wdown{d}(i),1,'last'));
            wcv_up_dur = t_10_wdown{d}(i) - wcv_down;
            wup_spike_avg(i,:) = jmm_smooth_1d(spike_vec(t_10_wdown{d}(i)-maxlag:t_10_wdown{d}(i)+maxlag),5);
%             wup_hist_avg(i,:) = hist_vec(t_10_wdown{d}(i)-maxlag:t_10_wdown{d}(i)+maxlag);
            if wcv_up_dur < maxlag
                wup_spike_avg(i,1:maxlag-wcv_up_dur) = nan;
%                 wup_hist_avg(i,maxlag+wcv_up_dur:end) = nan;
            end
            wup_wcv_avg(i,:) = wcv(t_10_wdown{d}(i)-maxlag:t_10_wdown{d}(i)+maxlag);
            wup_lf8_avg(i,:) = lf8(t_10_wdown{d}(i)-maxlag:t_10_wdown{d}(i)+maxlag);
        else
            wup_spike_avg(i,:) = nan;
            wup_wcv_avg(i,:) = nan;
%             wup_hist_avg(i,:) = nan;
            wup_lf8_avg(i,:) = nan;
        end

    end

    wup_spike_avg = 250*wup_spike_avg;

    
        m_wup_spike_avg(d,:) = nanmean(wup_spike_avg);
    u_wup_spike_avg(d,:) = m_wup_spike_avg(d,:)+2*nanstd(wup_spike_avg)./sqrt(sum(~isnan(wup_spike_avg)));
    d_wup_spike_avg(d,:) = m_wup_spike_avg(d,:)-2*nanstd(wup_spike_avg)./sqrt(sum(~isnan(wup_spike_avg)));
    %
%     z_spike(d,:) = m_wup_spike_avg(d,:)-nanmean(m_wup_spike_avg(d,:));
%     z_spike(d,:) = z_spike(d,:)/nanstd(m_wup_spike_avg(d,:));
    m_wup_wcv_avg(d,:) = nanmean(wup_wcv_avg);
    u_wup_wcv_avg(d,:) = m_wup_wcv_avg(d,:)+2*nanstd(wup_wcv_avg)/sqrt(num_wcv_ups);
    d_wup_wcv_avg(d,:) = m_wup_wcv_avg(d,:)-2*nanstd(wup_wcv_avg)/sqrt(num_wcv_ups);
    m_wup_lf8_avg(d,:) = nanmean(wup_lf8_avg);
    u_wup_lf8_avg(d,:) = m_wup_lf8_avg(d,:)+2*nanstd(wup_lf8_avg)/sqrt(num_wcv_ups);
    d_wup_lf8_avg(d,:) = m_wup_lf8_avg(d,:)-2*nanstd(wup_lf8_avg)/sqrt(num_wcv_ups);

    
    
    
%% get mean spike rate in each up state
% spike_rate = zeros(length(up_trans{d}),1);
% spike_count = spike_rate;
% for i = 1:length(up_trans{d})
%    
%     cur_spikes = find(spike_id > up_trans{d}(i) & spike_id < down_trans{d}(i));
%     spike_rate(i) = length(cur_spikes)/up_state_dur{d}(i);
%     spike_count(i)  = length(cur_spikes);
%     
% end

% %% get time since last spike for each up state transition
% tsls = zeros(length(up_trans{d}),1);
% for i = 1:length(up_trans{d})
%     
%    prev_spike = spike_id(find(spike_id < up_trans{d}(i),1,'last'));
%    if ~isempty(prev_spike)
%    tsls(i) = up_trans{d}(i)-prev_spike;
%    else
%        tsls(i) = nan;
%    end
% end
    
% plot(spike_rate(1:end-1),rltau_wdown{d}(2:end),'o')
% 
% figure
% plot(spike_count(1:end-1),rltau_wdown{d}(2:end),'o')

% pause

    
%% visualize results
    
% [tvals,tau_sort] = sort(rltau_wdown{d});
%     
%     subplot(2,1,1)
%     pcolor(wup_spike_avg(tau_sort,:));shading flat
%     subplot(2,1,2)
%     pcolor(wup_wcv_avg(tau_sort,:));shading flat

    m_spk(d,:) = nanmean(wup_spike_avg);
    u_spk = m_spk(d,:) + 2*nanstd(wup_spike_avg)/sqrt(num_wcv_ups);
    l_spk = m_spk(d,:) - 2*nanstd(wup_spike_avg)/sqrt(num_wcv_ups);
    figure
    subplot(2,1,1)
    plot(lags/Fsd,m_spk(d,:),'k','linewidth',2)
    hold on
    plot(lags/Fsd,u_spk,'k--')
    plot(lags/Fsd,l_spk,'k--')
    subplot(2,1,2)
    plot(lags/Fsd,mean(wup_wcv_avg))
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wcv_down_spike_avg' f_names{d}];
    print('-dpng',t_names)
    close all
% 
%     close all
% tsls = tsls';
% 
%  figure 
% good_pts = ~isnan(rltau_wup{d}) & ~isnan(tsls);
% [a,b] = corrcoef(rltau_wup{d}(good_pts),tsls(good_pts));
% plot(rltau_wup{d},tsls,'.')
% title(sprintf('tau_v_tsls C = %0.2g  P = %0.2g',a(2,1),b(2,1))) 
% 
% good_pts = ~isnan(rltau_wup{d}) & ~isnan(up_state_dur{d});
% [a,b] = corrcoef(rltau_wup{d}(good_pts),up_state_dur{d}(good_pts));
% figure
% plot(rltau_wup{d},up_state_dur{d},'.')
% title(sprintf('tau v next up dur C = %0.2g  P = %0.2g',a(2,1),b(2,1))) 
% % plot(rltau_wup{d},spike_rate,'o')
% % 
% % plot(rltau_wup{d}(2:end),spike_rate(1:end-1),'o')
% figure
% temp_r = rltau_wup{d}(2:end-1);
% temp_d = down_state_dur{d}(1:end-1);
% good_pts = ~isnan(temp_r) & ~isnan(temp_d);
% [a,b] = corrcoef(temp_r(good_pts),temp_d(good_pts));
% plot(temp_r,temp_d,'.')
% title(sprintf('tau v prev down state dur C = %0.2g  P = %0.2g',a(2,1),b(2,1))) 
% 
% figure
% temp_u = up_state_dur{d}(1:end-1);
% temp_d = down_state_dur{d};
% good_pts = ~isnan(temp_u) & ~isnan(temp_d);
% [a,b] = corrcoef(temp_u(good_pts),temp_d(good_pts));
% plot(temp_u,temp_d,'.')
% title(sprintf('up v next down C = %0.2g  P = %0.2g',a(2,1),b(2,1))) 
% 
% figure
% temp_u = up_state_dur{d}(2:end-1);
% temp_d = down_state_dur{d}(1:end-1);
% good_pts = ~isnan(temp_u) & ~isnan(temp_d);
% [a,b] = corrcoef(temp_u,temp_d);
% plot(temp_u(good_pts),temp_d(good_pts),'.')
% title(sprintf('up v prev down C = %0.2g  P = %0.2g',a(2,1),b(2,1))) 
% 
% figure
% temp_u = up_state_dur{d}(1:end-1);
% temp_u2 = up_state_dur{d}(2:end);
% good_pts = ~isnan(temp_u) & ~isnan(temp_u2);
% [a,b] = corrcoef(temp_u,temp_u2);
% plot(temp_u(good_pts),temp_u2(good_pts),'.')
% title(sprintf('up v prev up  C = %0.2g  P = %0.2g',a(2,1),b(2,1))) 
% 
% 
%     close all
%     
end


    m_s = nanmean(m_wup_spike_avg);
u_s = m_s+2*nanstd(m_wup_spike_avg)/sqrt(17);
d_s = m_s-2*nanstd(m_wup_spike_avg)/sqrt(17);
%



m_w = mean(m_wup_wcv_avg);
u_w = m_w+2*std(m_wup_wcv_avg)/sqrt(17);
d_w = m_w-2*std(m_wup_wcv_avg)/sqrt(17);
%
m_8 = mean(m_wup_lf8_avg);
u_8 = m_8+2*std(m_wup_lf8_avg)/sqrt(17);
d_8 = m_8-2*std(m_wup_lf8_avg)/sqrt(17);
%
subplot(2,1,1)
plot(lags/Fsd,m_s,'k','linewidth',2)
hold on
plot(lags/Fsd,u_s,'k--')
plot(lags/Fsd,d_s,'k--')
subplot(2,1,2)
plot(lags/Fsd,m_w,'linewidth',2)
hold on
plot(lags/Fsd,m_8,'r','linewidth',2)
legend('WCV Average','LF8 Average')
plot(lags/Fsd,u_w,'--')
plot(lags/Fsd,d_w,'--')
plot(lags/Fsd,u_8,'r--')
plot(lags/Fsd,d_8,'r--')
