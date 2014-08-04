clear all

load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28
load C:\WC_Germany\JMM_analysis_pyr\t_90_data
load C:\WC_Germany\JMM_analysis_pyr\wcv_up_tau

Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

lags = -2*Fsd:2*Fsd;
maxlag = 2*Fsd;

for d = 1:length(dir_array)

    cd(dir_array{d});
    pwd

    if ~exist('spike_time.mat')
        load spike_time_br
    else
        load spike_time
    end
    load used_data wcv_minus_spike lf8 wcv

    wcv = filtfilt(b,a,wcv_minus_spike);
    wcv = downsample(wcv,dsf);
    wcv = zscore(wcv);
    wcv = wcv';

    lf8 = filtfilt(b,a,lf8);
    lf8 = downsample(lf8,dsf);
    lf8 = zscore(lf8);
    lf8 = lf8';

    spike_id = round(spkid/8);
    t_90_wup{d} = round(t_90_wup{d}/8);

    %% calculate overall average up and down lfp
    %    ov_lfp_ups = zeros(length(up_trans8{d}),length(lags));
    %    for i = 1:length(up_trans8{d})
    %
    %        if up_trans8{d}(i) > maxlag & up_trans8{d}(i) < length(wcv) - maxlag
    %            ov_lfp_ups(i,:) = lf8(up_trans8{d}(i)-maxlag:up_trans8{d}(i)+maxlag);
    %        else
    %            ov_lfp_ups(i,:) = nan;
    %        end
    %
    %    end
    %    ov_lfp_downs = zeros(length(down_trans8{d}),length(lags));
    %    for i = 1:length(down_trans8{d})
    %        if down_trans8{d}(i) > maxlag & down_trans8{d}(i) < length(wcv) - maxlag
    %            ov_lfp_downs(i,:) = lf8(down_trans8{d}(i)-maxlag:down_trans8{d}(i)+maxlag);
    %        else
    %            ov_lfp_downs(i,:) = nan;
    %        end
    %    end

    %    m_lfp_up = nanmean(ov_lfp_ups);
    %    u_lfp_up = m_lfp_up+2*nanstd(ov_lfp_ups)/sqrt(length(up_trans8{d}));
    %    l_lfp_up = m_lfp_up-2*nanstd(ov_lfp_ups)/sqrt(length(up_trans8{d}));
    %
    %    m_lfp_down = nanmean(ov_lfp_downs);
    %    u_lfp_down = m_lfp_down+2*nanstd(ov_lfp_downs)/sqrt(length(down_trans8{d}));
    %    l_lfp_down = m_lfp_down-2*nanstd(ov_lfp_downs)/sqrt(length(down_trans8{d}));



%     %find lfp up transitions that occur during a wcv up state (wcv must have
%     %already been in up state for at least 0.5 s)
%     good_ups = [];
%     for i = 1:length(up_trans8{d})
% 
%         %find previous wcv up transition
%         prev_up = find(up_trans{d} < up_trans8{d}(i) - 0.3*Fsd,1,'last');
% 
%         if up_trans8{d}(i) < down_trans{d}(prev_up)
%             good_ups = [good_ups i];
%         end
% 
%     end
% 
%     %find lfp down transitions that occur during wcv up state
%     good_downs = [];
%     for i = 1:length(down_trans8{d})
% 
%         prev_up = find(up_trans{d} < down_trans8{d}(i) - 0.3*Fsd,1,'last');
% 
%         if down_trans8{d}(i) <= down_trans{d}(prev_up)
%             good_downs = [good_downs i];
%         end
% 
%     end

%%
%     %find wcv up transitions that occur during a lfp up state 
%     
%     good_ups = [];
%     for i = 1:length(up_trans{d})
% 
%         %find previous wcv up transition
%         prev_up = find(up_trans8{d} < up_trans{d}(i),1,'last');
% 
%         if up_trans{d}(i) < down_trans8{d}(prev_up)
%             good_ups = [good_ups i];
%         end
% 
%     end
% 
%     %find wcv down transitions that occur during lfp up state
%     good_downs = [];
%     for i = 1:length(down_trans{d})
% 
%         prev_up = find(up_trans8{d} < down_trans{d}(i),1,'last');
% 
%         if down_trans{d}(i) <= down_trans8{d}(prev_up)
%             good_downs = [good_downs i];
%         end
% 
%     end
%     
%     %find wcv up transitions that occur during a lfp down state
%        wup_ldown = [];
%        for i = 1:length(up_trans{d})
%     
%            %find previous wcv up transition
%            prev_down = find(down_trans8{d} < up_trans{d}(i),1,'last');
%     
%            if prev_down < length(up_trans8{d})
%                if up_trans{d}(i) < up_trans8{d}(prev_down+1)
%                    wup_ldown = [wup_ldown i];
%                end
%            end
%     
%        end
%        
%            %find wcv down transitions that occur during a lfp down state
%        wdown_ldown = [];
%        for i = 1:length(down_trans{d})
%     
%            prev_down = find(down_trans8{d} < down_trans{d}(i),1,'last');
%     
%            if prev_down < length(up_trans8{d})
%                if down_trans{d}(i) < up_trans8{d}(prev_down+1)
%                    wdown_ldown = [wdown_ldown i];
%                end
%            end
%     
%        end

%%
    
    
    
    %     %find lfp down transitions that occur during wcv up state and dont
    %     %trigger a down state trans in wcv
    %     ineff_downs = [];
    %     for i = 1:length(down_trans8{d})
    %
    %         prev_up = find(up_trans{d} < down_trans8{d}(i) - 0.3*Fsd,1,'last');
    %
    %         if down_trans8{d}(i) <= down_trans{d}(prev_up)-2*Fsd
    %             ineff_downs = [ineff_downs i];
    %         end
    %
    %     end






       %find lfp up transitions that occur during a wcv down state
       lup_wdown = [];
       lup_wdown_wind = [];
       for i = 1:length(up_trans8{d})
    
           %find previous wcv up transition
           prev_down = find(down_trans{d} < up_trans8{d}(i),1,'last');
    
           if prev_down < length(up_trans{d})
               if up_trans8{d}(i) < up_trans{d}(prev_down+1)
                   lup_wdown = [lup_wdown i];
                   lup_wdown_wind = [lup_wdown_wind prev_down+1];
               end
           end
    
       end

    %       %find lfp up transitions that occur during a wcv down state and don't
    %       %trigger an wcv up state within 1 s
    %    lup_wdown_i = [];
    %    for i = 1:length(up_trans8{d})
    %
    %        %find previous wcv up transition
    %        prev_down = find(down_trans{d} < up_trans8{d}(i),1,'last');
    %
    %        if prev_down < length(up_trans{d})
    %            if up_trans8{d}(i) < up_trans{d}(prev_down+1) - 1*Fsd
    %                lup_wdown_i = [lup_wdown_i i];
    %            end
    %        end
    %
    %    end


    %find lfp down transitions that occur during a wcv down state
    %    ldown_wdown = [];
    %    for i = 1:length(down_trans8{d})
    %
    %        find previous wcv up transition
    %        prev_down = find(down_trans{d} < down_trans8{d}(i),1,'last');
    %
    %        if prev_down < length(up_trans{d})
    %            if down_trans8{d}(i) < up_trans{d}(prev_down+1)
    %                ldown_wdown = [ldown_wdown i];
    %            end
    %        end
    %
    %    end


    %     bad_ups = setdiff([1:length(up_trans8{d})],good_ups);
    %     bad_downs = setdiff([1:length(down_trans8{d})],good_downs);

%         cor_up_trans = up_trans{d}(good_ups);
%         cor_down_trans = down_trans{d}(good_downs);

    %     %transitions that occur during a wcv down state
    %     cor_bup_trans = up_trans8{d}(bad_ups);
    %     cor_bdown_trans = down_trans8{d}(bad_downs);

    %     cor_idown_trans = down_trans8{d}(ineff_downs);

        cor_lup_wdown_trans = up_trans8{d}(lup_wdown);
    %     cor_ldown_wdown_trans = down_trans8{d}(ldown_wdown);
    %     cor_ilup_wdown_trans = up_trans8{d}(lup_wdown_i);
%         cor_wup_ldown_trans = up_trans{d}(wup_ldown);
%         cor_wdown_ldown_trans = down_trans{d}(wdown_ldown);

    %initialize trig spike counts
    %     up_trig_spk = zeros(length(cor_up_trans),length(lags));
    %     down_trig_spk = zeros(length(cor_down_trans),length(lags));
    %     wcv_up_trig = zeros(length(cor_up_trans),length(lags));
    %     lf8_up_trig = zeros(length(cor_up_trans),length(lags));
    %     wcv_down_trig = zeros(length(cor_down_trans),length(lags));
    %     lf8_down_trig = zeros(length(cor_down_trans),length(lags));
    %     wcv_bup_trig = zeros(length(cor_bup_trans),length(lags));
    %     wcv_bdown_trig = zeros(length(cor_bdown_trans),length(lags));

    %% calculate spike triggered lfp averages

    % %find spikes that occur during a lfp down state
    % spike_ldown = [];
    %
    % bad_spikes = find(spike_id< up_trans8{d}(1));
    % spike_id(bad_spikes) = [];
    %
    % for i = 1:length(spike_id)
    %
    %    %find previous lfp down transition
    %    prev_down = find(down_trans8{d} < spike_id(i),1,'last');
    %
    %    if prev_down < length(up_trans8{d})
    %       if spike_id(i) < up_trans8{d}(prev_down+1)
    %           spike_ldown = [spike_ldown i];
    %       end
    %    end
    %
    %
    % end
    %
    % %find spikes that occur during a lfp up state
    % spike_lup = [];
    %
    % for i = 1:length(spike_id)
    %
    %     %find previous lfp up transition
    %    prev_up = find(up_trans8{d} < spike_id(i),1,'last');
    %
    %    if prev_up < length(down_trans8{d})
    %       if spike_id(i) < down_trans8{d}(prev_up)
    %           spike_lup = [spike_lup i];
    %       end
    %    end
    %
    %
    % end
    %
    %
    % spk_trig_lfp_ldown = zeros(length(spike_ldown),length(lags));
    % spk_trig_lfp_lup = zeros(length(spike_lup),length(lags));
    %
    % for i = 1:length(spike_ldown)
    %
    %     cur_spike = spike_id(spike_ldown(i));
    %     prev_down = find(down_trans8{d}<cur_spike,1,'last');
    %     next_up = find(up_trans8{d} > cur_spike,1,'first');
    %     spk_trig_lfp_ldown(i,:) = lf8(cur_spike - maxlag:cur_spike+maxlag);
    %     if cur_spike - maxlag < down_trans8{d}(prev_down)
    %        spk_trig_lfp_ldown(i,1:(down_trans8{d}(prev_down) - cur_spike + maxlag )) = nan;
    %     end
    %     if cur_spike + maxlag > up_trans8{d}(next_up)
    %        spk_trig_lfp_ldown(i,up_trans8{d}(next_up)-cur_spike+maxlag:end) = nan;
    %     end
    %
    % end
    %
    % for i = 1:length(spike_lup)
    %
    %     cur_spike = spike_id(spike_lup(i));
    %     prev_up = find(up_trans8{d}<cur_spike,1,'last');
    %     next_down = find(down_trans8{d} > cur_spike,1,'first');
    %     spk_trig_lfp_lup(i,:) = lf8(cur_spike-maxlag:cur_spike+maxlag);
    %     if cur_spike-maxlag < up_trans8{d}(prev_up)
    %         spk_trig_lfp_lup(i,1:(up_trans8{d}(prev_up) - cur_spike+maxlag)) = nan;
    %     end
    %     if cur_spike + maxlag > down_trans8{d}(next_down);
    %         spk_trig_lfp_lup(i,down_trans8{d}(next_down)-cur_spike+maxlag:end) = nan;
    %     end
    %
    % end


    % %find number of real valued points in each matrix
    % num_used_pts = sum(~isnan(spk_trig_lfp_ldown));
    % m_spk_trig_lfp_ldown = nanmean(spk_trig_lfp_ldown);
    % u_spk_trig_lfp_ldown = m_spk_trig_lfp_ldown+2*nanstd(spk_trig_lfp_ldown)./sqrt(num_used_pts);
    % l_spk_trig_lfp_ldown = m_spk_trig_lfp_ldown-2*nanstd(spk_trig_lfp_ldown)./sqrt(num_used_pts);
    %
    % num_used_pts = sum(~isnan(spk_trig_lfp_lup));
    % m_spk_trig_lfp_lup = nanmean(spk_trig_lfp_lup);
    % u_spk_trig_lfp_lup = m_spk_trig_lfp_lup+2*nanstd(spk_trig_lfp_lup)./sqrt(num_used_pts);
    % l_spk_trig_lfp_lup = m_spk_trig_lfp_lup-2*nanstd(spk_trig_lfp_lup)./sqrt(num_used_pts);
    %
    % figure
    % plot(lags/Fsd,m_spk_trig_lfp_ldown,'linewidth',2)
    % hold on
    % plot(lags/Fsd,u_spk_trig_lfp_ldown,'--')
    % plot(lags/Fsd,l_spk_trig_lfp_ldown,'--')
    % xlabel('Time Lag (s)','FontSize',14)
    % ylabel('Average LFP amplitude (Zscore)','FontSize',14)
    % title('Spike Triggered LFP (LFP DOWN) Average')
    % xlim([-0.2 0.2])
    % grid
    % t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\spike_trig_lfp_down_zoom' f_names{d}];
    % print('-dpng',t_names)
    % close all
    %
    % figure
    % plot(lags/Fsd,m_spk_trig_lfp_lup,'linewidth',2)
    % hold on
    % plot(lags/Fsd,u_spk_trig_lfp_lup,'--')
    % plot(lags/Fsd,l_spk_trig_lfp_lup,'--')
    % xlabel('Time Lag (s)','FontSize',14)
    % ylabel('Average LFP amplitude (Zscore)','FontSize',14)
    % title('Spike Triggered LFP (LFP UP) Average')
    % xlim([-0.2 0.2])
    % grid
    % t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\spike_trig_lfp_up_zoom' f_names{d}];
    % print('-dpng',t_names)
    % close all


%     %% calculate wcv up triggered spike average
%     num_wcv_ups = length(t_90_wup{d});
%     wup_spike_avg = zeros(num_wcv_ups,length(lags));
%     wup_wcv_avg = zeros(num_wcv_ups,length(lags));
%     wup_lf8_avg = zeros(num_wcv_ups,length(lags));
%     for i = 1:num_wcv_ups
%         if t_90_wup{d}(i) > maxlag & t_90_wup{d}(i) < length(wcv)-maxlag
%             cur_hist = hist(spike_id-t_90_wup{d}(i),lags);
%             cur_hist(1) = 0;
%             cur_hist(end) = 0;
%             wcv_down = down_trans{d}(find(down_trans{d} > t_90_wup{d}(i),1,'first'));
%             wcv_up_dur = wcv_down-t_90_wup{d}(i);
%             if wcv_up_dur < maxlag
%                 cur_hist(maxlag+wcv_up_dur:end) = nan;
%             end
%                 wup_spike_avg(i,:) = jmm_smooth_1d(cur_hist,5);
% %             wup_spike_avg(i,:) = cur_hist;
%             wup_wcv_avg(i,:) = wcv(t_90_wup{d}(i)-maxlag:t_90_wup{d}(i)+maxlag);
%             wup_lf8_avg(i,:) = lf8(t_90_wup{d}(i)-maxlag:t_90_wup{d}(i)+maxlag);
%         else
%             wup_spike_avg(i,:) = nan;
%             wup_wcv_avg(i,:) = nan;
%             wup_lf8_avg(i,:) = nan;
%         end
% 
%     end
% 
%     wup_spike_avg = 250*wup_spike_avg;
    %
    % rel_down = zeros(num_wcv_ups,1);
    %
    % rel_down(2:end) = t_90_wup{d}(2:end)-down_trans{d}(1:end-1);

%     init_spikes = zeros(1,num_wcv_ups);
% spike_win = 50;
%     for i = 1:num_wcv_ups
% 
%         first_spike = find(spike_id > up_trans{d}(i)-50,1,'first');
% %         if spike_id(first_spike) > t_90_wup{d}(i)+50
% %             init_spikes(i) = length(find(spike_id > t_90_wup{d}(i) & spike_id < t_90_wup{d}(i)+spike_win));
% %         else
% 
%             if ~isempty(first_spike)
% 
%                 init_spikes(i) = 1+length(find(spike_id(first_spike+1:end) < spike_id(first_spike) + spike_win));
% 
%             else
%                 init_spikes(i) = nan;
%             end
% 
% %         end
%     end
% 
%     prev_downs = down_state_dur{d};
% % 
%     cur_tau = rltau_wup{d};
% 
% %     cur_tau(1) = [];
% % 
%     max_down = 3;
% % 
% % 
% tsls = zeros(num_wcv_ups,1);
% for i = 1:num_wcv_ups
%     
%     cur_u_tran = up_trans{d}(i);
%     prev_spike = find(spike_id < cur_u_tran,1,'last');
%     if~isempty(prev_spike)
%     tsls(i) = cur_u_tran-prev_spike;
%     else
%         tsls(i) = nan;
%     end
%     
% end
% 
% tsls = tsls';
% 
% 
%     z_tau = (cur_tau- nanmean(cur_tau))/nanstd(cur_tau);
%     good_pts = find(prev_downs < max_down & ~isnan(prev_downs) & ~isnan(cur_tau) & abs(z_tau) < 3);
%     % %
%     [a,b] = corrcoef(prev_downs(good_pts),cur_tau(good_pts));
%     plot(prev_downs(good_pts),cur_tau(good_pts),'o')
%     line([max_down max_down],[0 max(cur_tau(good_pts))],'Color','k','linewidt',2)
%     title(sprintf('Corr = %0.2g   P = %0.2g',a(2,1),b(2,1)))
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\up_slope\prev_down_v_slope_' f_names{d}];
%         print('-dpng',t_names)
%         close all


%     [dummy,tsort] = sort(cur_tau);
% 
%     good_pts = ~isnan(tsls) & ~isnan(cur_tau) & abs(z_tau) < 3;
%     [a,b] = corrcoef(cur_tau(good_pts),tsls(good_pts));
%     plot(cur_tau(good_pts),tsls(good_pts),'o')
%     title(sprintf('Corr = %0.2g   P = %0.2g',a(2,1),b(2,1)))
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\up_slope\tsls_v_slope_' f_names{d}];
%         print('-dpng',t_names)
%         close all

%     pcolor(wup_spike_avg(tsort,:));shading flat
% %     hold on
% %     plot(length(lags)/2-spike_lag(tsort),[1:length(spike_lag)],'r.')
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\up_slope\isi_v_slope_mat_' f_names{d}];
%         print('-dpng',t_names)
%         close all

%     m_wup_spike_avg(d,:) = nanmean(wup_spike_avg);
%     u_wup_spike_avg(d,:) = m_wup_spike_avg(d,:)+2*nanstd(wup_spike_avg)./sqrt(sum(~isnan(wup_spike_avg)));
%     d_wup_spike_avg(d,:) = m_wup_spike_avg(d,:)-2*nanstd(wup_spike_avg)./sqrt(sum(~isnan(wup_spike_avg)));
%     %
% %     z_spike(d,:) = wup_spik-nanmean(m_wup_spike_avg(d,:));
% %     z_spike(d,:) = z_spike(d,:)/nanstd(m_wup_spike_avg(d,:));
%     m_wup_wcv_avg(d,:) = nanmean(wup_wcv_avg);
%     u_wup_wcv_avg(d,:) = m_wup_wcv_avg(d,:)+2*nanstd(wup_wcv_avg)/sqrt(num_wcv_ups);
%     d_wup_wcv_avg(d,:) = m_wup_wcv_avg(d,:)-2*nanstd(wup_wcv_avg)/sqrt(num_wcv_ups);
%     m_wup_lf8_avg(d,:) = nanmean(wup_lf8_avg);
%     u_wup_lf8_avg(d,:) = m_wup_lf8_avg(d,:)+2*nanstd(wup_lf8_avg)/sqrt(num_wcv_ups);
%     d_wup_lf8_avg(d,:) = m_wup_lf8_avg(d,:)-2*nanstd(wup_lf8_avg)/sqrt(num_wcv_ups);
% 
    %     %for ineffective wup conditional lfp downs
    %     if length(cor_idown_trans) > 0
    %         wcv_idown_trig = zeros(length(cor_idown_trans),length(lags));
    %         idown_trig_spk = zeros(length(cor_idown_trans),length(lags));
    %
    %         for i = 1:length(cor_idown_trans)
    %             if cor_idown_trans(i) > maxlag & cor_idown_trans(i) < length(wcv) - maxlag
    %                 cur_hist = hist(spike_id-cor_idown_trans(i),lags);
    %                 cur_hist(1) = 0;
    %                 cur_hist(end) = 0;
    %                 idown_trig_spk(i,:) = cur_hist;
    %                 wcv_idown_trig(i,:) = wcv(cor_idown_trans(i)-maxlag:cor_idown_trans(i)+maxlag);
    %             else
    %                 idown_trig_spk(i,:) = nan;
    %                 wcv_idown_trig(i,:) = nan;
    %             end
    %         end
    %     end


    %     for wcv down conditional lfp ups
        if length(cor_lup_wdown_trans) > 0
            lf8_wup_ldown_trig = zeros(length(cor_lup_wdown_trans),length(lags));
            wcv_wup_ldown_trig = zeros(length(cor_lup_wdown_trans),length(lags));
            for i = 1:length(cor_lup_wdown_trans)
                if cor_lup_wdown_trans(i) > maxlag & cor_lup_wdown_trans(i) < length(wcv) - maxlag
                    wcv_wup_ldown_trig(i,:) = wcv(cor_lup_wdown_trans(i)-maxlag:cor_lup_wdown_trans(i)+maxlag);
                    lf8_wup_ldown_trig(i,:) = lf8(cor_lup_wdown_trans(i)-maxlag:cor_lup_wdown_trans(i)+maxlag);
                else
                    wcv_wup_ldown_trig(i,:) = nan;
                    lf8_wup_ldown_trig(i,:) = nan;
                end
            end
        end
[dummy,up_sort] = sort(up_state_dur{d}(lup_wdown_wind));

subplot(2,1,1)
pcolor(lags,1:size(wcv_wup_ldown_trig,1),wcv_wup_ldown_trig(up_sort,:));shading flat
subplot(2,1,2)
pcolor(lags,1:size(lf8_wup_ldown_trig,1),lf8_wup_ldown_trig(up_sort,:));shading flat
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\lup_wdown_updursort' f_names{d}];
%         print('-dpng',t_names)
%         close all

        lag_prev_lup = zeros(size(wcv_wup_ldown_trig,1),1);
        lag_prev_lup = (up_trans{d}(lup_wdown_wind)-cor_lup_wdown_trans)/Fsd;
        
        scatter_with_cor(lag_prev_lup,up_state_dur{d}(lup_wdown_wind));
        t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\lup_wdown_updursort_lagcor' f_names{d}];
        print('-dpng',t_names)
        close all

    %      %for ineffective wcv down conditional lfp ups
    %     if length(cor_ilup_wdown_trans) > 0
    %         wcv_ilup_wdown_trig = zeros(length(cor_ilup_wdown_trans),length(lags));
    %         for i = 1:length(cor_ilup_wdown_trans)
    %             if cor_ilup_wdown_trans(i) > maxlag & cor_ilup_wdown_trans(i) < length(wcv) - maxlag
    %                 wcv_ilup_wdown_trig(i,:) = wcv(cor_ilup_wdown_trans(i)-maxlag:cor_ilup_wdown_trans(i)+maxlag);
    %             else
    %                 wcv_ilup_wdown_trig(i,:) = nan;
    %             end
    %         end
    %     end


%     for wcv down conditional lfp downs
%         if length(cor_wdown_ldown_trans) > 0
%             wcv_wdown_ldown_trig = zeros(length(cor_wdown_ldown_trans),length(lags));
%             lf8_wdown_ldown = wcv_wdown_ldown_trig;
%             for i = 1:length(cor_wdown_ldown_trans)
%                 if cor_wdown_ldown_trans(i) > maxlag & cor_wdown_ldown_trans(i) < length(wcv) - maxlag
%                     wcv_wdown_ldown_trig(i,:) = wcv(cor_wdown_ldown_trans(i)-maxlag:cor_wdown_ldown_trans(i)+maxlag);
%                     lf8_wdown_ldown(i,:) = lf8(cor_wdown_ldown_trans(i)-maxlag:cor_wdown_ldown_trans(i)+maxlag);
%                 else
%                     wcv_wdown_ldown_trig(i,:) = nan;
%                     lf8_wdown_ldown(i,:) = nan;
%                 end
%             end
%         end
% 
%         for i = 1:length(cor_up_trans)
%             if cor_up_trans(i) > maxlag & cor_up_trans(i) < length(wcv)-maxlag
%     
%     %         cur_hist = hist(spike_id-cor_up_trans(i),lags);
%     %         cur_hist(1) = 0;
%     %         cur_hist(end) = 0;
%     %         up_trig_spk(i,:) = jmm_smooth_1d(cur_hist,10);
%             wcv_up_trig(i,:) = wcv(cor_up_trans(i)-maxlag:cor_up_trans(i)+maxlag);
%             lf8_up_trig(i,:) = lf8(cor_up_trans(i)-maxlag:cor_up_trans(i)+maxlag);
%     
%             else
%     %             up_trig_spk(i,:) = nan;
%                 wcv_up_trig(i,:) = nan;
%                 lf8_up_trig(i,:) = nan;
%             end
%         end
% 
%         for i = 1:length(cor_down_trans)
%             if cor_down_trans(i) > maxlag & cor_down_trans(i) < length(wcv) - maxlag
%     %         cur_hist = hist(spike_id-cor_down_trans(i),lags);
%     %         cur_hist(1) = 0;
%     %         cur_hist(end) = 0;
%     %         down_trig_spk(i,:) = cur_hist;
%             wcv_down_trig(i,:) = wcv(cor_down_trans(i)-maxlag:cor_down_trans(i)+maxlag);
%             lf8_down_trig(i,:) = lf8(cor_down_trans(i)-maxlag:cor_down_trans(i)+maxlag);
%     
%             else
%     %             down_trig_spk(i,:) = nan;
%                 wcv_down_trig(i,:) = nan;
%                 lf8_down_trig(i,:) = nan;
%             end
%         end
% 
    %     for i = 1:length(cor_bup_trans)
    %         if cor_bup_trans(i) > maxlag & cor_bup_trans(i) < length(wcv) - maxlag
    %         wcv_bup_trig(i,:) = wcv(cor_bup_trans(i)-maxlag:cor_bup_trans(i)+maxlag);
    %         end
    %     end
    %
    %     for i = 1:length(cor_bdown_trans)
    %         if cor_bdown_trans(i) > maxlag & cor_bdown_trans < length(wcv) - maxlag
    %         wcv_bdown_trig(i,:) = wcv(cor_bdown_trans(i)-maxlag:cor_bdown_trans(i)+maxlag);
    %         end
    %     end

    %% plot wcv up triggered spikes
    % wup_spike_avg(d,:) = wup_spike_avg(d,:)*250;
    % subplot(2,1,1)
    % plot(lags/Fsd,m_wup_spike_avg(d,:),'k','linewidth',2)
    % up_lim = max(m_wup_spike_avg(d,:))+0.01;
    % hold on
    % plot(lags/Fsd,u_wup_spike_avg(d,:),'k--')
    % plot(lags/Fsd,d_wup_spike_avg(d,:),'k--')
    % ylim([0 up_lim])
    % title('Spike Average')
    % subplot(2,1,2)
    % plot(lags/Fsd,m_wup_wcv_avg(d,:),'linewidth',2)
    % hold on
    % plot(lags/Fsd,u_wup_wcv_avg(d,:),'--')
    % plot(lags/Fsd,d_wup_wcv_avg(d,:),'--')
    % plot(lags/Fsd,m_wup_lf8_avg(d,:),'r','linewidth',2)
    % plot(lags/Fsd,u_wup_lf8_avg(d,:),'r--')
    % plot(lags/Fsd,d_wup_lf8_avg(d,:),'r--')
    % legend('WCV','LF8')
    % title('WCV/LF8 Average')
    %     t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wup_spike_avg_t_90' f_names{d}];
    %     print('-dpng',t_names)
    %     close all


    %     subplot(3,1,1)
    %     pcolor(lags/Fsd,[1:size(wup_wcv_avg,1)],wup_wcv_avg);shading flat
    %     title('WCV')
    %     subplot(3,1,2)
    %     pcolor(lags/Fsd,[1:size(wup_wcv_avg,1)],wup_lf8_avg);shading flat
    %         title('LF8')
    %     subplot(3,1,3)
    %     hold on
    % for i = 1:size(wup_wcv_avg,1)
    %     spike_locs = find(wup_spike_avg(i,:)>0);
    % % plot(spike_locs/Fsd,ones(size(spike_locs))*i,'.','MarkerSize',6)
    % sm_spikes(i,:) = jmm_smooth_1d(wup_spike_avg(i,:),10);
    % end
    %     pcolor(lags/Fsd,[1:size(wup_wcv_avg,1)],sm_spikes);shading flat
    %     ylim([0 size(wup_wcv_avg,1)])
    %     title('Spike Rate')
    %         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wup_spikes\' f_names{d}];
    %     print('-dpng',t_names)
    %     close all
    %     clear sm_spikes


    % % %% plot all wcv up lfp up states
%     if size(wcv_up_trig,1) > 1
%     Fig = figure(1)
%         clf
%         set(Fig,'PaperUnits','centimeters');
%         set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%         set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%         subplot(2,1,1)
%         pcolor(lags/Fsd,[1:size(wcv_up_trig,1)],wcv_up_trig);shading flat
%         colorbar
%         subplot(2,1,2)
%         pcolor(lags/Fsd,[1:size(lf8_up_trig,1)],lf8_up_trig);shading flat
%         colorbar
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wtrig\wup_trig_lup_' f_names{d}];
%         print('-dpng',t_names)
%         close all
%     end
%     %% plot all wcv up lfp down states
%     if size(wcv_down_trig,1) > 1
%     Fig = figure(1)
%         clf
%         set(Fig,'PaperUnits','centimeters');
%         set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%         set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%         subplot(2,1,1)
%         pcolor(lags/Fsd,[1:size(wcv_down_trig,1)],wcv_down_trig);shading flat
%         colorbar
%         subplot(2,1,2)
%         pcolor(lags/Fsd,[1:size(wcv_down_trig,1)],lf8_down_trig);shading flat
%         colorbar
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wtrig\wdown_trig_lup_' f_names{d}];
%         print('-dpng',t_names)
%         close all
%     end
    % %% plot all wcv up lfp idown states
    % % if size(wcv_idown_trig,1) > 1
    % % Fig = figure(1)
    % %     clf
    % %     set(Fig,'PaperUnits','centimeters');
    % %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    % %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    % %     pcolor(lags/Fsd,[1:size(wcv_idown_trig)],wcv_idown_trig);shading flat
    % %     colorbar
    % %     t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wup_ildown\' f_names{d}];
    % %     print('-dpng',t_names)
    % %     close all
    % % end
    %
    %
    %% plot all wcv down lfp down states
%     if size(wcv_wdown_ldown_trig,1) > 1
%     Fig = figure(1)
%         clf
%         set(Fig,'PaperUnits','centimeters');
%         set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%         set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%         subplot(2,1,1)
%         pcolor(lags/Fsd,[1:size(wcv_wdown_ldown_trig,1)],wcv_wdown_ldown_trig);shading flat
%         colorbar
%         subplot(2,1,2)
%         pcolor(lags/Fsd,[1:size(wcv_wdown_ldown_trig,1)],lf8_wdown_ldown);shading flat
%         colorbar
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wtrig\wdown_trig_ldown_' f_names{d}];
%         print('-dpng',t_names)
%         close all
%     end
%     %
%     %% plot all wcv down lfp up states
%     if size(wcv_wup_ldown_trig,1) > 1
%     Fig = figure(1)
%         clf
%         set(Fig,'PaperUnits','centimeters');
%         set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%         set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%         subplot(2,1,1)
%         pcolor(lags/Fsd,[1:size(wcv_wup_ldown_trig,1)],wcv_wup_ldown_trig);shading flat
%         colorbar
%         subplot(2,1,2)
%         pcolor(lags/Fsd,[1:size(wcv_wup_ldown_trig,1)],lf8_wup_ldown_trig);shading flat
%         colorbar
%         t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wtrig\wup_trig_ldown_' f_names{d}];
%         print('-dpng',t_names)
%         close all
%     end
    %
    % % %% plot averages for spikes
    % %
    % %     %convert spike trig avgs to units of rate
    % %     up_trig_spk = up_trig_spk*250;
    % %     down_trig_spk = down_trig_spk*250;
    % %     idown_trig_spk = idown_trig_spk*250;
    % %
    % %     Fig = figure(1)
    % %     clf
    % %     set(Fig,'PaperUnits','centimeters');
    % %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    % %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    % %     plot(lags/Fsd,jmm_smooth_1d(mean(up_trig_spk),5),'linewidth',2)
    % %     hold on
    % %     plot(lags/Fsd,jmm_smooth_1d(mean(down_trig_spk),5),'r','linewidth',2)
    % %     plot(lags/Fsd,jmm_smooth_1d(mean(idown_trig_spk),5),'k','linewidth',2)
    % %     legend('Up Triggered','Down Triggered','I-Down Triggered')
    % %     t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\spike_' f_names{d}];
    % %     print('-dpng',t_names)
    % %     close all
    % %
    % %
    % %% PLot averages for wcv up state
    %     num_wup_trig = size(wcv_up_trig,1);
    %     m_wup_trig = nanmean(wcv_up_trig);
    %     u_wup_trig = m_wup_trig+2*nanstd(wcv_up_trig)/sqrt(num_wup_trig);
    %     l_wup_trig = m_wup_trig-2*nanstd(wcv_up_trig)/sqrt(num_wup_trig);
    %
    % %     m_lup_trig = mean(lf8_up_trig);
    % %     u_lup_trig = m_lup_trig + 2*std(lf8_up_trig)/sqrt(num_wup_trig);
    % %     l_lup_trig = m_lup_trig - 2*std(lf8_up_trig)/sqrt(num_wup_trig);
    % %
    %     num_wdown_trig = size(wcv_down_trig,1);
    %     m_wdown_trig = nanmean(wcv_down_trig);
    %     u_wdown_trig = m_wdown_trig + 2*nanstd(wcv_down_trig)/sqrt(num_wdown_trig);
    %     l_wdown_trig = m_wdown_trig - 2*nanstd(wcv_down_trig)/sqrt(num_wdown_trig);
    %
    % %     num_widown_trig = size(wcv_idown_trig,1);
    % %     m_widown_trig = nanmean(wcv_idown_trig,1);
    % %     u_widown_trig = m_widown_trig+2*nanstd(wcv_idown_trig)/sqrt(num_widown_trig);
    % %     l_widown_trig = m_widown_trig-2*nanstd(wcv_idown_trig)/sqrt(num_widown_trig);
    %
    % %     Fig = figure(1)
    % %     clf
    % %     set(Fig,'PaperUnits','centimeters');
    % %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    % %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    % %     plot(lags/Fsd,m_wup_trig,'linewidth',2)
    % %     hold on
    % %     plot(lags/Fsd,m_wdown_trig,'r','linewidth',2)
    % %     plot(lags/Fsd,m_widown_trig,'k','linewidth',2)
    % %     legend('WCV UP - LFP UP','WCV UP - LFP DOWN','WCV UP - I- LFP DOWN')
    % %     plot(lags/Fsd,u_wup_trig,'--')
    % %     plot(lags/Fsd,l_wup_trig,'--')
    % %     plot(lags/Fsd,u_wdown_trig,'r--')
    % %     plot(lags/Fsd,l_wdown_trig,'r--')
    % %     plot(lags/Fsd,u_widown_trig,'k--')
    % %     plot(lags/Fsd,l_widown_trig,'k--')
    % %     title(sprintf('%0.3g  %0.3g  %0.3g',num_wup_trig,num_wdown_trig,num_widown_trig))
    % %     t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wcv_up_' f_names{d}];
    % %     print('-dpng',t_names)
    % %     close all
    %
    % %% for wcv down states
    %
    % num_lup_wdown = size(wcv_lup_wdown_trig,1);
    % m_lup_wdown = nanmean(wcv_lup_wdown_trig);
    % u_lup_wdown = m_lup_wdown+2*nanstd(wcv_lup_wdown_trig)/sqrt(num_lup_wdown);
    % l_lup_wdown = m_lup_wdown-2*nanstd(wcv_lup_wdown_trig)/sqrt(num_lup_wdown);
    %
    % % lfp_m_lup_wdown = mean(lf8_lup_wdown_trig);
    % % lfp_u_lup_wdown = lfp_m_lup_wdown+2*std(lf8_lup_wdown_trig)/sqrt(num_lup_wdown);
    % % lfp_l_lup_wdown = lfp_m_lup_wdown-2*std(lf8_lup_wdown_trig)/sqrt(num_lup_wdown);
    %
    % num_ldown_wdown = size(wcv_ldown_wdown_trig,1);
    % m_ldown_wdown = nanmean(wcv_ldown_wdown_trig);
    % u_ldown_wdown = m_ldown_wdown + 2*nanstd(wcv_ldown_wdown_trig)/sqrt(num_ldown_wdown);
    % l_ldown_wdown = m_ldown_wdown - 2*nanstd(wcv_ldown_wdown_trig)/sqrt(num_ldown_wdown);
    %
    % % num_ilup_wdown = size(wcv_ilup_wdown_trig,1);
    % % m_ilup_wdown = nanmean(wcv_ilup_wdown_trig);
    % % u_ilup_wdown = m_ilup_wdown + 2*nanstd(wcv_ilup_wdown_trig)/sqrt(num_ilup_wdown);
    % % l_ilup_wdown = m_ilup_wdown - 2*nanstd(wcv_ilup_wdown_trig)/sqrt(num_ilup_wdown);
    % %
    % %     Fig = figure(1)
    % %     clf
    % %     set(Fig,'PaperUnits','centimeters');
    % %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    % %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    % %     plot(lags/Fsd,m_lup_wdown,'linewidth',2)
    % %     hold on
    % %     plot(lags/Fsd,m_ldown_wdown,'r','linewidth',2)
    % %     plot(lags/Fsd,m_ilup_wdown,'k','linewidth',2)
    % %     legend('WCV DOWN - LFP UP','WCV DOWN - LFP DOWN','WCV DOWN - I- LFP UP')
    % %     plot(lags/Fsd,u_lup_wdown,'--')
    % %     plot(lags/Fsd,l_lup_wdown,'--')
    % %     plot(lags/Fsd,u_ldown_wdown,'r--')
    % %     plot(lags/Fsd,l_ldown_wdown,'r--')
    % %     plot(lags/Fsd,u_ilup_wdown,'k--')
    % %     plot(lags/Fsd,l_ilup_wdown,'k--')
    % %     title(sprintf('%0.3g  %0.3g   %0.3g',num_lup_wdown,num_ldown_wdown,num_ilup_wdown))
    % %     t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\wcv_down_' f_names{d}];
    % %     print('-dpng',t_names)
    % %     close all
    %
    %
    %
    %     Fig = figure(1)
    %     clf
    %     set(Fig,'PaperUnits','centimeters');
    %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    %     plot(lags/Fsd,m_wup_trig,'linewidth',2)
    %     hold on
    %     plot(lags/Fsd,m_lup_wdown,'r','linewidth',2)
    %     plot(lags/Fsd,m_lfp_up,'k','linewidth',2)
    %    	legend('WCV UP cond, LFP UP','WCV DOWN cond, LFP UP','LFP UP')
    %     plot(lags/Fsd,u_wup_trig,'--')
    %     plot(lags/Fsd,l_wup_trig,'--')
    %     plot(lags/Fsd,u_lup_wdown,'r--')
    %     plot(lags/Fsd,l_lup_wdown,'r--')
    %     plot(lags/Fsd,u_lfp_up,'k--')
    %     plot(lags/Fsd,l_lfp_up,'k--')
    %     grid
    %     t_names = ['C:\WC_Germany\Layer5\State_Trig\LFP_ups\' f_names{d}];
    %     print('-dpng',t_names)
    %     close all
    %
    %     Fig = figure(1)
    %     clf
    %     set(Fig,'PaperUnits','centimeters');
    %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    %     plot(lags/Fsd,m_ldown_wdown,'linewidth',2)
    %     hold on
    %     plot(lags/Fsd,m_wdown_trig,'r','linewidth',2)
    %     plot(lags/Fsd,m_lfp_down,'k','linewidth',2)
    %    	legend('WCV DOWN cond, LFP DOWN','WCV UP cond, LFP DOWN','LFP DOWN')
    %     plot(lags/Fsd,u_ldown_wdown,'--')
    %     plot(lags/Fsd,l_ldown_wdown,'--')
    %     plot(lags/Fsd,u_wdown_trig,'r--')
    %     plot(lags/Fsd,l_wdown_trig,'r--')
    %     plot(lags/Fsd,u_lfp_down,'k--')
    %     plot(lags/Fsd,l_lfp_down,'k--')
    %     grid
    %     t_names = ['C:\WC_Germany\Layer5\State_Trig\LFP_downs\' f_names{d}];
    %     print('-dpng',t_names)
    %     close all


end

% m_s = nanmean(z_spike);
% u_s = m_s+2*nanstd(z_spike)/sqrt(17);
% d_s = m_s-2*nanstd(z_spike)/sqrt(17);
%
% m_s = nanmean(m_wup_spike_avg);
% u_s = m_s+2*nanstd(m_wup_spike_avg)/sqrt(17);
% d_s = m_s-2*nanstd(m_wup_spike_avg)/sqrt(17);
% 
% 
% 
% m_w = mean(m_wup_wcv_avg);
% u_w = m_w+2*std(m_wup_wcv_avg)/sqrt(17);
% d_w = m_w-2*std(m_wup_wcv_avg)/sqrt(17);
% %
% m_8 = mean(m_wup_lf8_avg);
% u_8 = m_8+2*std(m_wup_lf8_avg)/sqrt(17);
% d_8 = m_8-2*std(m_wup_lf8_avg)/sqrt(17);
% %
% subplot(2,1,1)
% plot(lags/Fsd,m_s,'k','linewidth',2)
% hold on
% plot(lags/Fsd,u_s,'k--')
% plot(lags/Fsd,d_s,'k--')
% subplot(2,1,2)
% plot(lags/Fsd,m_w,'linewidth',2)
% hold on
% plot(lags/Fsd,m_8,'r','linewidth',2)
% legend('WCV Average','LF8 Average')
% plot(lags/Fsd,u_w,'--')
% plot(lags/Fsd,d_w,'--')
% plot(lags/Fsd,u_8,'r--')
% plot(lags/Fsd,d_8,'r--')
