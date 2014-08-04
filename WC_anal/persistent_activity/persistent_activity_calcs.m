clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
% load C:\WC_Germany\Persistent_activity\sigmoid_fit\sig_fit_data
load C:\WC_Germany\Persistent_activity\lf8_period_f_data2
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;
Fsd = 2016/dsf;
backlag = 5*Fsd;
forwardlag = 10*Fsd
lags = -backlag:forwardlag;
% dlags = -backlag:backlag;

for d = 1:length(dir_array)
    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8_f = filtfilt(b,a,lf8);
    wcv_f = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8_f,dsf);
    wcv_d = downsample(wcv_f,dsf);

    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);

    wcv_up_ctrig_mat{d} = zeros(length(synch_ups{d}),length(lags));
    lf8_up_ctrig_mat{d} = zeros(length(synch_ups{d}),length(lags));
%     wcv_down_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
%     lf8_down_ctrig_mat = zeros(length(synch_ups{d}),length(lags));

    %initialize
    secondary_lfp_states = cell(length(synch_ups{d}),1);
    mp_down_state = zeros(length(synch_ups{d}),1);
    lfp_pup = zeros(1,length(synch_ups{d}));
    lfp_pdown = lfp_pup;
    lfp_nup = lfp_pup;
    lfp_pupdur = lfp_pup;
    lfp_period_dur{d} = lfp_pup;
    mean_lfp_up = lfp_pup;
    bad_rows = [];
    
    lfp_up_lag{d} = zeros(size(synch_ups{d}));
    
    for i = 1:length(synch_ups{d})

        %find duration of MP up state in units of LFP periods
        lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i)))-lf8_period_f{d}(up_trans{d}(synch_ups{d}(i)));

        %find first lfp down state after MP up transition
        first_in_lfp_down = find(down_trans8{d} > up_trans{d}(synch_ups{d}(i)),1,'first');

        if ~isempty(first_in_lfp_down)
            %find all lfp up transitions that occur before mp comes back down
            secondary_lfp_states{i} = find(down_trans8{d} < down_trans{d}(synch_ups{d}(i)) ...
                & up_trans8{d} > down_trans8{d}(first_in_lfp_down));

%             %if there are secondary lfp up states, set up portions to nans
%             %for visualization purposes
%             if ~isempty(secondary_lfp_states{i})
%                 for q = 1:length(secondary_lfp_states{i})
%                     lf8_d(up_trans8{d}(secondary_lfp_states{i}(q)):down_trans8{d}(secondary_lfp_states{i}(q))) = nan;
%                 end
%             end

        end

        %acquire MP up related lfp transitions
        [dummy,near_lfp_up] = min(abs(up_trans{d}(synch_ups{d}(i)) - up_trans8{d}(synch_ups8{d})));
        lfp_up_lag{d}(i) = up_trans{d}(synch_ups{d}(i))-up_trans8{d}(synch_ups8{d}(near_lfp_up));
        
        prev_lfp_up = find(up_trans8{d}<up_trans{d}(synch_ups{d}(i)),1,'last');
        prev_lfp_down = find(down_trans8{d}<down_trans{d}(synch_ups{d}(i)),1,'last');
        next_lfp_up = find(up_trans8{d}>down_trans{d}(synch_ups{d}(i)),1,'first');
        if ~isempty(prev_lfp_up)
            up_lag = up_trans{d}(synch_ups{d}(i))-up_trans8{d}(prev_lfp_up);
            cor_lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i))-up_lag)-lf8_period_f{d}(up_trans{d}(synch_ups{d}(i))-up_lag);
            lfp_pup(i) = up_trans8{d}(prev_lfp_up);
            lfp_pupdur(i) = up_state_dur8{d}(prev_lfp_up);
        else
            lfp_pup(i) = nan;
            lfp_pupdur(i) = nan;
        end
        if ~isempty(prev_lfp_down)
            lfp_pdown(i) = down_trans8{d}(prev_lfp_down);
        else
            lfp_pdown(i) = nan;
        end
        if ~isempty(next_lfp_up)
            lfp_nup(i) = up_trans8{d}(next_lfp_up);
        else
            lfp_nup(i) = nan;
        end

        %calculate average LFP up state duration for this MP up state
        if isempty(secondary_lfp_states{i}) & ~isempty(prev_lfp_up)
            mean_lfp_up(i) = up_state_dur8{d}(prev_lfp_up);
        elseif ~isempty(prev_lfp_up)
            mean_lfp_up(i) = mean(up_state_dur8{d}([secondary_lfp_states{i} prev_lfp_up]));
        else
            mean_lfp_up(i) = nan;
        end
        
        %check if MP down transition occured during and LFP down state
        if ~isempty(prev_lfp_down) & prev_lfp_down < length(up_trans8{d})
            mp_down_state(i) = up_trans8{d}(prev_lfp_down+1) > down_trans{d}(synch_ups{d}(i));
        else
            mp_down_state(i) = nan;
        end

        %calculate mup triggered averages
        if up_trans{d}(synch_ups{d}(i)) > backlag & up_trans{d}(synch_ups{d}(i)) ...
                < length(wcv_d)-forwardlag
            wcv_up_ctrig_mat{d}(i,:) = wcv_d(up_trans{d}(synch_ups{d}(i))-backlag: ...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
            lf8_up_ctrig_mat{d}(i,:) = lf8_d(up_trans{d}(synch_ups{d}(i))-backlag: ...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
        else
            wcv_up_ctrig_mat{d}(i,:) = nan;
            lf8_up_ctrig_mat{d}(i,:) = nan;
            bad_rows = [bad_rows i];
        end

                %calculate mdown triggered averages
%         if down_trans{d}(synch_ups{d}(i)) > backlag & down_trans{d}(synch_ups{d}(i)) ...
%                 < length(wcv_d)-forwardlag
%             wcv_down_ctrig_mat(i,:) = wcv_d(down_trans{d}(synch_ups{d}(i))-backlag: ...
%                 down_trans{d}(synch_ups{d}(i))+forwardlag);
%             lf8_down_ctrig_mat(i,:) = lf8_d(down_trans{d}(synch_ups{d}(i))-backlag: ...
%                 down_trans{d}(synch_ups{d}(i))+forwardlag);
%         else
%             wcv_down_ctrig_mat(i,:) = nan;
%             lf8_down_ctrig_mat(i,:) = nan;
%         end
%         
    end

    tsld = (down_trans{d}(synch_ups{d})-lfp_pdown)/Fsd; %time from mp down since last lfp down state
    ttnu = (lfp_nup-down_trans{d}(synch_ups{d}))/Fsd; %time from mp down to next lfp up state
    tslu = (up_trans{d}(synch_ups{d})-lfp_pup)/Fsd; %time from mp up since last lfp up

    mp_up_dur_lfp_up_units{d} = up_state_dur{d}(synch_ups{d})./mean_lfp_up;

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

%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    %
%     subplot(2,1,1)
%     pcolor(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,wcv_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))...
%         /length(synch_ups{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     xlim([-1 1])
%     subplot(2,1,2)
%     pcolor(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d})...
%         ,lf8_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     xlim([-1 1])
%     tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\zoom_mp_up_trig_avg' f_names{d}];
%     print('-dpng',tname);
%     close

    %calculate fraction of MP up states with secondary states
    sec_lfp_check{d} = zeros(length(synch_ups{d}),1);
    num_sec_states{d} = zeros(length(synch_ups{d}),1);

    for i = 1:length(synch_ups{d})

        sec_lfp_check{d}(i) = ~isempty(secondary_lfp_states{i});
        num_sec_states{d}(i) = length(secondary_lfp_states{i});

    end

%     mp_down_frac(d) = nansum(mp_down_state)/length(mp_down_state);


%     %find outliers in tslu and tsld
%     z_tslu = tslu - nanmean(tslu);
%     z_tslu = z_tslu/nanstd(z_tslu);
%     tslu(abs(z_tslu) > 2) = nan;
%     z_tsld = tsld - nanmean(tsld);
%     z_tsld = z_tsld/nanstd(z_tsld);
%     tsld(abs(z_tsld) > 2) = nan;
%        plot(tslu,tsld,'o')
%         tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\tsld_v_tslu_' f_names{d}];
%     print('-dpng',tname);
%     close
   
 
%     up_locking_cv(d) = nanstd(tslu)/nanmean(tslu);
%     down_locking_cv(d) = nanstd(tsld)/nanmean(tsld);

        pers_fract(d) = length(find(cor_lfp_period_dur{d}>1))/length(cor_lfp_period_dur{d})

    hist_range = linspace(0,6,100);
    mp_lperiod_dur(d,:) = hist(lfp_period_dur{d},hist_range);
    cor_mp_lperiod_dur(d,:) = hist(cor_lfp_period_dur{d},hist_range);
%     stairs(hist_range,mp_lperiod_dur(d,:))
%     hold on
%     stairs(hist_range,cor_mp_lperiod_dur(d,:),'r')
%     xlim([0 5])
%     tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\lfp_fract_period_half_periods_with_cor' f_names{d}];
%     print('-dpng',tname);
%     close

% hist_range2 = linspace(0,5,200);
% hist(mp_up_dur_lfp_up_units{d},hist_range2)
% xlim([0 5])
%     tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\mp_up_dur_lfpunits_half_periods' f_names{d}];
%     print('-dpang',tname);
%     close

    clear near_lfp_ind
    clear prev_lfp_down

end

% save C:\WC_Germany\Persistent_activity\persistent_analysis\wcv_up_ctrig_mat_data wcv_up_ctrig_mat lf8_up_ctrig_mat lags Fsd