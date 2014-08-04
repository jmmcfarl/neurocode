clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data_new_method
load C:\WC_Germany\Persistent_activity\lf8_period_f_data_new_method
load C:\WC_Germany\Persistent_activity\UDS_synch_state_dur\UDS_synch_state_dur_data_new_method

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;
Fsd = 2016/dsf;
backlag = 5*Fsd;
forwardlag = 5*Fsd;
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
    
    lfp_period_dur{d} = zeros(1,length(synch_ups{d}));
    cor_lfp_period_dur{d} = zeros(1,length(synch_ups{d}));
    mean_lfp_up = zeros(1,length(synch_ups{d}));
    bad_rows = [];
    
    lfp_up_lag{d} = zeros(size(synch_ups{d}));
    
    for i = 1:length(synch_ups{d})

        %find duration of MP up state in units of LFP periods
        lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i)))-lf8_period_f{d}(up_trans{d}(synch_ups{d}(i)));

        %find first lfp down state after MP up transition
        first_in_lfp_down = find(down_trans8{d} > up_trans{d}(synch_ups{d}(i)),1,'first');


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

    end

 %for down triggered averages
     wcv_down_ctrig_mat{d} = zeros(length(synch_downs{d}),length(lags));
    lf8_down_ctrig_mat{d} = zeros(length(synch_downs{d}),length(lags));
        bad_rows = [];
        
    for i = 1:length(synch_downs{d})

        

        %calculate mup triggered averages
        if down_trans{d}(synch_downs{d}(i)) > backlag & down_trans{d}(synch_downs{d}(i)) ...
                < length(wcv_d)-forwardlag
            wcv_down_ctrig_mat{d}(i,:) = wcv_d(down_trans{d}(synch_downs{d}(i))-backlag: ...
                down_trans{d}(synch_downs{d}(i))+forwardlag);
            lf8_down_ctrig_mat{d}(i,:) = lf8_d(down_trans{d}(synch_downs{d}(i))-backlag: ...
                down_trans{d}(synch_downs{d}(i))+forwardlag);
        else
            wcv_down_ctrig_mat{d}(i,:) = nan;
            lf8_down_ctrig_mat{d}(i,:) = nan;
            bad_rows = [bad_rows i];
        end

    end

    
    
    mp_up_dur_lfp_up_units{d} = up_state_dur{d}(synch_ups{d})./mean_lfp_up;
    
    m_wcv_up_ctrig_mat(d,:) = nanmean(wcv_up_ctrig_mat{d});
    u_wcv_up_ctrig_mat(d,:) = m_wcv_up_ctrig_mat(d,:)+nanstd(wcv_up_ctrig_mat{d})/sqrt(length(synch_ups{d}));
    l_wcv_up_ctrig_mat(d,:) = m_wcv_up_ctrig_mat(d,:)-nanstd(wcv_up_ctrig_mat{d})/sqrt(length(synch_ups{d}));

    m_lf8_up_ctrig_mat(d,:) = nanmean(lf8_up_ctrig_mat{d});
    u_lf8_up_ctrig_mat(d,:) = m_lf8_up_ctrig_mat(d,:)+nanstd(lf8_up_ctrig_mat{d})/sqrt(length(synch_ups{d}));
    l_lf8_up_ctrig_mat(d,:) = m_lf8_up_ctrig_mat(d,:)-nanstd(lf8_up_ctrig_mat{d})/sqrt(length(synch_ups{d}));

    m_wcv_down_ctrig_mat(d,:) = nanmean(wcv_down_ctrig_mat{d});
    u_wcv_down_ctrig_mat(d,:) = m_wcv_down_ctrig_mat(d,:)+nanstd(wcv_down_ctrig_mat{d})/sqrt(length(synch_downs{d}));
    l_wcv_down_ctrig_mat(d,:) = m_wcv_down_ctrig_mat(d,:)-nanstd(wcv_down_ctrig_mat{d})/sqrt(length(synch_downs{d}));

    m_lf8_down_ctrig_mat(d,:) = nanmean(lf8_down_ctrig_mat{d});
    u_lf8_down_ctrig_mat(d,:) = m_lf8_down_ctrig_mat(d,:)+nanstd(lf8_down_ctrig_mat{d})/sqrt(length(synch_downs{d}));
    l_lf8_down_ctrig_mat(d,:) = m_lf8_down_ctrig_mat(d,:)-nanstd(lf8_down_ctrig_mat{d})/sqrt(length(synch_downs{d}));

    [dummy,up_order] = sort(up_state_dur{d}(synch_ups{d}));

    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    
%     subplot(2,1,1)
%     pcolor(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,wcv_up_ctrig_mat{d}(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))...
%         /length(synch_ups{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     xlim([-2 10])
%     subplot(2,1,2)
%     pcolor(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d})...
%         ,lf8_up_ctrig_mat{d}(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     xlim([-2 10])
%     tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\new_method_' f_names{d}];
%     print('-dpng',tname);
%     close


 
    pers_fract(d) = length(find(cor_lfp_period_dur{d}>1))/length(cor_lfp_period_dur{d})
    pers_fract_noshift(d) = length(find(lfp_period_dur{d}>1))/length(cor_lfp_period_dur{d})

    hist_range = linspace(0,6,100);
    mp_lperiod_dur(d,:) = hist(lfp_period_dur{d},hist_range);
    cor_mp_lperiod_dur(d,:) = hist(cor_lfp_period_dur{d},hist_range);
%     stairs(hist_range,mp_lperiod_dur(d,:))
%     hold on
%     stairs(hist_range,cor_mp_lperiod_dur(d,:),'r')
%     xlim([0 5])
%     tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\new_method_' f_names{d}];
%     print('-dpng',tname);
%     close


    clear near_lfp_ind
    clear prev_lfp_down

end

save C:\WC_Germany\Persistent_activity\persistent_analysis\wcv_up_ctrig_data_new_method ...
    wcv_up_ctrig_mat lf8_up_ctrig_mat lags Fsd *dur pers_fract* hist_range mp_up_dur_lfp_up_units ...
    *ctrig_mat