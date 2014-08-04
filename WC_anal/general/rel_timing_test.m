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

for d = 1:length(dir_array)
    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8_f = filtfilt(b,a,lf8);
    wcv_f = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8_f,dsf);
    wcv_d = downsample(wcv_f,dsf);

    datalen = length(lf8_d);
    

    %initialize
    lfp_period_dur{d} = zeros(1,length(synch_ups{d}));
    cor_lfp_period_dur{d} = zeros(1,length(synch_ups{d}));
    bad_rows = [];
    
    lfp_up_lag{d} = zeros(size(synch_ups{d}));
    
    tslu{d} = nan(size(synch_ups{d}));
    ttnu{d} = nan(size(synch_ups{d}));
    tsld{d} = nan(size(synch_ups{d}));
    ttnd{d} = nan(size(synch_ups{d}));
    ttnextup{d} = nan(size(synch_ups{d}));    
    
    
   last_up_dur{d} = nan(size(synch_ups{d}));
   last_down_dur{d} = nan(size(synch_ups{d}));
    
    for i = 1:length(synch_ups{d})

        cur_up_trans = up_trans{d}(synch_ups{d}(i));
        cur_down_trans = down_trans{d}(synch_ups{d}(i));
        
        [dummy,nearest_lfp_up] = min(abs(up_trans8{d}-cur_up_trans));
        if ~isempty(nearest_lfp_up)
        ttnu{d}(i) = (cur_up_trans-up_trans8{d}(nearest_lfp_up))/Fsd;
        end
        
        previous_lfp_up = find(up_trans8{d} < cur_up_trans,1,'last');
        if ~isempty(previous_lfp_up)
        tslu{d}(i) = (cur_up_trans-up_trans8{d}(previous_lfp_up))/Fsd;
        last_up_dur{d}(i) = up_state_dur8{d}(previous_lfp_up);
        end
        
        [dummy,nearest_lfp_down] = min(abs(down_trans8{d}-cur_down_trans));
        if ~isempty(nearest_lfp_down)
        ttnd{d}(i) = (cur_down_trans-down_trans8{d}(nearest_lfp_down))/Fsd;
        end
        
        previous_lfp_down = find(down_trans8{d} < cur_down_trans,1,'last');
        if ~isempty(previous_lfp_down)
        tsld{d}(i) = (cur_down_trans - down_trans8{d}(previous_lfp_down))/Fsd;
        if length(down_state_dur8{d}) >= previous_lfp_down
        last_down_dur{d}(i) = down_state_dur8{d}(previous_lfp_down);
        else
            last_down_dur{d}(i) = nan;
            tsld{d}(i) = nan;
        end
        end
        
        next_up = find(up_trans8{d} > cur_down_trans,1,'first');
        if ~isempty(next_up)
            ttnextup{d}(i) = (up_trans8{d}(next_up)-cur_down_trans)/Fsd;
        end
        
        %find duration of MP up state in units of LFP periods
        lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i)))-lf8_period_f{d}(up_trans{d}(synch_ups{d}(i)));

        %acquire MP up related lfp transitions
        [dummy,near_lfp_up] = min(abs(up_trans{d}(synch_ups{d}(i)) - up_trans8{d}(synch_ups8{d})));
        lfp_up_lag{d}(i) = up_trans{d}(synch_ups{d}(i))-up_trans8{d}(synch_ups8{d}(near_lfp_up));
        
        prev_lfp_up = find(up_trans8{d}<up_trans{d}(synch_ups{d}(i)),1,'last');
        prev_lfp_down = find(down_trans8{d}<down_trans{d}(synch_ups{d}(i)),1,'last');
        next_lfp_up = find(up_trans8{d}>down_trans{d}(synch_ups{d}(i)),1,'first');
        if ~isempty(prev_lfp_up)
            up_lag = up_trans{d}(synch_ups{d}(i))-up_trans8{d}(prev_lfp_up);
            cor_lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i))-up_lag)-lf8_period_f{d}(up_trans{d}(synch_ups{d}(i))-up_lag);
        end
    end

    used_pts = find(abs(zscore(cor_lfp_period_dur{d})) < 3);
    round_cor_lfp_mean(d) = nanmean(floor(cor_lfp_period_dur{d}(used_pts)));
        used_pts = find(abs(zscore(lfp_period_dur{d})) < 3);
    round_lfp_mean(d) = nanmean(floor(lfp_period_dur{d}(used_pts)));
        round_cor_lfp_med(d) = nanmedian(floor(cor_lfp_period_dur{d}));
    round_lfp_med(d) = nanmedian(floor(lfp_period_dur{d}));

%     hist_range = linspace(0,6,100);
%     mp_lperiod_dur(d,:) = hist(lfp_period_dur{d},hist_range);
%     cor_mp_lperiod_dur(d,:) = hist(cor_lfp_period_dur{d},hist_range);
%     stairs(hist_range,mp_lperiod_dur(d,:))
%     hold on
%     stairs(hist_range,cor_mp_lperiod_dur(d,:),'r')
%     xlim([0 5])
%     tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\raw_compare_' f_names{d}];
%     print('-dpng',tname);
%     close

%     hist_range = linspace(0,5,200);
%     ttnu_binned(d,:) = hist(ttnu{d},hist_range);
%     ttnd_binned(d,:) = hist(ttnd{d},hist_range);
%     tslu_binned(d,:) = hist(tslu{d},hist_range);
%     tsld_binned(d,:) = hist(tsld{d},hist_range);
% %     stairs(hist_range,ttnu_binned(d,:))
%     stairs(hist_range,tslu_binned(d,:),'r')
%         hold on
% 
% %     stairs(hist_range,ttnd_binned(d,:),'k')
%     stairs(hist_range,tsld_binned(d,:),'c')
%     xlim([0 4.9])
%     pause
%     clf
% %     tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\raw_compare_' f_names{d}];
% %     print('-dpng',tname);
% %     close
% 
med_tslu(d) = nanmedian(tslu{d});
med_tsld(d) = nanmedian(tsld{d});
mean_tslu(d) = nanmean(tslu{d});
mean_tsld(d) = nanmean(tsld{d});
cv_tslu(d) = nanstd(tslu{d})/mean_tslu(d);
cv_tsld(d) = nanstd(tsld{d})/mean_tsld(d);

% hist_range = linspace(0,5,200);
% next_up_binned(d,:) = hist(ttnextup{d},hist_range);
% prev_down_binned(d,:) = hist(tsld{d},hist_range);
% stairs(hist_range,next_up_binned(d,:))
% hold on
% stairs(hist_range,prev_down_binned(d,:),'r')
% pause
% clf
synch_up_dur = up_state_dur{d}(synch_ups{d});
good_pts = find(~isnan(synch_up_dur) & ~isnan(tslu{d}) & ~isnan(tsld{d}));

no_pers = find(lfp_period_dur{d}(good_pts) < 1);
one_pers = find(lfp_period_dur{d}(good_pts) > 1 & lfp_period_dur{d}(good_pts) < 2);
two_pers = find(lfp_period_dur{d}(good_pts) > 2 & lfp_period_dur{d}(good_pts) < 3);
many_pers = find(lfp_period_dur{d}(good_pts) > 3);
pers = find(lfp_period_dur{d}(good_pts) > 1);


tsld_diff{d} = tsld{d}-tslu{d};
mean_diff(d) = nanmean(tsld_diff{d});
med_diff(d) = nanmedian(tsld_diff{d});
% [r1(d),p1(d)] = corr(synch_up_dur(good_pts)',tsld_diff') 
% [r2(d),p2(d)] = corr(synch_up_dur(good_pts)',tsld_diff','type','spearman')

% hist_range = linspace(0,3,100);
% np_binned = hist(tsld{d}(good_pts(no_pers)),hist_range);
% p_binned = hist(tsld{d}(good_pts(pers)),hist_range);
% [sigval(d)] = ranksum(tsld{d}(good_pts(no_pers)),tsld{d}(good_pts(pers)));
% stairs(hist_range,np_binned)
% hold on
% stairs(hist_range,p_binned,'r')
% title(['sig = ' num2str(sigval)])
% pause
% clf

% synch_up_dur = up_state_dur{d}(synch_ups{d});
% good_cases = find(tsld{d} < 3);
% scatter(tsld{d}(good_pts(no_pers)),synch_up_dur(good_pts(no_pers)),'.')
% hold on
% scatter(tsld{d}(good_pts(one_pers)),synch_up_dur(good_pts(one_pers)),'r.')
% scatter(tsld{d}(good_pts(two_pers)),synch_up_dur(good_pts(two_pers)),'g.')
% scatter(tsld{d}(good_pts(many_pers)),synch_up_dur(good_pts(many_pers)),'c.')

% scatter(tsld{d}(good_pts(no_pers)),tslu{d}(good_pts(no_pers)),'.')
% hold on
% scatter(tsld{d}(good_pts(one_pers)),tslu{d}(good_pts(one_pers)),'r.')
% scatter(tsld{d}(good_pts(two_pers)),tslu{d}(good_pts(two_pers)),'g.')
% scatter(tsld{d}(good_pts(many_pers)),tslu{d}(good_pts(many_pers)),'c.')

% scatter(tsld_diff(no_pers),synch_up_dur(good_pts(no_pers)),'.')
% hold on
% scatter(tsld_diff(one_pers),synch_up_dur(good_pts(one_pers)),'r.')
% scatter(tsld_diff(two_pers),synch_up_dur(good_pts(two_pers)),'g.')
% scatter(tsld_diff(many_pers),synch_up_dur(good_pts(many_pers)),'c.')

% pause
% clf

    clear near_lfp_ind
    clear prev_lfp_down

end

% save
% C:\WC_Germany\Persistent_activity\persistent_analysis\wcv_up_ctrig_mat_data wcv_up_ctrig_mat lf8_up_ctrig_mat lags Fsd