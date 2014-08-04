clear all
close all
load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

Fs = 2016;
niqf = Fs/2;
dsf = 8;
Fsd = Fs/dsf;

lcf = 0.05/niqf;
hcf = 40/niqf;
[blow,alow] = butter(2,[lcf hcf]);
maxlag = 10*Fsd;
lags = -maxlag:maxlag;
% lcf = 2/niqf;
% hcf = 10/niqf;
% [bmid,amid] = butter(2,[lcf hcf]);
%
% lcf = 100/niqf;
% hcf = 250/niqf;
% [bhigh,ahigh] = butter(2,[lcf hcf]);
flip_sign = [1 1 1 1 1 1 -1 1 -1 1 -1 -1 -1 1 1 -1 -1];

win = 20;


% for d = 1:length(dir_array)
d=14
    cd(dir_array{d})
    pwd

    load used_data lf2 lf3 lf5 lf6 lf7 lf8 wcv_minus_spike
    if exist('spike_time.mat') > 0
        load spike_time
    else
        load spike_time_br
    end
    %    lf2_theta = filtfilt(bmid,amid,lf2);
    %    lf2_ripple = filtfilt(bhigh,ahigh,lf2);
    lf2_uds = filtfilt(blow,alow,lf2);

    lf3_uds = filtfilt(blow,alow,lf3);

    lf5_uds = filtfilt(blow,alow,lf5);
    %     lf6_uds = filtfilt(blow,alow,lf6);
    lf7_uds = filtfilt(blow,alow,lf7);
    lf8_uds = filtfilt(blow,alow,lf8);
    wcv = filtfilt(blow,alow,wcv_minus_spike);

    %    lf2_theta = downsample(lf2_theta,dsf);
    %    lf2_ripple = downsample(lf2_ripple,dsf);
    lf2_uds = downsample(lf2_uds,dsf);
    lf3_uds = downsample(lf3_uds,dsf);
    lf5_uds = downsample(lf5_uds,dsf);
    %     lf6_uds = downsample(lf6_uds,dsf);
    lf7_uds = downsample(lf7_uds,dsf);
    lf8_uds = downsample(lf8_uds,dsf);
    wcv = downsample(wcv,dsf);

    %
    lf3_uds = lf3_uds+lf5_uds;
    
    
    spike_id = round(spkid/dsf);

    %    lf2_theta = zscore(lf2_theta);
    %    lf2_ripple = zscore(lf2_ripple);
    lf2_uds = zscore(lf2_uds);
    %     lf6_uds = zscore(lf6_uds);
    lf7_uds = zscore(lf7_uds);
    lf8_uds = zscore(lf8_uds);
    lf5_uds = zscore(lf5_uds);
    lf3_uds = zscore(lf3_uds);
    wcv = zscore(wcv);

    t_axis = (1:length(wcv))/Fsd;

    %% up trig avges

    lf2_lf8_up_tr_mat = zeros(length(synch_ups8{d}),length(lags));
    lf3_lf8_up_tr_mat = zeros(length(synch_ups8{d}),length(lags));
    lf5_lf8_up_tr_mat = zeros(length(synch_ups8{d}),length(lags));
    %     lf6_lf8_up_tr_mat = zeros(length(synch_ups{d}),length(lags));
    lf7_lf8_up_tr_mat = zeros(length(synch_ups8{d}),length(lags));

    lf8_lf8_up_tr_mat = zeros(length(synch_ups8{d}),length(lags));
    mp_up_tr_mat = zeros(length(synch_ups8{d}),length(lags));

    for i = 1:length(synch_ups8{d})

        if up_trans8{d}(synch_ups8{d}(i)) > maxlag & up_trans8{d}(synch_ups8{d}(i)) ...
                < length(wcv) - maxlag
            lf2_lf8_up_tr_mat(i,:) = lf2_uds(up_trans8{d}(synch_ups8{d}(i))-maxlag:...
                up_trans8{d}(synch_ups8{d}(i))+maxlag);
            lf3_lf8_up_tr_mat(i,:) = lf3_uds(up_trans8{d}(synch_ups8{d}(i))-maxlag:...
                up_trans8{d}(synch_ups8{d}(i))+maxlag);
            lf5_lf8_up_tr_mat(i,:) = lf5_uds(up_trans8{d}(synch_ups8{d}(i))-maxlag:...
                up_trans8{d}(synch_ups8{d}(i))+maxlag);
            %                         lf6_lf8_up_tr_mat(i,:) = lf6_uds(up_trans{d}(synch_ups{d}(i))-maxlag:...
            %                 up_trans{d}(synch_ups{d}(i))+maxlag);
            lf7_lf8_up_tr_mat(i,:) = lf7_uds(up_trans8{d}(synch_ups8{d}(i))-maxlag:...
                up_trans8{d}(synch_ups8{d}(i))+maxlag);
            lf8_lf8_up_tr_mat(i,:) = lf8_uds(up_trans8{d}(synch_ups8{d}(i))-maxlag:...
                up_trans8{d}(synch_ups8{d}(i))+maxlag);
            mp_up_tr_mat(i,:) = wcv(up_trans8{d}(synch_ups8{d}(i))-maxlag:...
                up_trans8{d}(synch_ups8{d}(i))+maxlag);

        else
            lf2_lf8_up_tr_mat(i,:) = nan;
            lf3_lf8_up_tr_mat(i,:) = nan;
            lf5_lf8_up_tr_mat(i,:) = nan;
            %                         lf6_lf8_up_tr_mat(i,:) = nan;
            lf7_lf8_up_tr_mat(i,:) = nan;
            lf8_lf8_up_tr_mat(i,:) = nan;
            mp_up_tr_mat(i,:) = nan;
        end

    end

    av_lf2_lf8_up_tr_mat(d,:) = nanmean(lf2_lf8_up_tr_mat);
    av_lf3_lf8_up_tr_mat(d,:) = nanmean(lf3_lf8_up_tr_mat);
    av_lf5_lf8_up_tr_mat(d,:) = nanmean(lf5_lf8_up_tr_mat);
    %         av_lf6_lf8_up_tr_mat(d,:) = nanmean(lf6_lf8_up_tr_mat);
    av_lf7_lf8_up_tr_mat(d,:) = nanmean(lf7_lf8_up_tr_mat);
    av_lf8_lf8_up_tr_mat(d,:) = nanmean(lf8_lf8_up_tr_mat);
    av_mp_up_tr_mat(d,:) = nanmean(mp_up_tr_mat);

    %% now for down trig avges

    lf2_lf8_down_tr_mat = zeros(length(synch_ups8{d}),length(lags));
    lf3_lf8_down_tr_mat = zeros(length(synch_ups8{d}),length(lags));
    lf5_lf8_down_tr_mat = zeros(length(synch_ups8{d}),length(lags));
    %         lf6_lf8_down_tr_mat = zeros(length(synch_ups{d}),length(lags));
    lf7_lf8_down_tr_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_lf8_down_tr_mat = zeros(length(synch_ups8{d}),length(lags));
    mp_down_tr_mat = zeros(length(synch_ups8{d}),length(lags));

    for i = 1:length(synch_ups8{d})

        if down_trans8{d}(synch_ups8{d}(i)) > maxlag & down_trans8{d}(synch_ups8{d}(i)) ...
                < length(wcv) - maxlag
            lf2_lf8_down_tr_mat(i,:) = lf2_uds(down_trans8{d}(synch_ups8{d}(i))-maxlag:...
                down_trans8{d}(synch_ups8{d}(i))+maxlag);
            lf3_lf8_down_tr_mat(i,:) = lf3_uds(down_trans8{d}(synch_ups8{d}(i))-maxlag:...
                down_trans8{d}(synch_ups8{d}(i))+maxlag);
            lf5_lf8_down_tr_mat(i,:) = lf5_uds(down_trans8{d}(synch_ups8{d}(i))-maxlag:...
                down_trans8{d}(synch_ups8{d}(i))+maxlag);
            %                         lf6_lf8_down_tr_mat(i,:) = lf6_uds(down_trans{d}(synch_ups{d}(i))-maxlag:...
            %                 down_trans{d}(synch_ups{d}(i))+maxlag);
            lf7_lf8_down_tr_mat(i,:) = lf7_uds(down_trans8{d}(synch_ups8{d}(i))-maxlag:...
                down_trans8{d}(synch_ups8{d}(i))+maxlag);
            lf8_lf8_down_tr_mat(i,:) = lf8_uds(down_trans8{d}(synch_ups8{d}(i))-maxlag:...
                down_trans8{d}(synch_ups8{d}(i))+maxlag);
            mp_down_tr_mat(i,:) = wcv(down_trans8{d}(synch_ups8{d}(i))-maxlag:...
                down_trans8{d}(synch_ups8{d}(i))+maxlag);

        else
            lf2_lf8_down_tr_mat(i,:) = nan;
            lf3_lf8_down_tr_mat(i,:) = nan;
            lf5_lf8_down_tr_mat(i,:) = nan;
            %                         lf6_lf8_down_tr_mat(i,:) = nan;
            lf7_lf8_down_tr_mat(i,:) = nan;
            lf8_lf8_down_tr_mat(i,:) = nan;
            mp_down_tr_mat(i,:) = nan;
        end

    end

    av_lf2_lf8_down_tr_mat(d,:) = nanmean(lf2_lf8_down_tr_mat);
    av_lf3_lf8_down_tr_mat(d,:) = nanmean(lf3_lf8_down_tr_mat);
    av_lf5_lf8_down_tr_mat(d,:) = nanmean(lf5_lf8_down_tr_mat);
    %         av_lf6_lf8_down_tr_mat(d,:) = nanmean(lf6_lf8_down_tr_mat);
    av_lf7_lf8_down_tr_mat(d,:) = nanmean(lf7_lf8_down_tr_mat);
    av_lf8_lf8_down_tr_mat(d,:) = nanmean(lf8_lf8_down_tr_mat);
    av_mp_down_tr_mat(d,:) = nanmean(mp_down_tr_mat);

%     %% create plots
%     
        close all
        figure
        plot(lags/Fsd,av_mp_up_tr_mat(d,:),'linewidth',2)
        hold on
            plot(lags/Fsd,av_lf2_lf8_up_tr_mat(d,:),'c','linewidth',2)
    plot(lags/Fsd,av_lf3_lf8_up_tr_mat(d,:),'k','linewidth',2)
            plot(lags/Fsd,av_lf5_lf8_up_tr_mat(d,:),'g','linewidth',2)
    %     plot(lags/Fsd,av_lf6_lf8_up_tr_mat(d,:),'y','linewidth',2)
        plot(lags/Fsd,av_lf7_lf8_up_tr_mat(d,:),'m','linewidth',2)
    plot(lags/Fsd,av_lf8_lf8_up_tr_mat(d,:),'r','linewidth',2)
        legend('MP','LF2','LF3','LF5','LF7','LF8')
        xlim([-5 5])
%         t_names = ['C:\WC_Germany\Persistent_activity\lf8_trig_all_lfp\up_' f_names{d}];
%         print('-dpng',t_names)
%         close all
%     
        figure
        plot(lags/Fsd,av_mp_down_tr_mat(d,:),'linewidth',2)
        hold on
            plot(lags/Fsd,av_lf2_lf8_down_tr_mat(d,:),'c','linewidth',2)
    plot(lags/Fsd,av_lf3_lf8_down_tr_mat(d,:),'k','linewidth',2)
            plot(lags/Fsd,av_lf5_lf8_down_tr_mat(d,:),'g','linewidth',2)
    %         plot(lags/Fsd,av_lf6_lf8_down_tr_mat(d,:),'y','linewidth',2)
            plot(lags/Fsd,av_lf7_lf8_down_tr_mat(d,:),'m','linewidth',2)
    plot(lags/Fsd,av_lf8_lf8_down_tr_mat(d,:),'r','linewidth',2)
        xlim([-5 5])
%         legend('MP','LF2','LF3','LF5','LF7','LF8')
%         t_names = ['C:\WC_Germany\Persistent_activity\lf8_trig_all_lfp\down_' f_names{d}];
%         print('-dpng',t_names)
%         close all
%     
%     
%         close all
%         figure
%             plot(lags/Fsd,av_mp_up_tr_mat(d,:),'linewidth',2)
%             hold on
%             plot(lags/Fsd,av_lf2_lf8_up_tr_mat(d,:),'c','linewidth',2)
%     plot(lags/Fsd,flip_sign(d)*av_lf3_lf8_up_tr_mat(d,:),'k','linewidth',2)
%             plot(lags/Fsd,av_lf5_lf8_up_tr_mat(d,:),'g','linewidth',2)
%     %     plot(lags/Fsd,av_lf6_lf8_up_tr_mat(d,:),'y','linewidth',2)
%         plot(lags/Fsd,av_lf7_lf8_up_tr_mat(d,:),'m','linewidth',2)
%     plot(lags/Fsd,av_lf8_lf8_up_tr_mat(d,:),'r','linewidth',2)
%         legend('MP','LF2','LF3','LF5','LF7','LF8')
%         xlim([-5 5])
%     
%         t_names = ['C:\WC_Germany\Persistent_activity\lf8_trig_all_lfp\flip_up_' f_names{d}];
%         print('-dpng',t_names)
%         close all
%     
%         figure
%                 plot(lags/Fsd,av_mp_down_tr_mat(d,:),'linewidth',2)
%         hold on
%     plot(lags/Fsd,av_lf2_lf8_down_tr_mat(d,:),'c','linewidth',2)
%     plot(lags/Fsd,flip_sign(d)*av_lf3_lf8_down_tr_mat(d,:),'k','linewidth',2)
%             plot(lags/Fsd,av_lf5_lf8_down_tr_mat(d,:),'g','linewidth',2)
%     %         plot(lags/Fsd,av_lf6_lf8_down_tr_mat(d,:),'y','linewidth',2)
%             plot(lags/Fsd,av_lf7_lf8_down_tr_mat(d,:),'m','linewidth',2)
%     plot(lags/Fsd,av_lf8_lf8_down_tr_mat(d,:),'r','linewidth',2)
%         xlim([-5 5])
%         legend('MP','LF2','LF3','LF5','LF7','LF8')
%     
%         t_names = ['C:\WC_Germany\Persistent_activity\lf8_trig_all_lfp\flip_down_' f_names{d}];
%         print('-dpng',t_names)
%         close all
    
    
%% create conditional LFP matrix
    lf8_cond_mat{d} = zeros(length(lags),300);
    lf3_cond_mat{d} = zeros(length(lags),300);

    bad_rows = [];
    for i = 1:size(lf8_lf8_up_tr_mat,1)
        if max(isnan(lf8_lf8_up_tr_mat(i,:))) > 0
            bad_rows = [bad_rows i];
        end
    end
    lf8_lf8_up_tr_mat(bad_rows,:) = [];
    lf3_lf8_up_tr_mat(bad_rows,:) = [];
    for i = 1:length(lags)

        [lf8_cond_mat{d}(i,:),cond_dist] = gpkde(squeeze(lf8_lf8_up_tr_mat(:,i)),-3,[-4,4,300]);

    end
    for i = 1:length(lags)

        [lf3_cond_mat{d}(i,:),cond_dist] = gpkde(squeeze(lf3_lf8_up_tr_mat(:,i)),-3,[-4,4,300]);

    end
    figure
subplot(2,1,1)
    pcolor(lags/Fsd,cond_dist,lf8_cond_mat{d}');shading flat
    hold on
    plot(lags/Fsd,av_lf8_lf8_up_tr_mat,'w')
    xlim([-5 5])
    ylim([-3 3])
    line([0 0],[-4 4],'Color','w','linewidth',2)
subplot(2,1,2)
    pcolor(lags/Fsd,cond_dist,lf3_cond_mat{d}');shading flat
    hold on
    plot(lags/Fsd,av_lf3_lf8_up_tr_mat,'w')
    xlim([-5 5])
    ylim([-3 3])
        line([0 0],[-4 4],'Color','w','linewidth',2)

% end