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
maxlag = 1*Fsd;
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


for d = 1:length(dir_array)
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

    %     lf5_uds = filtfilt(blow,alow,lf5);
    %     lf6_uds = filtfilt(blow,alow,lf6);
    %     lf7_uds = filtfilt(blow,alow,lf7);
    lf8_uds = filtfilt(blow,alow,lf8);
    wcv = filtfilt(blow,alow,wcv_minus_spike);

    %    lf2_theta = downsample(lf2_theta,dsf);
    %    lf2_ripple = downsample(lf2_ripple,dsf);
        lf2_uds = downsample(lf2_uds,dsf);
    lf3_uds = downsample(lf3_uds,dsf);
    %     lf5_uds = downsample(lf5_uds,dsf);
    %     lf6_uds = downsample(lf6_uds,dsf);
    %     lf7_uds = downsample(lf7_uds,dsf);
    lf8_uds = downsample(lf8_uds,dsf);
    wcv = downsample(wcv,dsf);


    spike_id = round(spkid/dsf);

    %    lf2_theta = zscore(lf2_theta);
    %    lf2_ripple = zscore(lf2_ripple);
        lf2_uds = zscore(lf2_uds);
    %     lf6_uds = zscore(lf6_uds);
    %     lf7_uds = zscore(lf7_uds);
    lf8_uds = zscore(lf8_uds);
    %     lf5_uds = zscore(lf5_uds);
    lf3_uds = zscore(lf3_uds);
    wcv = zscore(wcv);

    lf3_uds = lf3_uds*flip_sign(d);
    
    t_axis = (1:length(wcv))/Fsd;

    up_times = zeros(size(t_axis));
    
    for i = 1:length(synch_ups{d})
        
        up_times(up_trans{d}(synch_ups{d}(i)):down_trans{d}(synch_ups{d}(i))) = 1;
        
    end
    up_times8 = zeros(size(t_axis));
    for i = 1:length(synch_ups8{d})
       up_times8(up_trans8{d}(synch_ups8{d}(i)):down_trans8{d}(synch_ups8{d}(i))) = 1; 
    end
    up_times = logical(up_times);
    up_times8 = logical(up_times8);
    %calculate mp state conditional distribution of lf8 and lf3
    lf3_mp_up = lf3_uds(up_times);
    lf3_mp_down = lf3_uds(~up_times);
        lf2_mp_up = lf2_uds(up_times);
    lf2_mp_down = lf2_uds(~up_times);

    lf8_mp_up = lf8_uds(up_times);
    lf8_mp_down=  lf8_uds(~up_times);
    
    lf8_up_cond_mp_up(d) = length(find(up_times8 & up_times))/length(find(up_times));
    lf8_down_cond_mp_up(d) = 1-lf8_up_cond_mp_up(d);
    lf8_up_cond_mp_down(d) = length(find(up_times8 & ~up_times))/length(find(~up_times));
    lf8_down_cond_mp_down(d) = 1-lf8_up_cond_mp_down(d);
%     [lf3_mp_up_dist(d,:),lf3_mp_up_grid(d,:)] = gpkde(lf3_mp_up,-1);
%     [lf3_mp_down_dist(d,:),lf3_mp_down_grid(d,:)] = gpkde(lf3_mp_down,-1);
%         [lf2_mp_up_dist(d,:),lf2_mp_up_grid(d,:)] = gpkde(lf2_mp_up,-1);
%     [lf2_mp_down_dist(d,:),lf2_mp_down_grid(d,:)] = gpkde(lf2_mp_down,-1);
% 
%     [lf8_mp_up_dist(d,:),lf8_mp_up_grid(d,:)] = gpkde(lf8_mp_up,-1);
%     [lf8_mp_down_dist(d,:),lf8_mp_down_grid(d,:)] = gpkde(lf8_mp_down,-1);
% 
%     [lf3_ov_dist(d,:),lf3_ov_grid(d,:)] = gpkde(lf3_uds,-1);
%         [lf2_ov_dist(d,:),lf2_ov_grid(d,:)] = gpkde(lf3_uds,-1);
% 
%        [lf8_ov_dist(d,:),lf8_ov_grid(d,:)] = gpkde(lf8_uds,-1);
 
%     subplot(3,1,1)
%     plot(lf3_mp_up_grid(d,:),lf3_mp_up_dist(d,:))
%     hold on
%     plot(lf3_mp_down_grid(d,:),lf3_mp_down_dist(d,:),'r')
%         plot(lf3_ov_grid(d,:),lf3_ov_dist(d,:),'k')
%         xlim([-3 3])
%         if flip_sign(d) > 0
%             title('Flip LF3')
%         else
% title('Lf3')
%         end
% legend('MP up','MP down','overall')
%     subplot(3,1,2)
%     plot(lf2_mp_up_grid(d,:),lf2_mp_up_dist(d,:))
%     hold on
%     plot(lf2_mp_down_grid(d,:),lf2_mp_down_dist(d,:),'r')
%         plot(lf2_ov_grid(d,:),lf2_ov_dist(d,:),'k')
%                 xlim([-3 3])
%         if flip_sign(d) > 0
%             title('Flip LF2')
%         else
% title('Lf2')
%         end
% legend('MP up','MP down','overall')
% 
%     subplot(3,1,3)
%     plot(lf8_mp_up_grid(d,:),lf8_mp_up_dist(d,:))
%     hold on
%     plot(lf8_mp_down_grid(d,:),lf8_mp_down_dist(d,:),'r')
%            plot(lf8_ov_grid(d,:),lf8_ov_dist(d,:),'k')
%                    xlim([-3 3])
%         if flip_sign(d) > 0
%             title('Flip LF8')
%         else
% title('Lf8')
%         end
% legend('MP up','MP down','overall')
% t_names = ['C:\WC_Germany\Persistent_activity\lf3_lf8_cond_dist\' f_names{d}];
% print('-dpng',t_names)
% close
% 
% %now calculate spike triggered average
% spike_trig_mat_3 = zeros(length(spike_id),length(lags));
% spike_trig_mat_8 = zeros(length(spike_id),length(lags));
%     for i = 1:length(spike_id)
%         if spike_id(i) > maxlag & spike_id(i) < length(lf3_uds)-maxlag
%         spike_trig_mat_3(i,:) = lf3_uds(spike_id(i) - maxlag:spike_id(i)+maxlag);
%         spike_trig_mat_8(i,:) = lf8_uds(spike_id(i) - maxlag:spike_id(i)+maxlag);
%         else
%             spike_trig_mat_3(i,:) = nan;
%             spike_trig_mat_8(i,:) = nan;
%         end
%         
%     end
%     
%     spike_trig_lf3(d,:) = nanmean(spike_trig_mat_3);
%     spike_trig_lf8(d,:) = nanmean(spike_trig_mat_8);
%     
%     figure
%     plot(lags/Fsd,spike_trig_lf3(d,:))
%     hold on
%     plot(lags/Fsd,spike_trig_lf8(d,:),'r')
%     legend('Lf3','Lf8')
%     grid
%     if flip_sign(d) > 0
%         title('Flip')
%     else
%         title('no flip')
%     end
% t_names = ['C:\WC_Germany\Persistent_activity\lf3_lf8_cond_dist\spike_trig_' f_names{d}];
% print('-dpng',t_names)
% close
    
    
end