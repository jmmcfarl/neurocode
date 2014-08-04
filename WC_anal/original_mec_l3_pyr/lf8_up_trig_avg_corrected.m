clear all
close all
load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

Fs = 2016;
niqf = Fs/2;
dsf = 8;
Fsd = Fs/dsf;

lcf = 10/niqf;
hcf = 100/niqf;
[blow,alow] = butter(2,[lcf hcf]);
maxlag = 1*Fsd;
lags = -maxlag:maxlag;

backlag = 3*Fsd;
forwardlag = 8*Fsd
lags = -backlag:forwardlag;

% lcf = 2/niqf;
% hcf = 10/niqf;
% [bmid,amid] = butter(2,[lcf hcf]);
%
% lcf = 100/niqf;
% hcf = 250/niqf;
% [bhigh,ahigh] = butter(2,[lcf hcf]);
flip_sign = [1 1 1 1 1 1 -1 1 -1 1 -1 -1 -1 1 1 -1 -1];

win = 20;
% scale_fact = [.05 0 -0.01 0.15 0.1 0.15 0.52 0.28 0.53 0.08 0.75 0.81 0.62 0.08 0.14 0.7 0.35];

for d = 1:length(dir_array)
    

    cd(dir_array{d})
    pwd

    load used_data lf2 lf3 lf5 lf6 lf7 lf8 wcv_minus_spike
    wcv_f = filtfilt(blow,alow,wcv_minus_spike);
    hipp_lf = lf3 - lf2;
    hipp_lf = filtfilt(blow,alow,hipp_lf);
%     lf3_uds = filtfilt(blow,alow,lf3);
%     lf2_uds = filtfilt(blow,alow,lf2);
%     lf5_uds = filtfilt(blow,alow,lf5);
%     lf8_uds = filtfilt(blow,alow,lf8);
    wcv_d = downsample(wcv_f,dsf);
    hipp_lf = downsample(hipp_lf,dsf);
%     lf3_uds = downsample(lf3_uds,dsf);
%     lf8_uds = downsample(lf8_uds,dsf);
%     lf5_uds = downsample(lf5_uds,dsf);
%     lf2_uds = downsample(lf2_uds,dsf);
    
%     lf8_uds = zscore(lf8_uds);
%     lf3_uds = zscore(lf3_uds);
%     lf5_uds = zscore(lf5_uds);
    
%     a = corrcoef(lf3_uds,lf5_uds);
%     scale_fact = a(2,1);
%     lf3_t = lf3_uds-scale_fact*lf5_uds;
%     t_axis = (1:length(lf3_uds))/Fsd;
hipp_lf = jmm_smooth_1d(hipp_lf.^2,20);
hipp_lf = zscore(hipp_lf);
    wcv_d = zscore(wcv_d);
    
    clear *_f clear lf8 lf2 lf3 lf5 lf7
    
    
    wcv_up_ctrig_mat = zeros(length(up_trans{d}),length(lags));
%     lf8_up_ctrig_mat = wcv_up_ctrig_mat;
%     lf2_up_ctrig_mat = wcv_up_ctrig_mat;
%     lf3_up_ctrig_mat = wcv_up_ctrig_mat;
%         lf3t_up_ctrig_mat = wcv_up_ctrig_mat;
    hipp_up_ctrig_mat = wcv_up_ctrig_mat;
%     lf5_up_ctrig_mat = wcv_up_ctrig_mat;
%     lf7_up_ctrig_mat = wcv_up_ctrig_mat;
    

    for i = 1:length(up_trans{d})

        if up_trans{d}(i) > backlag & up_trans{d}(i) < length(wcv_d)-forwardlag
%             wcv_up_ctrig_mat(i,:) = wcv_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
%             lf8_up_ctrig_mat(i,:) = lf8_uds(up_trans8{d}(i)-backlag:up_trans8{d}(i)+forwardlag);
%             lf2_up_ctrig_mat(i,:) = lf2_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
%             lf3_up_ctrig_mat(i,:) = lf3_uds(up_trans8{d}(i)-backlag:up_trans8{d}(i)+forwardlag);
%             lf3t_up_ctrig_mat(i,:) = lf3_t(up_trans8{d}(i)-backlag:up_trans8{d}(i)+forwardlag);
              hipp_up_ctrig_mat(i,:) = hipp_lf(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
%             lf5_up_ctrig_mat(i,:) = lf5_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
%             lf7_up_ctrig_mat(i,:) = lf7_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);


        else
%             wcv_up_ctrig_mat(i,:) = nan;
%             lf8_up_ctrig_mat(i,:) = nan;
%             lf2_up_ctrig_mat(i,:) = nan;
%             lf3_up_ctrig_mat(i,:) = nan;
%                         lf3t_up_ctrig_mat(i,:) = nan;
                hipp_up_ctrig_mat(i,:) = nan;
%             lf5_up_ctrig_mat(i,:) = nan;
%             lf7_up_ctrig_mat(i,:) = nan;

        end

    end

%   
% 
    [dummy,up_order] = sort(up_state_dur{d});
% 


    d
    %
    %   figure
    %       subplot(2,1,1)
    %     pcolor(dlags/Fsd,1:length(down_trans{d}),wcv_down_ctrig_mat(up_order,:));shading flat
    %     subplot(2,1,2)
%     close all
        pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d}),hipp_up_ctrig_mat(up_order,:));shading flat
        hold on
    plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
    line([0 0],[0 1],'Color','k')
    caxis([0 3]);colorbar;
        tname = ['C:\WC_Germany\JMM_Analysis_Pyr\newest_mp_up_trig_lf8_lf3_lf2\high_freq_hipp_' f_names{d}];
        print('-dpng',tname);
        close
% %     %
%         pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d}),lf2_up_ctrig_mat(up_order,:));shading flat
%         hold on
%                 caxis([-1 3]);colorbar
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%         tname = ['C:\WC_Germany\JMM_Analysis_Pyr\newest_mp_up_trig_lf8_lf3_lf2\lf2_' f_names{d}];
%         print('-dpng',tname);
%         close
        
%         pcolor(lags/Fsd,(1:length(up_trans8{d}))/length(up_trans8{d}),lf3_up_ctrig_mat(up_order,:));shading flat
%         hold on
%                 caxis([-1 3]);colorbar
%     plot(up_state_dur8{d}(up_order),(1:length(up_trans8{d}))/length(up_trans8{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%         tname = ['C:\WC_Germany\JMM_Analysis_Pyr\newest_lf8_up_trig_lf8_lf3_lf2\high_freq_lf3_' f_names{d}];
%         print('-dpng',tname);
%         close
% 
%         pcolor(lags/Fsd,(1:length(up_trans8{d}))/length(up_trans8{d}),lf3t_up_ctrig_mat(up_order,:));shading flat
%         hold on
%                 caxis([-1 3]);colorbar
%     plot(up_state_dur8{d}(up_order),(1:length(up_trans8{d}))/length(up_trans8{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     caxis([-1 2]);colorbar
%         tname = ['C:\WC_Germany\JMM_Analysis_Pyr\newest_lf8_up_trig_lf8_lf3_lf2\high_freq_lf3_corrected_' f_names{d}];
%         print('-dpng',tname);
%         close

%         
%         pcolor(dlags/Fsd,1:length(down_trans{d}),lf5_down_ctrig_mat(up_order,:));shading flat
%         tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_down_ctrig_all\lf5_' f_names{d}];
%         print('-dpng',tname);
%         close
% 
%                 pcolor(dlags/Fsd,1:length(down_trans{d}),lf7_down_ctrig_mat(up_order,:));shading flat
%         tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_down_ctrig_all\lf7_' f_names{d}];
%         print('-dpng',tname);
%         close

%     clear near_lfp_ind
%     clear prev_lfp_down

%     clear wcv_d lf2_d lf3_d lf5_d lf7_d lf8_d
    

end