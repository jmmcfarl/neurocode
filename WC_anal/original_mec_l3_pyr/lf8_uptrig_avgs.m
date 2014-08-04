clear all

load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28
% load C:\WC_Germany\JMM_analysis_pyr\WCV_LFP_up_trans_sig_fit_data
load C:\WC_Germany\JMM_analysis_pyr\lf8_period_f_data

niqf = 2016/2;
lcf = 0.1/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);
% lcf2 = 4/niqf;
% hcf2 = 12/niqf;
dsf = 8;
Fsd = 2016/dsf;
backlag = 5*Fsd;
forwardlag = 5*Fsd
lags = -backlag:forwardlag;
dlags = -backlag:backlag;

for d = 1:length(dir_array)
    cd(dir_array{d})
    pwd

    load used_data lf8 lf2 lf3 lf5 lf7 wcv_minus_spike

    lf8_f = filtfilt(b,a,lf8);
    lf2_f = filtfilt(b,a,lf2);
    lf3_f = filtfilt(b,a,lf3);
    lf5_f = filtfilt(b,a,lf5);
    lf7_f = filtfilt(b,a,lf7);
    wcv_f = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8_f,dsf);
    lf2_d = downsample(lf2_f,dsf);
    lf3_d = downsample(lf3_f,dsf);
    lf5_d = downsample(lf5_f,dsf);
    lf7_d = downsample(lf7_f,dsf);
    wcv_d = downsample(wcv_f,dsf);

%     lf8_d = jmm_smooth_1d(lf8_d.^2,10);
%     lf2_d = jmm_smooth_1d(lf2_d.^2,10);
%      lf3_d = jmm_smooth_1d(lf3_d.^2,10);
%        lf5_d = jmm_smooth_1d(lf5_d.^2,10);
%     lf7_d = jmm_smooth_1d(lf7_d.^2,10);

    lf8_d = zscore(lf8_d);
    lf2_d = zscore(lf2_d);
    lf3_d = zscore(lf3_d);
    lf5_d = zscore(lf5_d);
    lf7_d = zscore(lf7_d);
    wcv_d = zscore(wcv_d);
    
    clear *_f clear lf8 lf2 lf3 lf5 lf7

    wcv_up_ctrig_mat = zeros(length(up_trans8{d}),length(lags));
    lf8_up_ctrig_mat = wcv_up_ctrig_mat;
    lf2_up_ctrig_mat = wcv_up_ctrig_mat;
    lf3_up_ctrig_mat = wcv_up_ctrig_mat;
    lf5_up_ctrig_mat = wcv_up_ctrig_mat;
    lf7_up_ctrig_mat = wcv_up_ctrig_mat;
    

    for i = 1:length(up_trans8{d})

        
        if up_trans8{d}(i) > backlag & up_trans8{d}(i) < length(wcv_d)-forwardlag
            wcv_up_ctrig_mat(i,:) = wcv_d(up_trans8{d}(i)-backlag:up_trans8{d}(i)+forwardlag);
            lf8_up_ctrig_mat(i,:) = lf8_d(up_trans8{d}(i)-backlag:up_trans8{d}(i)+forwardlag);
            lf2_up_ctrig_mat(i,:) = lf2_d(up_trans8{d}(i)-backlag:up_trans8{d}(i)+forwardlag);
            lf3_up_ctrig_mat(i,:) = lf3_d(up_trans8{d}(i)-backlag:up_trans8{d}(i)+forwardlag);
            lf5_up_ctrig_mat(i,:) = lf5_d(up_trans8{d}(i)-backlag:up_trans8{d}(i)+forwardlag);
            lf7_up_ctrig_mat(i,:) = lf7_d(up_trans8{d}(i)-backlag:up_trans8{d}(i)+forwardlag);


        else
            wcv_up_ctrig_mat(i,:) = nan;
            lf8_up_ctrig_mat(i,:) = nan;
            lf2_up_ctrig_mat(i,:) = nan;
            lf3_up_ctrig_mat(i,:) = nan;
            lf5_up_ctrig_mat(i,:) = nan;
            lf7_up_ctrig_mat(i,:) = nan;

        end

    end


    
    
wcv_up_ctrig_avg(d,:) = nanmean(wcv_up_ctrig_mat);
lf2_up_ctrig_avg(d,:) = nanmean(lf2_up_ctrig_mat);
lf3_up_ctrig_avg(d,:) = nanmean(lf3_up_ctrig_mat);
lf5_up_ctrig_avg(d,:) = nanmean(lf5_up_ctrig_mat);
lf7_up_ctrig_avg(d,:) = nanmean(lf7_up_ctrig_mat);
lf8_up_ctrig_avg(d,:) = nanmean(lf8_up_ctrig_mat);
  
% [dummy,up_order] = sort(up_state_dur8{d});

%     pcolor(lags/Fsd,(1:length(up_trans8{d}))/length(up_trans8{d}) ...
%         ,wcv_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur8{d}(up_order),(1:length(up_trans8{d}))/length(up_trans8{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_uptrig_avgs\wcv_' f_names{d}];
%     print('-dpng',tname);
% close
% 
%     pcolor(lags/Fsd,(1:length(up_trans8{d}))/length(up_trans8{d}) ...
%         ,lf8_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur8{d}(up_order),(1:length(up_trans8{d}))/length(up_trans8{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_uptrig_avgs\lf8_' f_names{d}];
%     print('-dpng',tname);
% close
% 
%     pcolor(lags/Fsd,(1:length(up_trans8{d}))/length(up_trans8{d}) ...
%         ,lf2_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur8{d}(up_order),(1:length(up_trans8{d}))/length(up_trans8{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_uptrig_avgs\lf2_' f_names{d}];
%     print('-dpng',tname);
% close 
% 
%         pcolor(lags/Fsd,(1:length(up_trans8{d}))/length(up_trans8{d}) ...
%         ,lf3_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur8{d}(up_order),(1:length(up_trans8{d}))/length(up_trans8{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_uptrig_avgs\lf3_' f_names{d}];
%     print('-dpng',tname);
%     close
%     
%     pcolor(lags/Fsd,(1:length(up_trans8{d}))/length(up_trans8{d}) ...
%         ,lf7_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur8{d}(up_order),(1:length(up_trans8{d}))/length(up_trans8{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_uptrig_avgs\lf7_' f_names{d}];
%     print('-dpng',tname);
% close

    d
    %
    clear near_lfp_ind
    clear prev_lfp_down

    clear wcv_d lf2_d lf3_d lf5_d lf7_d lf8_d
    
end

cd C:\WC_Germany\JMM_analysis_pyr
save lf8_up_trig_avgs *_avg