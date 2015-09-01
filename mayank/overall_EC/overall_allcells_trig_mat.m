clear all
close all

load G:\WC_Germany\overall_EC\overall_EC_dir
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

drive_letter = 'G';

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 10/niqf;
[b,a] = butter(2,[lcf hcf]);
lcf2 = 1/niqf;
hcf2 = 10/niqf;
[b2,a2] = butter(2,[lcf2 hcf2]);
lcf3 = 15/niqf;
hcf3 = 80/niqf;
[b3,a3] = butter(2,[lcf3 hcf3]);
dsf = 8;
Fsd = Fs/dsf;
pow_smooth = round(Fsd*0.05);

backlag = 4*Fsd;
forwardlag = 10*Fsd;
lags = (-backlag:forwardlag)/Fsd;


for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name,'_atr');
    
    load used_data lf8 wcv_minus_spike lf3 
        
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_ft = filtfilt(b2,a2,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    lf3_f = filtfilt(b,a,lf3);
%     lf2 = lf2/sess_data(d).gains(2);
%     lf5 = lf5/sess_data(d).gains(5);
%     lf2_r = lf2-lf5;
    wcv_f = downsample(wcv_f,dsf)/sess_data(d).gains(1);
    wcv_ft = downsample(wcv_ft,dsf)/sess_data(d).gains(1);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);
    lf3_f = downsample(lf3_f,dsf)/sess_data(d).gains(3);
%     lf2_r = filtfilt(b,a,lf2_r);
%     lf2_r = downsample(lf2_r,dsf);

%     lf2_hf = downsample(filtfilt(b3,a3,lf2),dsf);
%     lf2_hf = zscore(sqrt(jmm_smooth_1d_cor(lf2_hf.^2,pow_smooth)));
    
%     bad_lf2r(d) = max(isnan(lf2_r));
%     if bad_lf2r(d) == 0 && exist('ec_hmm_state_seq2r.mat')
    
    wcv_f = zscore(wcv_f);
    wcv_ft = zscore(wcv_ft);
%     lf2_f = zscore(lf2_f);
%     lf2_r = zscore(lf2_r);
    lf3_f = zscore(lf3_f);
    lf8_f = zscore(lf8_f);
    
    t_axis = (1:length(lf8_f))/Fsd;
    
    %% extract up and down transition times for MP and LF8
    load ec_hmm_state_seq
    load ec_hmm_state_seq8
%     load ec_hmm_state_seq2r
    mp_state_seq_c =  hmm_bbstate_seq;
    lf8_state_seq_c = hmm_bbstate_seq8;
%     lf2r_state_seq_c = hmm_bbstate_seq2r;
    [new_mp_seg_inds] = round(resample_uds_seg_inds(hmm.UDS_segs,50.4,Fsd,length(t_axis)));
    [new_lf8_seg_inds] = round(resample_uds_seg_inds(hmm8.UDS_segs,50.4,Fsd,length(t_axis)));
%     [new_lf2r_seg_inds] = round(resample_uds_seg_inds(hmm2r.UDS_segs,50.4,Fsd,length(t_axis)));
    
    mp_state_seq = nan(size(t_axis));
    lf8_state_seq = nan(size(t_axis));
%     lf2r_state_seq = nan(size(t_axis));
    
    mp_utrans = [];
    mp_dtrans = [];
    for n = 1:hmm.Nsegs
        mp_state_seq(new_mp_seg_inds(n,1):new_mp_seg_inds(n,2)) = mp_state_seq_c{n};
        cur_mp_utrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 1 & mp_state_seq_c{n}(2:end) == 2);
        cur_mp_dtrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 2 & mp_state_seq_c{n}(2:end) == 1);
        cur_mp_dtrans(cur_mp_dtrans < cur_mp_utrans(1)) = [];
        cur_mp_utrans(cur_mp_utrans > cur_mp_dtrans(end)) = [];
        mp_utrans = [mp_utrans; cur_mp_utrans];
        mp_dtrans = [mp_dtrans; cur_mp_dtrans];
    end
    
    lf8_utrans = [];
    lf8_dtrans = [];
    for n = 1:hmm8.Nsegs
        lf8_state_seq(new_lf8_seg_inds(n,1):new_lf8_seg_inds(n,2)) = lf8_state_seq_c{n};
        cur_lf8_utrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 1 & lf8_state_seq_c{n}(2:end) == 2);
        cur_lf8_dtrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 2 & lf8_state_seq_c{n}(2:end) == 1);
        cur_lf8_dtrans(cur_lf8_dtrans < cur_lf8_utrans(1)) = [];
        cur_lf8_utrans(cur_lf8_utrans > cur_lf8_dtrans(end)) = [];
        lf8_utrans = [lf8_utrans; cur_lf8_utrans];
        lf8_dtrans = [lf8_dtrans; cur_lf8_dtrans];
    end
% 
%     lf2r_utrans = [];
%     lf2r_dtrans = [];
%     for n = 1:hmm2r.Nsegs
%         lf2r_state_seq(new_lf2r_seg_inds(n,1):new_lf2r_seg_inds(n,2)) = lf2r_state_seq_c{n};
%         cur_lf2r_utrans = new_lf2r_seg_inds(n,1) + find(lf2r_state_seq_c{n}(1:end-1) == 1 & lf2r_state_seq_c{n}(2:end) == 2);
%         cur_lf2r_dtrans = new_lf2r_seg_inds(n,1) + find(lf2r_state_seq_c{n}(1:end-1) == 2 & lf2r_state_seq_c{n}(2:end) == 1);
%         cur_lf2r_dtrans(cur_lf2r_dtrans < cur_lf2r_utrans(1)) = [];
%         cur_lf2r_utrans(cur_lf2r_utrans > cur_lf2r_dtrans(end)) = [];
%         lf2r_utrans = [lf2r_utrans; cur_lf2r_utrans];
%         lf2r_dtrans = [lf2r_dtrans; cur_lf2r_dtrans];
%     end
% 
    n_mp_ups = length(mp_utrans);
%     n_lf2r_ups = length(lf2r_utrans);
    n_lf8_ups = length(lf8_utrans);
    
    %% initialize
    mp_utrig_mp_mat = nan(n_mp_ups,length(lags));
    mp_utrig_mpt_mat = nan(n_mp_ups,length(lags));
    mp_utrig_lf8_mat = nan(n_mp_ups,length(lags));
    mp_utrig_lf3_mat = nan(n_mp_ups,length(lags));
%     mp_utrig_lf2hf_mat = nan(n_mp_ups,length(lags));
%     mp_utrig_lf2r_mat = nan(n_mp_ups,length(lags));
    mp_dtrig_mp_mat = nan(n_mp_ups,length(lags));
    mp_dtrig_lf8_mat = nan(n_mp_ups,length(lags));
    mp_dtrig_lf3_mat = nan(n_mp_ups,length(lags));
%     mp_dtrig_lf2_mat = nan(n_mp_ups,length(lags));
%     mp_dtrig_lf2r_mat = nan(n_mp_ups,length(lags));
    
%     calculate mp utrigs
    for i = 1:n_mp_ups       
        if mp_utrans(i) > backlag && length(wcv_f) - mp_utrans(i) > forwardlag            
            mp_utrig_mp_mat(i,:) = wcv_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
            mp_utrig_mpt_mat(i,:) = wcv_ft(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
            mp_utrig_lf8_mat(i,:) = lf8_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
            mp_utrig_lf3_mat(i,:) = lf3_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);            
%             mp_utrig_lf2hf_mat(i,:) = lf2_hf(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);            
%             mp_utrig_lf2r_mat(i,:) = lf2_r(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);            
        end       
    end
    
%     calculate mp dtrigs
    for i = 1:n_mp_ups     
        if mp_dtrans(i) > backlag && length(wcv_f) - mp_dtrans(i) > forwardlag          
            mp_dtrig_mp_mat(i,:) = wcv_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
            mp_dtrig_lf8_mat(i,:) = lf8_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
            mp_dtrig_lf3_mat(i,:) = lf3_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);            
%             mp_dtrig_lf2_mat(i,:) = lf2_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);            
%             mp_dtrig_lf2r_mat(i,:) = lf2_r(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);            
        end        
    end
%      
    mp_utrig_mp(d,:) = nanmean(mp_utrig_mp_mat);
    mp_utrig_mpt(d,:) = nanmean(mp_utrig_mpt_mat);
    mp_utrig_lf8(d,:) = nanmean(mp_utrig_lf8_mat);
    mp_utrig_lf3(d,:) = nanmean(mp_utrig_lf3_mat);
%     mp_utrig_lf2hf(d,:) = nanmean(mp_utrig_lf2hf_mat);
%     mp_utrig_lf2r(d,:) = nanmean(mp_utrig_lf2r_mat);
    mp_dtrig_mp(d,:) = nanmean(mp_dtrig_mp_mat);
    mp_dtrig_lf8(d,:) = nanmean(mp_dtrig_lf8_mat);
    mp_dtrig_lf3(d,:) = nanmean(mp_dtrig_lf3_mat);
%     mp_dtrig_lf2(d,:) = nanmean(mp_dtrig_lf2_mat);
%     mp_dtrig_lf2r(d,:) = nanmean(mp_dtrig_lf2r_mat);
    
    mp_updur = (mp_dtrans-mp_utrans)/Fsd;
    mp_downdur = (mp_utrans(2:end)-mp_dtrans(1:end-1))/Fsd;
    
    [dummy,up_order] = sort(mp_downdur);    
%     
% %      plot mp up trig matrices
%     Fig= figure('visible','off');
%     set(Fig,'PaperUnits','centimeters');
%     set(Fig, 'PaperSize', [15 20]);
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
%     subplot(3,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-2.3 2.]);colorbar
%     xlim([-2 6])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('MP')
%     subplot(3,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mpt_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 6])
%     line([0 0],[0 1],'Color','k')    
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('Cortical LFP')
%     subplot(3,1,3)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf8_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.7 1.7]);colorbar
%     xlim([-2 6])
%     line([0 0],[0 1],'Color','k')    
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('Hippocampal LFP')
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\mp_utrig_theta_' s_name];
%     print('-dpng',t_names);
%     close

     
%          Fig= figure('visible','off');
%     set(Fig,'PaperUnits','centimeters');
%     set(Fig, 'PaperSize', [15 20]);
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
%     subplot(3,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-2.3 2.]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('MP')
%     subplot(3,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf8_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 2.5]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')    
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('Cortical LFP')
%     subplot(3,1,3)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf2hf_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.7 1.7]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')    
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('Hippocampal LFP')
%     t_names = ['G:\WC_Germany\overall_EC\trig_mats\mp_utrig_hf_' s_name];
%     print('-dpng',t_names);
%     close
% 
%     figure('visible','off')
%     subplot(2,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-3 3]);colorbar
%     xlim([-2 10])
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf8_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-3 3]);colorbar
%     xlim([-2 10])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\mp_utrig_lf8_' s_name];
%     print('-dpng',t_names);
%     close
% 
%     figure('visible','off')
%     subplot(2,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 10])
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf3_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 10])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\mp_utrig_lf3_' s_name];
%     print('-dpng',t_names);
%     close
% 
%     figure('visible','off')
%         subplot(2,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 10])
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf2_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 10])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\mp_utrig_lf2_' s_name];
%     print('-dpng',t_names);
%     close
% 
%      figure('visible','off')
%        subplot(2,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 10])
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf2r_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 10])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\mp_utrig_lf2r_' s_name];
%     print('-dpng',t_names);
%     close

%% plot mp down trig matrices

n_mp_ups = n_mp_ups - 1;
    figure('visible','off')
    subplot(2,1,1)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_mp_mat(up_order,:));
    shading flat; colorbar; hold on
    caxis([-3 3]);colorbar
    xlim([-4 4])
    line([0 0],[0 1],'Color','k')
    subplot(2,1,2)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_lf8_mat(up_order,:));
    shading flat; colorbar; hold on
    caxis([-3 3]);colorbar
    xlim([-4 4])
    line([0 0],[0 1],'Color','k')    
    t_names = ['G:\WC_Germany\overall_EC\trig_mats\mp_dtrig_lf8_' s_name];
    print('-dpng',t_names);
    close
% 
    figure('visible','off')
    subplot(2,1,1)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_mp_mat(up_order,:));
    shading flat; colorbar; hold on
    caxis([-1.5 1.5]);colorbar
    xlim([-4 4])
    line([0 0],[0 1],'Color','k')
    subplot(2,1,2)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_lf3_mat(up_order,:));
    shading flat; colorbar; hold on
    caxis([-1.5 1.5]);colorbar
    xlim([-4 4])
    line([0 0],[0 1],'Color','k')    
    t_names = ['G:\WC_Germany\overall_EC\trig_mats\mp_dtrig_lf3_' s_name];
    print('-dpng',t_names);
    close
% 
%     figure('visible','off')
%         subplot(2,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     caxis([-1.5 1.5]);colorbar
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_lf2_mat(up_order,:));
%     shading flat; colorbar; hold on
%     caxis([-1.5 1.5]);colorbar
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\mp_dtrig_lf2_' s_name];
%     print('-dpng',t_names);
%     close
% 
%     figure('visible','off')
%         subplot(2,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     caxis([-1.5 1.5]);colorbar
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_lf2r_mat(up_order,:));
%     shading flat; colorbar; hold on
%     caxis([-1.5 1.5]);colorbar
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\mp_dtrig_lf2r_' s_name];
%     print('-dpng',t_names);
%     close
% 
    
        %% For LF2r triggered analysis initialize
%     lf2r_utrig_mp_mat = nan(n_lf2r_ups,length(lags));
%     lf2r_utrig_lf8_mat = nan(n_lf2r_ups,length(lags));
% %     lf2r_utrig_lf3_mat = nan(n_lf2r_ups,length(lags));
% %     lf2r_utrig_lf2_mat = nan(n_lf2r_ups,length(lags));
%     lf2r_utrig_lf2r_mat = nan(n_lf2r_ups,length(lags));
% %     lf2r_dtrig_mp_mat = nan(n_lf2r_ups,length(lags));
% %     lf2r_dtrig_lf8_mat = nan(n_lf2r_ups,length(lags));
% %     lf2r_dtrig_lf3_mat = nan(n_lf2r_ups,length(lags));
% %     lf2r_dtrig_lf2_mat = nan(n_lf2r_ups,length(lags));
% %     lf2r_dtrig_lf2r_mat = nan(n_lf2r_ups,length(lags));
%     
% %     calculate mp utrigs
%     for i = 1:n_lf2r_ups       
%         if lf2r_utrans(i) > backlag && length(wcv_f) - lf2r_utrans(i) > forwardlag            
%             lf2r_utrig_mp_mat(i,:) = wcv_f(lf2r_utrans(i)-backlag:lf2r_utrans(i)+forwardlag);
%             lf2r_utrig_lf8_mat(i,:) = lf8_f(lf2r_utrans(i)-backlag:lf2r_utrans(i)+forwardlag);
% %             lf2r_utrig_lf3_mat(i,:) = lf3_f(lf2r_utrans(i)-backlag:lf2r_utrans(i)+forwardlag);            
% %             lf2r_utrig_lf2_mat(i,:) = lf2_f(lf2r_utrans(i)-backlag:lf2r_utrans(i)+forwardlag);            
%             lf2r_utrig_lf2r_mat(i,:) = lf2_r(lf2r_utrans(i)-backlag:lf2r_utrans(i)+forwardlag);            
%         end       
%     end
%     
% %     %calculate mp dtrigs
% %     for i = 1:n_lf2r_ups     
% %         if lf2r_dtrans(i) > backlag && length(wcv_f) - lf2r_dtrans(i) > forwardlag          
% %             lf2r_dtrig_mp_mat(i,:) = wcv_f(lf2r_dtrans(i)-backlag:lf2r_dtrans(i)+forwardlag);
% %             lf2r_dtrig_lf8_mat(i,:) = lf8_f(lf2r_dtrans(i)-backlag:lf2r_dtrans(i)+forwardlag);
% %             lf2r_dtrig_lf3_mat(i,:) = lf3_f(lf2r_dtrans(i)-backlag:lf2r_dtrans(i)+forwardlag);            
% %             lf2r_dtrig_lf2_mat(i,:) = lf2_f(lf2r_dtrans(i)-backlag:lf2r_dtrans(i)+forwardlag);            
% %             lf2r_dtrig_lf2r_mat(i,:) = lf2_r(lf2r_dtrans(i)-backlag:lf2r_dtrans(i)+forwardlag);            
% %         end        
% %     end
% %      
%     lf2r_utrig_mp(d,:) = nanmean(lf2r_utrig_mp_mat);
%     lf2r_utrig_lf8(d,:) = nanmean(lf2r_utrig_lf8_mat);
% %     lf2r_utrig_lf3(d,:) = nanmean(lf2r_utrig_lf3_mat);
% %     lf2r_utrig_lf2(d,:) = nanmean(lf2r_utrig_lf2_mat);
%     lf2r_utrig_lf2r(d,:) = nanmean(lf2r_utrig_lf2r_mat);
% %     lf2r_dtrig_mp(d,:) = nanmean(lf2r_dtrig_mp_mat);
% %     lf2r_dtrig_lf8(d,:) = nanmean(lf2r_dtrig_lf8_mat);
% %     lf2r_dtrig_lf3(d,:) = nanmean(lf2r_dtrig_lf3_mat);
% %     lf2r_dtrig_lf2(d,:) = nanmean(lf2r_dtrig_lf2_mat);
% %     lf2r_dtrig_lf2r(d,:) = nanmean(lf2r_dtrig_lf2r_mat);
%     
%     lf2r_updur = (lf2r_dtrans-lf2r_utrans)/Fsd;
%     lf2r_downdur = (lf2r_utrans(2:end)-lf2r_dtrans(1:end-1))/Fsd;
%     
%     [dummy,up_order] = sort(lf2r_updur);    
%     
%     %% plot lf2r up trig matrices
%     Fig= figure('visible','off');
%     set(Fig,'PaperUnits','centimeters');
%     set(Fig, 'PaperSize', [15 20]);
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
%     subplot(3,1,1)
%     pcolor(lags,(1:n_lf2r_ups)/n_lf2r_ups,lf2r_utrig_lf2r_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(lf2r_updur(up_order),(1:n_lf2r_ups)/n_lf2r_ups,'w','linewidth',2)
%     caxis([-2 2.5]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since hipp Up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('Hippocampal LFP','fontsize',14)
%     subplot(3,1,2)
%     pcolor(lags,(1:n_lf2r_ups)/n_lf2r_ups,lf2r_utrig_lf8_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(lf2r_updur(up_order),(1:n_lf2r_ups)/n_lf2r_ups,'w','linewidth',2)
%     caxis([-1.7 2.5]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')    
%     xlabel('Time since hipp Up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('Cortical LFP','fontsize',14)
%     subplot(3,1,3)
%     pcolor(lags,(1:n_lf2r_ups)/n_lf2r_ups,lf2r_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(lf2r_updur(up_order),(1:n_lf2r_ups)/n_lf2r_ups,'w','linewidth',2)
%     caxis([-1.5 1.7]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')    
%     xlabel('Time since hipp Up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('MP','fontsize',14)
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\lf2r_utrig_all_' s_name];
%     print('-dpng',t_names);
%     close
% % % 
% %     figure('visible','off')
% %     subplot(2,1,1)
% %     pcolor(lags,(1:n_lf2r_ups)/n_lf2r_ups,lf2r_utrig_lf2r_mat(up_order,:));
% %     shading flat; colorbar; hold on
% %     plot(lf2r_updur(up_order),(1:n_lf2r_ups)/n_lf2r_ups,'w','linewidth',2)
% %     caxis([-2 2.5]);colorbar
% %     xlim([-2 10])
% %     line([0 0],[0 1],'Color','k')
% %     xlabel('Time since hipp Up-transition (s)','fontsize',14)
% %     ylabel('Up state number','fontsize',14)
% %     title('Hippocampal LFP','fontsize',14)
% %     subplot(2,1,2)
% %     pcolor(lags,(1:n_lf2r_ups)/n_lf2r_ups,lf2r_utrig_mp_mat(up_order,:));
% %     shading flat; colorbar; hold on
% %     plot(lf2r_updur(up_order),(1:n_lf2r_ups)/n_lf2r_ups,'w','linewidth',2)
% %     caxis([-2 2.5]);colorbar
% %     xlim([-2 10])
% %     line([0 0],[0 1],'Color','k')    
% %     xlabel('Time since hipp Up-transition (s)','fontsize',14)
% %     ylabel('Up state number','fontsize',14)
% %     title('MP','fontsize',14)
% %     t_names = ['F:\WC_Germany\overall_EC\trig_mats\lf2r_utrig_mp_' s_name];
% %     print('-dpng',t_names);
% %     close

%% plot lf2r down trig matrices
%     figure('visible','off')
%     subplot(2,1,1)
%     pcolor(lags,(1:n_lf2r_ups)/n_lf2r_ups,lf2r_dtrig_lf2r_mat(up_order,:));
%     shading flat; colorbar; hold on
%     caxis([-3 3]);colorbar
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_lf2r_ups)/n_lf2r_ups,lf2r_dtrig_lf8_mat(up_order,:));
%     shading flat; colorbar; hold on
%     caxis([-3 3]);colorbar
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\lf2r_dtrig_lf8_' s_name];
%     print('-dpng',t_names);
%     close
% 
%     figure('visible','off')
%     subplot(2,1,1)
%     pcolor(lags,(1:n_lf2r_ups)/n_lf2r_ups,lf2r_dtrig_lf2r_mat(up_order,:));
%     shading flat; colorbar; hold on
%     caxis([-1.5 1.5]);colorbar
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_lf2r_ups)/n_lf2r_ups,lf2r_dtrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     caxis([-1.5 1.5]);colorbar
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\lf2r_dtrig_mp_' s_name];
%     print('-dpng',t_names);
%     close

    
    
%             %% For LF8 triggered analysis initialize
%     lf8_utrig_mp_mat = nan(n_lf8_ups,length(lags));
%     lf8_utrig_lf8_mat = nan(n_lf8_ups,length(lags));
% %     lf8_utrig_lf3_mat = nan(n_lf8_ups,length(lags));
%     lf8_utrig_lf2hf_mat = nan(n_lf8_ups,length(lags));
% %     lf8_utrig_lf2r_mat = nan(n_lf8_ups,length(lags));
% %     lf8_dtrig_mp_mat = nan(n_lf8_ups,length(lags));
% %     lf8_dtrig_lf8_mat = nan(n_lf8_ups,length(lags));
% %     lf8_dtrig_lf3_mat = nan(n_lf8_ups,length(lags));
% %     lf8_dtrig_lf2_mat = nan(n_lf8_ups,length(lags));
% %     lf8_dtrig_lf2r_mat = nan(n_lf8_ups,length(lags));
%     
%     %calculate mp utrigs
%     for i = 1:n_lf8_ups       
%         if lf8_utrans(i) > backlag && length(wcv_f) - lf8_utrans(i) > forwardlag            
%             lf8_utrig_mp_mat(i,:) = wcv_f(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);
%             lf8_utrig_lf8_mat(i,:) = lf8_f(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);
% %             lf8_utrig_lf3_mat(i,:) = lf3_f(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);            
%             lf8_utrig_lf2hf_mat(i,:) = lf2_hf(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);            
% %             lf8_utrig_lf2r_mat(i,:) = lf2_r(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);            
%         end       
%     end
%     
% %     %calculate mp dtrigs
% %     for i = 1:n_lf8_ups     
% %         if lf8_dtrans(i) > backlag && length(wcv_f) - lf8_dtrans(i) > forwardlag          
% %             lf8_dtrig_mp_mat(i,:) = wcv_f(lf8_dtrans(i)-backlag:lf8_dtrans(i)+forwardlag);
% %             lf8_dtrig_lf8_mat(i,:) = lf8_f(lf8_dtrans(i)-backlag:lf8_dtrans(i)+forwardlag);
% %             lf8_dtrig_lf3_mat(i,:) = lf3_f(lf8_dtrans(i)-backlag:lf8_dtrans(i)+forwardlag);            
% %             lf8_dtrig_lf2_mat(i,:) = lf2_f(lf8_dtrans(i)-backlag:lf8_dtrans(i)+forwardlag);            
% %             lf8_dtrig_lf2r_mat(i,:) = lf2_r(lf8_dtrans(i)-backlag:lf8_dtrans(i)+forwardlag);            
% %         end        
% %     end
%      
%     lf8_utrig_mp(d,:) = nanmean(lf8_utrig_mp_mat);
%     lf8_utrig_lf8(d,:) = nanmean(lf8_utrig_lf8_mat);
% %     lf8_utrig_lf3(d,:) = nanmean(lf8_utrig_lf3_mat);
%     lf8_utrig_lf2hf(d,:) = nanmean(lf8_utrig_lf2hf_mat);
% %     lf8_utrig_lf2r(d,:) = nanmean(lf8_utrig_lf2r_mat);
% %     lf8_dtrig_mp(d,:) = nanmean(lf8_dtrig_mp_mat);
% %     lf8_dtrig_lf8(d,:) = nanmean(lf8_dtrig_lf8_mat);
% %     lf8_dtrig_lf3(d,:) = nanmean(lf8_dtrig_lf3_mat);
% %     lf8_dtrig_lf2(d,:) = nanmean(lf8_dtrig_lf2_mat);
% %     lf8_dtrig_lf2r(d,:) = nanmean(lf8_dtrig_lf2r_mat);
%     
%     lf8_updur = (lf8_dtrans-lf8_utrans)/Fsd;
%     lf8_downdur = (lf8_utrans(2:end)-lf8_dtrans(1:end-1))/Fsd;
%     
%     [dummy,up_order] = sort(lf8_updur);    
    
%     %% plot lf8 up trig matrices
%     figure('visible','off')
%     subplot(3,1,1)
%     pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_lf8_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
%     caxis([-3 3]);colorbar
%     xlim([-2 3])
%     line([0 0],[0 1],'Color','k')
%     subplot(3,1,2)
%     pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
%     caxis([-3 3]);colorbar
%     xlim([-2 3])
%     line([0 0],[0 1],'Color','k')    
%     subplot(3,1,2)
%     pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_lf2hf_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
%     caxis([-3 3]);colorbar
%     xlim([-2 3])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['G:\WC_Germany\overall_EC\trig_mats\lf8_utrig_hf_' s_name];
%     print('-dpng',t_names);
%     close
% 
%     figure('visible','off')
%     subplot(2,1,1)
%     pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_lf8_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 3])
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_lf2r_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 3])
%     line([0 0],[0 1],'Color','k')    
%     t_names = ['F:\WC_Germany\overall_EC\trig_mats\lf8_utrig_lf2r_' s_name];
%     print('-dpng',t_names);
%     close


    end
% end

clear *mat
cd G:\WC_Germany\persistent_2010\
