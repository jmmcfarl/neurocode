clear all
% 
load C:\WC_Germany\JMM_Analysis_pyr\sim_record_UDS_dur_1
dsf = 2;
params.Fs = 2016/dsf;
params.err = [1 .01];
params.fpass = [0 400];
params.tapers = [1 1];
bad_f = [50 100 150 200 250 300 350 400];
% W = 0.04;
win = 100;
% niqf = 2016/2;
% hcf = 450/niqf;
% [b,a] = butter(2,hcf,'low');

movingwin = [0.5 0.5];
    
    load used_data lf15 lf8
      
%     T = floor(length(lf8)/2016);
%     TW = floor(T*W);
%     
%     params.tapers = [TW 100];
%     
%     
%     disp(sprintf('T = %0.5g K = %0.5g',T,(2*TW-1)))
  
    
    
    
%     lf15 = filtfilt(b,a,lf15);
%     lf8 = filtfilt(b,a,lf8);   

    down_15 = downsample(lf15,dsf);
    down_8 = downsample(lf8,dsf);
    down_15 = zscore(down_15);
    down_8 = zscore(down_8);
   
    up_trans8 = up_trans8*4;
    down_trans8 = down_trans8*4;
    up_trans = up_trans*4;
    down_trans = down_trans*4;
    
%     up_trans8(850:end) = [];
%     down_trans8(850:end) = [];
%     up_trans(820:end) = [];
%     down_trans(820:end) = [];
    
    
%     [C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencysegc(down_15,down_8,win,params);
%     [S8,f,varS,dummy,S8err]=mtspectrumsegc(down_8,win,params);
%     [S15,f,varS,dummy,S15err] = mtspectrumsegc(down_15,win,params);
up_markers8 = [up_trans8' down_trans8'];
up_markers15 = [up_trans' down_trans'];
down_markers8 = [down_trans8(1:end-1)' up_trans8(2:end)'];
down_markers15 = [down_trans(1:end-1)' down_trans(2:end)'];
%     [S8_up,f,Serr8_up]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, up_markers8);
%      [S8_down,f,Serr8_down]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, down_markers8);
%     [S15_up,f,Serr15_up]= mtspectrumc_unequal_length_trials(down_15, movingwin, params, up_markers15);
%      [S15_down,f,Serr15_down]= mtspectrumc_unequal_length_trials(down_15, movingwin, params, down_markers15);
     [S8_up,f]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, up_markers8);
     [S8_down,f]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, down_markers8);
    [S15_up,f]= mtspectrumc_unequal_length_trials(down_15, movingwin, params, up_markers15);
     [S15_down,f]= mtspectrumc_unequal_length_trials(down_15, movingwin, params, down_markers15);
 
     exclude_f = [];
     for i = 1:length(bad_f)
        exclude_f = [exclude_f find(abs(f-bad_f(i)) < 5)];
     end
     S8_up(exclude_f) = nan;
     S8_down(exclude_f) = nan;
     S15_up(exclude_f) = nan;
     S15_down(exclude_f) = nan;
     
        subplot(2,1,1)
    plot(f,10*log10(S8_up),'linewidth',2)
    hold on
    plot(f,10*log10(S8_down),'r','linewidth',2)
    legend('Up','Down')
%     plot(f,10*log10(Serr8_up(1,:)),'--')
%     plot(f,10*log10(Serr8_up(2,:)),'--')
%     plot(f,10*log10(Serr8_down(1,:)),'r--')
%     plot(f,10*log10(Serr8_down(2,:)),'r--')
    title('LFP')
    subplot(2,1,2)
    plot(f,10*log10(S15_up),'linewidth',2)
    hold on
    plot(f,10*log10(S15_down),'r','linewidth',2)
    legend('Up','Down')
%     plot(f,10*log10(Serr15_up(1,:)),'--')
%     plot(f,10*log10(Serr15_up(2,:)),'--')
%     plot(f,10*log10(Serr15_down(1,:)),'r--')
%     plot(f,10*log10(Serr15_down(2,:)),'r--')
    title('LF15')
    
    pow_diff_8 = 10*log10(S8_up)-10*log10(S8_down);
    pow_diff_15 = 10*log10(S15_up)-10*log10(S15_down);
    
    figure
    subplot(2,1,1)
    plot(f,10*log10(S8_up)-10*log10(S8_down),'linewidth',2)
    hold on
    plot(f,10*log10(S15_up)-10*log10(S15_down),'r','linewidth',2)
    legend('LFP','LF15')
    title('Up-Down power')
    subplot(2,1,2)
    plot(f,pow_diff_8-pow_diff_15,'k','linewidth',2)

%% COMPARE SIM DATA TO OVERALL MP
%     load C:\WC_Germany\JMM_Analysis_pyr\overall_spec
%     
% Swm = mean(10*log10(Sw));
% Swe = std(10*log10(Sw))/sqrt(17);
% S8m = mean(10*log10(S8));
% S8e = std(10*log10(S8))/sqrt(18);
% 
% Swl = Swm-2*Swe;
% Swu = Swm+2*Swe;
% S8l = S8m-2*S8e;
% S8u = S8m+2*S8e;
% 
% Cm = mean(C);
% Ce = std(C)/sqrt(17);
% Cu = Cm+2*Ce;
% Cl = Cm-2*Ce;
% 
% plot(f,S8m,'r','linewidth',2)
% hold on
% plot(f,Swm,'b','linewidth',2)
% load C:\WC_Germany\JMM_Analysis_pyr\sim_record_B_spec_data
% plot(f,10*log10(S8),'k','linewidth',2)
% plot(f,10*log10(S15),'c','linewidth',2)
% legend('Average LFP8','Avg MP','Session LFP8','Session LFP15')
% plot(f,10*log10(S8err(1,:)),'k--')
% plot(f,10*log10(S8err(2,:)),'k--')
% plot(f,10*log10(S15err(1,:)),'c--')
% plot(f,10*log10(S15err(2,:)),'c--')
% xlabel('Frequency (Hz)','FontSize',14)
% ylabel('Power (dB)','FontSize',14)
% xlim([0 1])
% % figure(1)
% % plot(f,Cm,'linewidth',2)
% % hold on
% % plot(f,Cu,'--')
% % plot(f,Cl,'--')
% % xlabel('Frequency (Hz)','FontSize',14)
% % ylabel('Coherency','FontSize',14)
% % xlim([0 1])
% % av_conf = mean(confC);
% % min_conf = min(confC);
% % max_conf = max(confC);
% % line([0 10],[av_conf av_conf],'Color','k')
% % line([0 10],[min_conf min_conf],'Color','k','linestyle','--')
% % line([0 10],[max_conf max_conf],'Color','k','linestyle','--')
% % 
% % %now load sim recording data
% % load C:\WC_Germany\JMM_Analysis_pyr\sim_record_A_spec_data
% % plot(f,C,'r','linewidth',2)
% % plot(f,Cerr(1,:),'r--')
% % plot(f,Cerr(2,:),'r--')
% % line([0 10],[confC confC],'Color','r')
% % load C:\WC_Germany\JMM_Analysis_pyr\sim_record_B_spec_data
% % plot(f,C,'g','linewidth',2)
% % plot(f,Cerr(1,:),'g--')
% % plot(f,Cerr(2,:),'g--')
% % line([0 10],[confC confC],'Color','c')
% % xlim([0 2])
% % ylim([0 1])
% 
% 
% 
% 
% 
% 
% %% Do the same plot as above, but smooth everything
% smooth_fac = 10;
% 
%     load C:\WC_Germany\JMM_Analysis_pyr\overall_spec
%     
% Swm = mean(10*log10(Sw));
% Swe = std(10*log10(Sw))/sqrt(17);
% S8m = mean(10*log10(S8));
% S8e = std(10*log10(S8))/sqrt(18);
% 
% Swl = Swm-2*Swe;
% Swu = Swm+2*Swe;
% S8l = S8m-2*S8e;
% S8u = S8m+2*S8e;
% 
% Cm = mean(C);
% Ce = std(C)/sqrt(17);
% Cu = Cm+2*Ce;
% Cl = Cm-2*Ce;
% 
% Cm = jmm_smooth_1d(Cm,smooth_fac);
% Cu = jmm_smooth_1d(Cu,smooth_fac);
% Cl = jmm_smooth_1d(Cl,smooth_fac);
% 
% % figure(1)
% % plot(f,Cm,'linewidth',2)
% % hold on
% % plot(f,Cu,'--')
% % plot(f,Cl,'--')
% % xlabel('Frequency (Hz)','FontSize',14)
% % ylabel('Coherency','FontSize',14)
% % xlim([0 1])
% % av_conf = mean(confC);
% % min_conf = min(confC);
% % max_conf = max(confC);
% % line([0 10],[av_conf av_conf],'Color','k')
% % line([0 10],[min_conf min_conf],'Color','k','linestyle','--')
% % line([0 10],[max_conf max_conf],'Color','k','linestyle','--')
% 
% %now load sim recording data
% load C:\WC_Germany\JMM_Analysis_pyr\sim_record_A_spec_data
% 
% C = jmm_smooth_1d(C,smooth_fac);
% Cerr(1,:) = jmm_Smooth_1d(Cerr(1,:),smooth_fac);
% Cerr(2,:) = jmm_smooth_1d(Cerr(2,:),smooth_fac);
% 
% plot(f,C,'r','linewidth',2)
% hold on
% plot(f,Cerr(1,:),'r--')
% plot(f,Cerr(2,:),'r--')
% line([0 30],[confC confC],'Color','r')
% load C:\WC_Germany\JMM_Analysis_pyr\sim_record_B_spec_data
% 
% C = jmm_smooth_1d(C,smooth_fac);
% Cerr(1,:) = jmm_Smooth_1d(Cerr(1,:),smooth_fac);
% Cerr(2,:) = jmm_smooth_1d(Cerr(2,:),smooth_fac);
% 
% plot(f,C,'g','linewidth',2)
% plot(f,Cerr(1,:),'g--')
% plot(f,Cerr(2,:),'g--')
% line([0 30],[confC confC],'Color','k')
% xlim([0 30])
% ylim([0 1])
% xlabel('Frequency (Hz)','FontSize',14)
% ylabel('Coherency','FontSize',14)
