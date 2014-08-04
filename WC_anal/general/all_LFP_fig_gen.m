load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


figure
subplot(2,1,1)
plot(corr_grid,stell_corr2_amp_hist,'linewidth',2)
hold on
plot(corr_grid,stell_corr3_amp_hist,'r','linewidth',2)
% plot(corr_grid,stell_corr4_amp_hist,'g','linewidth',2)
plot(corr_grid,stell_corr5_amp_hist,'c','linewidth',2)
plot(corr_grid,stell_corr7_amp_hist,'k','linewidth',2)
plot(corr_grid,stell_corr8_amp_hist,'y','linewidth',2)
title('Stellate Peak XCorr Distributions')
legend('2','3','5','7','8')
subplot(2,1,2)
plot(corr_grid,pyr_corr2_amp_hist,'linewidth',2)
hold on
plot(corr_grid,pyr_corr3_amp_hist,'r','linewidth',2)
% plot(corr_grid,pyr_corr4_amp_hist,'g','linewidth',2)
plot(corr_grid,pyr_corr5_amp_hist,'c','linewidth',2)
plot(corr_grid,pyr_corr7_amp_hist,'k','linewidth',2)
plot(corr_grid,pyr_corr8_amp_hist,'y','linewidth',2)
title('Pyramidal Peak XCorr Distributions')
legend('2','3','5','7','8')



figure
subplot(2,1,1)
plot(corr_grid,stell_corr2_coeff_hist,'linewidth',2)
hold on
plot(corr_grid,stell_corr3_coeff_hist,'r','linewidth',2)
% plot(corr_grid,stell_corr4_coeff_hist,'g','linewidth',2)
plot(corr_grid,stell_corr5_coeff_hist,'c','linewidth',2)
plot(corr_grid,stell_corr7_coeff_hist,'k','linewidth',2)
plot(corr_grid,stell_corr8_coeff_hist,'y','linewidth',2)
title('Stellate Corr Coeff Distributions')
legend('2','3','5','7','8')
subplot(2,1,2)
plot(corr_grid,pyr_corr2_coeff_hist,'linewidth',2)
hold on
plot(corr_grid,pyr_corr3_coeff_hist,'r','linewidth',2)
% plot(corr_grid,pyr_corr4_coeff_hist,'g','linewidth',2)
plot(corr_grid,pyr_corr5_coeff_hist,'c','linewidth',2)
plot(corr_grid,pyr_corr7_coeff_hist,'k','linewidth',2)
plot(corr_grid,pyr_corr8_coeff_hist,'y','linewidth',2)
title('Pyramidal Corr Coeff Distributions')
legend('2','3','5','7','8')



figure
subplot(2,1,1)
plot(corr_grid,stell_corr2_off_hist,'linewidth',2)
hold on
plot(corr_grid,stell_corr3_off_hist,'r','linewidth',2)
% plot(corr_grid,stell_corr4_off_hist,'g','linewidth',2)
plot(corr_grid,stell_corr5_off_hist,'c','linewidth',2)
plot(corr_grid,stell_corr7_off_hist,'k','linewidth',2)
plot(corr_grid,stell_corr8_off_hist,'y','linewidth',2)
title('Stellate Corr off Distributions')
legend('2','3','5','7','8')
subplot(2,1,2)
plot(corr_grid,pyr_corr2_off_hist,'linewidth',2)
hold on
plot(corr_grid,pyr_corr3_off_hist,'r','linewidth',2)
% plot(corr_grid,pyr_corr4_off_hist,'g','linewidth',2)
plot(corr_grid,pyr_corr5_off_hist,'c','linewidth',2)
plot(corr_grid,pyr_corr7_off_hist,'k','linewidth',2)
plot(corr_grid,pyr_corr8_off_hist,'y','linewidth',2)
title('Pyramidal Corr off Distributions')
legend('2','3','5','7','8')



figure
plot(dci_grid,total_lf2_hist,'linewidth',2)
hold on
plot(dci_grid,total_lf3_hist,'r','linewidth',2)
% plot(dci_grid,total_lf4_hist,'g','linewidth',2)
plot(dci_grid,total_lf5_hist,'c','linewidth',2)
plot(dci_grid,total_lf7_hist,'k','linewidth',2)
plot(dci_grid,total_lf8_hist,'y','linewidth',2)
plot(dci_grid,pyr_wcv_dci_hist,'linewidth',4)
plot(dci_grid,stell_wcv_dci_hist,'r','linewidth',4)
title('LFP DCI Distributions')
legend('2','3','5','7','8','Pyr','Stell')


figure
plot(period_grid,fgsmooth(pyr_wcv_period_hist,3),'linewidth',4)
hold on
plot(period_grid,fgsmooth(stell_wcv_period_hist,3),'r','linewidth',4)
plot(period_grid,fgsmooth(total_lf2_period_hist,3),'linewidth',2)
plot(period_grid,fgsmooth(total_lf3_period_hist,3),'r','linewidth',2)
% plot(period_grid,fgsmooth(total_lf4_period_hist,3),'g','linewidth',2)
plot(period_grid,fgsmooth(total_lf5_period_hist,3),'c','linewidth',2)
plot(period_grid,fgsmooth(total_lf7_period_hist,3),'k','linewidth',2)
plot(period_grid,fgsmooth(total_lf8_period_hist,3),'y','linewidth',2)
title('LFP Period Distributions')
legend('Pyr','Stell','2','3','5','7','8')

