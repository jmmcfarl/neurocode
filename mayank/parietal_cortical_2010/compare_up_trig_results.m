clear all
close all
cd G:\WC_Germany\parietal_cortical_2010
load up_trig_data

figure(1)
nfreqs = linspace(wfreqs(end),wfreqs(1),1000);
avg_trig_tfd_i = shiftdim(interp1(wfreqs,shiftdim(avg_trig8_tfd,1),nfreqs),2);
imagesc(lags(used_lags)/Fsd,nfreqs,squeeze(mean(avg_trig_tfd_i)));shading flat
line([0 0],[5 80],'Color','w')
colorbar

figure(2)
errorbar(lags/Fsd,mean(trig8_lf8),std(trig8_lf8)/sqrt(41)), hold on
hold on
xlim([-0.5 0.5])

figure(3)
errorbar(wfreqs,nanmean(freq_up_lag),nanstd(freq_up_lag)/sqrt(41)), hold on
xlim([5 80])

load up8_trig_data
figure(2)
errorbar(lags/Fsd,mean(trig8_lf8),std(trig8_lf8)/sqrt(41),'r')

figure(3)
errorbar(wfreqs,nanmean(freq_up_lag),nanstd(freq_up_lag)/sqrt(41),'r')

avg_trig_tfd_i = shiftdim(interp1(wfreqs,shiftdim(avg_trig8_tfd,1),nfreqs),2);
figure
imagesc(lags(used_lags)/Fsd,nfreqs,squeeze(mean(avg_trig_tfd_i)));shading flat
line([0 0],[5 80],'Color','w')
colorbar
