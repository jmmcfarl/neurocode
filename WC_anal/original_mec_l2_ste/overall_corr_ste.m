%%calculate overall correlograms and up triggered MP averages

clear all

% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\correl_data
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\lag_data
load C:\WC_Germany\JMM_Analysis_ste\dir_tree_ste

dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
hcf = 40/niqf;
lcf = 0.05/niqf;
[b,a] = butter(2,[lcf hcf]);
maxlag = round(10*Fsd);
tlags = [-maxlag:maxlag]/Fsd;

for d = 1:length(dir_array)
    
       disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd
    
    load used_data wcv_minus_spike lf8
 
    wcv = filtfilt(b,a,wcv_minus_spike);
    wcv_d = downsample(wcv,dsf);
    wcv_d = zscore(wcv_d);
    
    lf8 = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8,dsf);
    lf8_d = zscore(lf8_d);
        
    
    %calculate overall correlograms
    [l_acorr(d,:),lags] = xcov(lf8_d,maxlag,'coeff');
    [w_acorr(d,:),lags] = xcov(wcv_d,maxlag,'coeff');
    [x_corr(d,:),lags] = xcov(wcv_d,lf8_d,maxlag,'coeff');
    
    
    %%create figures
       Fig = figure(1)
       clf
       set(Fig,'PaperUnits','centimeters');
       set(gcf, 'PaperSize', [20 40]);% paper size is in [width height] format
        subplot(3,1,1)
        plot(lags/Fsd,x_corr(d,:),'linewidth',2)
        xlim([-5 5])
        title('Cross Correlogram')
        grid
        subplot(3,1,2)
        plot(lags/Fsd,w_acorr(d,:),'linewidth',2)
        title('WCV autocorrelogram')
        grid
        xlim([-5 5])
        subplot(3,1,3)
        plot(lags/Fsd,l_acorr(d,:),'linewidth',2)
        grid
        xlim([-5 5])
        title('LFP autocorrelogram')
        tname = ['C:\WC_Germany\JMM_Analysis_ste\time_domain_overall\' f_names{d}];
        print('-dpng',tname);
        close all
        
            
end



save C:\WC_Germany\JMM_Analysis_ste\overall_td_anal l_acorr w_acorr x_corr lags


% % %create overall plots
% 
m_l_acorr = mean(l_acorr);
u_l_acorr = m_l_acorr+2*std(l_acorr)/sqrt(17);
l_l_acorr = m_l_acorr-2*std(l_acorr)/sqrt(17);

m_w_acorr = mean(w_acorr);
u_w_acorr = m_w_acorr+2*std(w_acorr)/sqrt(17);
l_w_acorr = m_w_acorr-2*std(w_acorr)/sqrt(17);

m_x_corr = mean(x_corr);
u_x_corr = m_x_corr+2*std(x_corr)/sqrt(17);
l_x_corr = m_x_corr-2*std(x_corr)/sqrt(17);


%plot l_acorr
plot(lags/Fsd,m_l_acorr,'linewidth',2)
hold on
plot(lags/Fsd,u_l_acorr,'--')
plot(lags/Fsd,l_l_acorr,'--')
xlim([0 10])
ylim([-0.4 0.4])
xlabel('Time (s)','FontSize',14)
ylabel('Correlation Coefficient','FontSize',14)
grid

%plot w_acorr
plot(lags/Fsd,m_w_acorr,'linewidth',2)
hold on
plot(lags/Fsd,u_w_acorr,'--')
plot(lags/Fsd,l_w_acorr,'--')
xlim([0 10])
ylim([-0.3 0.3])
xlabel('Time (s)','FontSize',14)
ylabel('Correlation Coefficient','FontSize',14)
grid

%plot xcorr
plot(lags/Fsd,m_x_corr,'linewidth',2)
hold on
plot(lags/Fsd,u_x_corr,'--')
plot(lags/Fsd,l_x_corr,'--')
xlim([-5 5])
line([0 0],[-0.4 0.4],'Color','k')
ylim([-0.4 0.8])
xlabel('Time (s)','FontSize',14)
ylabel('Correlation Coefficient','FontSize',14)
grid

