%%calculate overall correlograms and up triggered MP averages

clear all

% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\correl_data
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\lag_data
load C:\WC_Germany\JMM_Analysis_ste\dir_tree_ste
load C:\WC_Germany\JMM_Analysis_ste\UDS_dur_run_hist\data
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
    
    up_trans8{d}(up_trans8{d}<maxlag) = [];
    up_trans8{d}(up_trans8{d}>length(wcv_d)-maxlag) = [];
    
    num_8ups = length(up_trans8{d});
    
    up_trans_mat = zeros(num_8ups,2*maxlag+1);
    up_trans_mat8 = up_trans_mat;
    
    for t = 1:num_8ups
       
        up_trans_mat(t,:) = wcv_d(up_trans8{d}(t)-maxlag:up_trans8{d}(t)+maxlag);
        up_trans_mat8(t,:) = lf8_d(up_trans8{d}(t)-maxlag:up_trans8{d}(t)+maxlag);
    end
    
    up_trig_wcv(d,:) = mean(up_trans_mat);
    up_trig_lfp(d,:) = mean(up_trans_mat8);
    
    down_trans8{d}(down_trans8{d}<maxlag) = [];
    down_trans8{d}(down_trans8{d}>length(wcv_d)-maxlag) = [];
    
    num_8downs = length(down_trans8{d});
    
    down_trans_mat = zeros(num_8downs,2*maxlag+1);
    down_trans_mat8 = down_trans_mat;
    
    for t = 1:num_8downs
       
        down_trans_mat(t,:) = wcv_d(down_trans8{d}(t)-maxlag:down_trans8{d}(t)+maxlag);
        down_trans_mat8(t,:) = lf8_d(down_trans8{d}(t)-maxlag:down_trans8{d}(t)+maxlag);
        
    end
    
    down_trig_wcv(d,:) = mean(down_trans_mat);
    down_trig_lfp(d,:) = mean(down_trans_mat8);
    
    
       
    
%     %%create figures
%        Fig = figure(1)
%        clf
%        set(Fig,'PaperUnits','centimeters');
%        set(gcf, 'PaperSize', [30 30]);% paper size is in [width height] format
%         subplot(2,2,1)
%         plot(lags/Fsd,x_corr(d,:),'linewidth',2)
%         title('Cross Correlogram')
%         subplot(2,2,2)
%         plot(lags/Fsd,w_acorr(d,:),'linewidth',2)
%         title('WCV autocorrelogram')
%         subplot(2,2,3)
%         plot(lags/Fsd,l_acorr(d,:),'linewidth',2)
%         title('LFP autocorrelogram')
%         subplot(2,2,4)
%         plot(tlags,up_trig_wcv(d,:),'linewidth',2)
%         title('LFP Up Triggered MP')
%         tname = ['C:\WC_Germany\JMM_Analysis_pyr\time_domain_overall\long_per' f_names{d}];
%         print('-dpng',tname);
%         close all

%     %%create figures
%        Fig = figure(1)
%        clf
%        set(Fig,'PaperUnits','centimeters');
%        set(gcf, 'PaperSize', [30 30]);% paper size is in [width height] format
%         plot(tlags,down_trig_wcv(d,:),'linewidth',2)
%         tname = ['C:\WC_Germany\JMM_Analysis_pyr\time_domain_overall\down_trig' f_names{d}];
%         print('-dpng',tname);

        
        
            %%create figures
       Fig = figure(1)
       clf
       set(Fig,'PaperUnits','centimeters');
       set(gcf, 'PaperSize', [30 30]);% paper size is in [width height] format
        subplot(2,1,1)
        plot(tlags,up_trig_wcv(d,:),'linewidth',2)
        hold on
        plot(tlags,up_trig_lfp(d,:),'r','linewidth',2)
        title('Up Trig')
        xlim([-5 5])
        subplot(2,1,2)
        plot(tlags,down_trig_wcv(d,:),'linewidth',2)
        hold on
        plot(tlags,down_trig_lfp(d,:),'r','linewidth',2)
        title('Down Trig')
        xlim([-5 5])
        tname = ['C:\WC_Germany\JMM_Analysis_ste\uta_dta\' f_names{d}];
        print('-dpng',tname);
        close all

        
                    %%create figures
       Fig = figure(1)
       clf
       set(Fig,'PaperUnits','centimeters');
       set(gcf, 'PaperSize', [30 30]);% paper size is in [width height] format
        subplot(2,1,1)
        plot(tlags,up_trig_wcv(d,:),'linewidth',2)
        hold on
        plot(tlags,up_trig_lfp(d,:),'r','linewidth',2)
        title('Up Trig')
        xlim([-1 1]);grid
        subplot(2,1,2)
        plot(tlags,down_trig_wcv(d,:),'linewidth',2)
        hold on
        plot(tlags,down_trig_lfp(d,:),'r','linewidth',2)
        title('Down Trig')
        xlim([-1 2]);grid
        tname = ['C:\WC_Germany\JMM_Analysis_ste\uta_dta\zoom_' f_names{d}];
        print('-dpng',tname);
        close all

        
    clear up_trans_mat down_trans_mat
    
    
end

save C:\WC_Germany\JMM_Analysis_ste\Uta_Dta up_trig* down_trig* tlags


% %% Plot figures
% 
m_uta = mean(up_trig_wcv);
u_uta = m_uta+2*std(up_trig_wcv)/sqrt(17);
l_uta = m_uta-2*std(up_trig_wcv)/sqrt(17);
% 

m_dta = mean(down_trig_wcv);
u_dta = m_dta+2*std(down_trig_wcv)/sqrt(17);
l_dta = m_dta-2*std(down_trig_wcv)/sqrt(17);

m_uta8 = mean(up_trig_lfp);
u_uta8 = m_uta8+2*std(up_trig_lfp)/sqrt(17);
l_uta8 = m_uta8-2*std(up_trig_lfp)/sqrt(17);
% 

m_dta8 = mean(down_trig_lfp);
u_dta8 = m_dta8+2*std(down_trig_lfp)/sqrt(17);
l_dta8 = m_dta8-2*std(down_trig_lfp)/sqrt(17);
% 
% 
% 
%plot uta
plot(tlags,m_uta,'linewidth',2)
hold on
plot(tlags,u_uta,'--')
plot(tlags,l_uta,'--')
plot(tlags,m_uta8,'r','linewidth',2)
plot(tlags,u_uta8,'r--')
plot(tlags,l_uta8,'r--')
xlim([-5 5])
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (z-score)','FontSize',14)
grid

% xlim([-5 5])
% line([0.85 0.85],[-0.8 0.8],'Color','r')
% line([0 0],[-0.8 0.8],'Color','k')
% line([-1.175 -1.175],[-0.8 0.8],'Color','r')

%plot dta
plot(tlags,m_dta,'linewidth',2)
hold on
plot(tlags,u_dta,'--')
plot(tlags,l_dta,'--')
plot(tlags,m_dta8,'r','linewidth',2)
hold on
plot(tlags,u_dta8,'r--')
plot(tlags,l_dta8,'r--')
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (z-score)','FontSize',14)

xlim([-5 5])
grid
