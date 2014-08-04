clear all

load C:\WC_Germany\JMM_Analysis_ste\dir_tree_ste

dsf = 8;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:length(dir_array)
    
    cd(dir_array{d})
    pwd
    
    load used_data lf8 wcv_minus_spike
    
    lf8 = filtfilt(b,a,lf8);
    wcv = filtfilt(b,a,wcv_minus_spike);
    lf8_d = downsample(lf8,dsf);
    wcv_d = downsample(wcv,dsf);
    
    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);
    
    [lf8_dist(d,:),gridv] = gpkde(lf8_d,-3,[-4 4]);
    [wcv_dist(d,:),gridv] = gpkde(wcv_d,-3,[-4 4]);
    
    
        
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    plot(gridv,lf8_dist(d,:),'r','linewidth',2)
    hold on
    plot(gridv,wcv_dist(d,:),'linewidth',2)
    legend('LFP','WCV')
    tnames = ['C:\WC_Germany\JMM_Analysis_ste\hist\' f_names{d}];
    print(tnames,'-dpng')
    close all
    
    
end
    
save C:\WC_Germany\JMM_Analysis_ste\overall_hist_data lf8_dist wcv_dist gridv


% %% Create overall plot
m_8 = mean(lf8_dist);
u_8 = m_8+2*std(lf8_dist)/sqrt(17);
l_8 = m_8-2*std(lf8_dist)/sqrt(17);

m_w = mean(wcv_dist);
u_w = m_w+2*std(wcv_dist)/sqrt(17);
l_w = m_w-2*std(wcv_dist)/sqrt(17);
% 
%    
plot(gridv,m_8,'r','linewidth',2)
hold on
plot(gridv,m_w,'linewidth',2)
legend('LFP','WCV')
plot(gridv,u_8,'r--')
plot(gridv,l_8,'r--')
plot(gridv,u_w,'--')
plot(gridv,l_w,'--')
xlim([-2.5 4])
xlabel('Amplitude (z-score)','FontSize',14)
ylabel('Probability','FontSize',14)