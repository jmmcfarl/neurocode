clear all

load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update

dsf =8;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:17
    
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
    
end
    
save C:\WC_Germany\JMM_Analysis_pyr\overall_hist_data lf8_dist wcv_dist gridv

m_8 = mean(lf8_dist);
u_8 = m_8+2*std(lf8_dist)/sqrt(17);
l_8 = m_8-2*std(lf8_dist)/sqrt(17);

m_w = mean(wcv_dist);
u_w = m_w+2*std(wcv_dist)/sqrt(17);
l_w = m_w-2*std(wcv_dist)/sqrt(17);

   
plot(gridv,m_8,'r','linewidth',2)
hold on
plot(gridv,u_8,'r--')
plot(gridv,l_8,'r--')

plot(gridv,m_w,'linewidth',2)
plot(gridv,u_w,'--')
plot(gridv,l_w,'--')
xlim([-3.5 3.5])
xlabel('Amplitude (z-score)','FontSize',14)
ylabel('Probability','FontSize',14)