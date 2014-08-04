clear all

load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update

niqf = 2016/2;
hcf = 40/niqf;
lcf = 0.05/niqf;
[b,a] = butter(2,[lcf hcf]);
dsf = 8;
Fsd = 2016/dsf;
maxlag = 10*Fsd;

for d = 1:length(dir_array)
    
    cd(dir_array{d})
    pwd
    
    load used_data lf8 wcv_minus_spike
    
    lf8_f = filtfilt(b,a,lf8);
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    
    lf8_d = downsample(lf8_f,dsf);
    wcv_d = downsample(wcv_f,dsf);
    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);
    
    [w_acorr,lags] = xcov(wcv_d,maxlag,'coeff');
    [l_acorr,lags] = xcov(lf8_d,maxlag,'coeff');
    [xcor,lags] = xcov(wcv_d,lf8_d,maxlag,'coeff');
    
    plot(lags/Fsd,w_acorr,'linewidth',2)
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\time_domain_overall\wacorr\z_' f_names{d}];
    ylim([-0.3 0.3])
    print('-dpng',t_names)
    close all
    
    plot(lags/Fsd,l_acorr,'linewidth',2)
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\time_domain_overall\lacorr\z_' f_names{d}];
    ylim([-0.5 0.5])
    print('-dpng',t_names)
    close all

    plot(lags/Fsd,xcor,'linewidth',2)
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\time_domain_overall\xcorr\z_' f_names{d}];
    print('-dpng',t_names)
    close all
    
    
   
end