%% look at correlation between up state and down staet durations
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
maxlag = 300;

for d = 1:17
    
%     p = corrcoef(avg_up_dur{d},avg_up_dur8{d});
%     up_cor(d) = p(2,1);
%     
%     figure
%     plot(avg_up_dur{d},avg_up_dur8{d},'o')
%     title(num2str(up_cor(d)))
%     tname = ['C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_dur_scatter_' f_names{d}];
%     print('-dpng',tname)
%     close all
%     
%     p = corrcoef(max_up_dur{d},max_up_dur8{d});
%     up_max_cor(d) = p(2,1);
% 
%     figure
%         plot(max_up_dur{d},max_up_dur8{d},'o')
%     title(num2str(up_max_cor(d)))
%     tname = ['C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_max_dur_scatter_' f_names{d}];
%     print('-dpng',tname)
% close all


[up_dur_acorr(d,:),lags] = xcov(avg_up_dur{d},maxlag,'coeff');
[up_dur8_acorr(d,:),lags] = xcov(avg_up_dur8{d},maxlag,'coeff');
[up_dur_xcor(d,:),lags] = xcov(avg_up_dur{d},avg_up_dur8{d},maxlag,'coeff');

    figure
    subplot(2,1,1)
plot(lags,up_dur_acorr(d,:))
hold on
plot(lags,up_dur8_acorr(d,:),'r')
subplot(2,1,2)
plot(lags,up_dur_xcor(d,:))
    tname = ['C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_dur_cor_' f_names{d}];
    print('-dpng',tname)


    
    
end
    