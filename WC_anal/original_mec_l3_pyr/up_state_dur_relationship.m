clear all

load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28
Fsd = 2016/8;
for d = 1:length(dir_array)
    
    prev_up_dur{d} = zeros(1,length(up_trans{d}));
    time_since_prev{d} = prev_up_dur{d};
   for i = 1:length(up_trans{d})
       
       prev_lfp_up = find(up_trans8{d} < up_trans{d}(i),1,'last');
       if ~isempty(prev_lfp_up)
           prev_up_dur{d}(i) = up_state_dur8{d}(prev_lfp_up);
           time_since_prev{d}(i) = (up_trans{d}(i)-up_trans8{d}(prev_lfp_up))/Fsd;
       else
           prev_up_dur{d}(i) = nan;
       end
       
       
   end
   prev_up_dur{d}(prev_up_dur{d}>5) = nan;
   [a,b] = corrcoef(prev_up_dur{d}(~isnan(prev_up_dur{d})),up_state_dur{d}(~isnan(prev_up_dur{d})));
   
   up_dur_cor(d) = a(2,1);
   up_dur_cor_p(d) = b(2,1);
   
    max_lfp = max(prev_up_dur{d});
    plot(prev_up_dur{d},up_state_dur{d},'o')
    line([0 max_lfp],[0 max_lfp],'Color','k')
    xlabel('Previous LFP Up State Duration','FontSize',14)
    ylabel('WCV Up State Duration','FontSize',14)
    title(sprintf('%0.3g  %0.3g',up_dur_cor(d),up_dur_cor_p(d)))
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\Up_Duration_cor\' f_names{d}]
    print('-dpng',t_names)
    close all
    
%     plot(time_since_prev{d},up_state_dur{d},'o')
%     
%     plot(up_trans{d}/Fsd,up_state_dur{d},'o')
%     hold on
%     plot(up_trans{d}/Fsd,prev_up_dur{d},'ro')
end