%% analyze up slopes
clear all
load C:\WC_Germany\JMM_analysis_pyr\WCV_LFP_up_trans_sig_fit_data

tau_range = linspace(0,log(2e5),50);

for d = 1:17
    
   rat_good_ups_wcv(d) = length(~isnan(rltau_wup{d}))/length(rltau_wup{d});
   rat_good_ups_lf8(d) = length(~isnan(rltau_lup{d}))/length(rltau_lup{d});
   
   wcv_up_dist(d,:) = hist(log(rltau_wup{d}),tau_range);
   lfp_up_dist(d,:) = hist(log(rltau_lup{d}),tau_range);
   wcv_up_dist(d,:) = wcv_up_dist(d,:)/sum(wcv_up_dist(d,:));
   lfp_up_dist(d,:) = lfp_up_dist(d,:)/sum(lfp_up_dist(d,:));
   
end