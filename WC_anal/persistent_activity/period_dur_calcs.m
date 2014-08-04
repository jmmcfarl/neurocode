clear all

load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data 
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

for i = 1:17
    
    synch_states = intersect(synch_ups{i},synch_downs{i});
    
    period_durs{i} = up_state_dur{i}(synch_states)+down_state_dur{i}(synch_states);
    
    mean_period_dur(i) = nanmean(period_durs{i});
    med_period_dur(i) = nanmedian(period_durs{i});
    
    mean_down_dur(i) = nanmean(down_state_dur{i}(synch_states));
    med_down_dur(i) = nanmedian(down_state_dur{i}(synch_states));
    
end


%%

load('C:\WC_Germany\Persistent_activity\persistent_analysis\quant_data_12_05.mat')

for i = 1:17
    
    good_ups = cor_lfp_period_dur{i}(cor_lfp_period_dur{i} > 0);
    
    mean_up_lp(i) = nanmean(good_ups);
    med_up_lp(i) = nanmedian(good_ups);
    
end