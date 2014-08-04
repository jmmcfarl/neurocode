clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds
all_mp_down_states = [];
all_lfp_down_states = [];

for d = 1:17
   
    mean_lfp_up(d) = nanmean(up_state_dur8{d}(synch_ups8{d}));
    mean_lfp_down(d) = nanmean(down_state_dur8{d}(synch_downs8{d}));
    mean_wcv_up(d) = nanmean(up_state_dur{d}(synch_ups{d}));
    mean_wcv_down(d) = nanmean(down_state_dur{d}(synch_downs{d}));
    median_wcv_up(d) = nanmedian(up_state_dur{d}(synch_ups{d}));
    mean_dc_lfp(d) = mean_lfp_up(d)/(mean_lfp_up(d)+mean_lfp_down(d));
    mean_dc_wcv(d) = mean_wcv_up(d)/(mean_wcv_up(d)+mean_wcv_down(d));
    all_mp_down_states = [all_mp_down_states down_state_dur{d}];
    all_lfp_down_states = [all_lfp_down_states down_state_dur8{d}];
end

save C:\WC_Germany\Persistent_activity\mean_state_dur_data mean* median*