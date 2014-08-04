clear all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data

for d = 1:28
    d
    up_dur_rel_lfp_near{d} = zeros(length(synch_ups{d})-1,1);
    up_dur_rel_lfp_prev{d} = zeros(length(synch_ups{d})-1,1);
    up_dur_rel_lfp_near_cycle{d} = zeros(length(synch_ups{d})-1,1);
    bad_prev_ups = [];
    
    for i = 1:length(synch_ups{d})
        [dummy,near_lfp_up(i)] = min(abs(up_trans{d}(synch_ups{d}(i))- up_trans8{d}(synch_ups8{d}))); 
        temp = find(up_trans8{d}(synch_ups8{d}) < up_trans{d}(synch_ups{d}(i)),1,'last');
        up_dur_rel_lfp_near{d}(i) = up_state_dur{d}(i)/up_state_dur8{d}(synch_ups8{d}(near_lfp_up(i)));
        if synch_ups8{d}(near_lfp_up(i)) < length(down_state_dur8{d})
            up_dur_rel_lfp_near_cycle{d}(i) = up_state_dur{d}(i)/(up_state_dur8{d}(synch_ups8{d}(near_lfp_up(i)))+down_state_dur8{d}(synch_ups8{d}(near_lfp_up(i))));
        else
            up_dur_rel_lfp_near_cycle{d}(i) = nan;
        end
        if ~isempty(temp)
            prev_lfp_up(i) = temp;
            up_dur_rel_lfp_prev{d}(i) = up_state_dur{d}(i)/up_state_dur8{d}(synch_ups8{d}(prev_lfp_up(i)));
        else
            up_dur_rel_lfp_prev{d}(i) = nan;
            prev_lfp_up(i) = nan;
            bad_prev_ups = [bad_prev_ups i];
        end
    end

    mean_up_dur_rel_lfp_near(d) = nanmean(up_dur_rel_lfp_near{d});
    mean_up_dur_rel_lfp_prev(d) = nanmean(up_dur_rel_lfp_prev{d});
    fract_longer_up_near(d) = length(find(up_dur_rel_lfp_near{d}>1))/length(up_dur_rel_lfp_near{d});
    fract_longer_up_near_cycle(d) = length(find(up_dur_rel_lfp_near_cycle{d} > 1))/length(up_dur_rel_lfp_near_cycle{d});
    fract_longer_up_prev(d) = length(find(up_dur_rel_lfp_prev{d}>1))/length(up_dur_rel_lfp_prev{d});
    
    [dummy,p_val_near(d)] = ttest(up_state_dur{d}(synch_ups{d}),up_state_dur8{d}(near_lfp_up));
    synch_ups{d}(bad_prev_ups) = [];
    prev_lfp_up(bad_prev_ups) = [];
    [dummy,p_val_prev(d)] = ttest(up_state_dur{d}(synch_ups{d}),up_state_dur8{d}(prev_lfp_up));
    clear near_lfp_up prev_lfp_up
end

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
mec_cells(18) = []; %get rid of 04-07

