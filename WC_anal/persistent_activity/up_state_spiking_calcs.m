clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

dsf = 8;
Fsd = 2016/dsf;

for d = 1:length(dir_array)
    cd(dir_array{d})
pwd
     if ~exist('spike_time.mat')
        load spike_time_br
    else
        load spike_time
     end
load used_data lf8

    spike_ids = round(spkid/dsf);
     
    up_rate{d} = zeros(size(synch_ups{d}));
    for i = 1:length(synch_ups{d})
        cur_spikes = find(spike_ids > up_trans{d}(synch_ups{d}(i)) & ...
            spike_ids < down_trans{d}(synch_ups{d}(i)));
        up_rate{d}(i) = length(cur_spikes)/up_state_dur{d}(synch_ups{d}(i));
    end
    
    mean_up_rate(d) = nanmean(up_rate{d});
    std_up_rate(d) = nanstd(up_rate{d});
    mean_rate(d) = length(spike_ids)/length(lf8)*2016;
    
    clear lf8
    
    
end

