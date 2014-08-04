clear all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data

dsf = 8;
Fsd = 2016/dsf;

for d = 1:28
    cd(dir_array{d})
pwd
load spike_time_jmm
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

lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
mec_cells = 1:21;
mec_cells(21) = [];


lec_up_rate = mean_up_rate(lec_cells);
lec_rate = mean_rate(lec_cells);
mec_up_rate = mean_up_rate(mec_cells);
mec_rate = mean_rate(mec_cells);


