clear all
close all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data

maxdur = 50;
up_range = [0.3 maxdur];
down_range = [0.3 maxdur];
numBins = 50;

for d = 1:length(dir_array)
    
    avg_up_dur(d) = nanmean(synch_up_dur{d});
    avg_up_dur8(d) = nanmean(synch_up_dur8{d});
    avg_down_dur(d) = nanmean(synch_down_dur{d});
    avg_down_dur8(d) = nanmean(synch_down_dur8{d});
    
    [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d}(synch_up_dur{d} <= maxdur),up_range,numBins);
    [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d}(synch_up_dur8{d} <= maxdur),up_range,numBins);
        
    [down_hist(d,:),down_grid] = log_hist(synch_down_dur{d}(synch_down_dur{d} <= maxdur),down_range,numBins);
    [down_hist8(d,:),down_grid] = log_hist(synch_down_dur8{d}(synch_down_dur8{d} <= maxdur),down_range,numBins);
      
end

save C:\WC_Germany\persistent_revised\state_dur_hist_data


