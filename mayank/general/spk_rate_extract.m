load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat

for d = 1:length(dir_array)
    
    cd(dir_array{d})
    pwd
    
    if exist('spike_time.mat')
        load spike_time   
    else 
        load spike_time_br
    end
    
    load sync_times synct
    
    num_spikes = length(spkid);
    record_dur = max(synct)/1e6;
    
    av_rate(d) = num_spikes/record_dur;
    
    
    clear num_spikes record_dur synct spkid spkamp
    
end

save E:\WC_Germany\jmm_Analysis\av_rate av_rate