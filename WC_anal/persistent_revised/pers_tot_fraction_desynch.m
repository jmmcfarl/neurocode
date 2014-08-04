clear all
close all
cd C:\WC_Germany\persistent_revised\desynch_detect
load desynch_times
load('C:\WC_Germany\persistent_revised\pers_revised_dir.mat')
Fs = 2016;
for d = 1:28
    cd(dir_array{d})
    load used_data lf8
    num_desynch_times(d) = length(desynch_start_times{d});
    
    tot_desynch_time(d) = 0;
    for i = 1:length(desynch_start_times{d})
        tot_desynch_time(d) = tot_desynch_time(d) + (desynch_stop_times{d}(i)-desynch_start_times{d}(i));
    end
    
    tot_time(d) = length(lf8)/Fs;
  
    
end
cd C:\WC_Germany\persistent_revised\
