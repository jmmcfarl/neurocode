clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir

Fs = 2016;

for d = 1:length(sess_data)
    
    cd(sess_data(d).directory)
    pwd
    
    load spike_time_jmm
    
    load used_data lf8
    
    total_dur = length(lf8)/Fs;
    
    mean_rate(d) = length(spkid)/total_dur;
    
end

cd C:\WC_Germany\Cortical_analysis
save cortical_mean_rate_data mean_rate