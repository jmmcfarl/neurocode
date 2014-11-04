clear all
close all
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd_dist.mat
%%
for d = 89:length(combined_dir)
    cd(combined_dir{d})
    pwd
    
    if exist('./Sc1.ntt','file')
        amp_threshold = 30;
        max_overlap = 0.5;
        [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
        
        save mua_data3 mua_times mua_amps mua_widths avg_waveform std_waveform
    end
end