clear all
close all

load G:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010
addpath('G:\Code\WC_anal\general\')
sess_data(1:44) = [];
for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
    
    sync_time_jmm(); %synchronize MP and LFP timestamps (different sample frequencies)

    wcv_to_spktim(); %extract spike times from MP recordings 
    
    hip_wc_lfp_spk_shift_jmm(); %preprocess aligned MP and LFP signals.  Despike MP signal via linear interpolation
    
end

    
    