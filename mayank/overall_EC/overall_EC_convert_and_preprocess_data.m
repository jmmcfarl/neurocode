clear all
close all

load G:\WC_Germany\overall_EC\overall_allcells_dir.mat
addpath('G:\Code\WC_anal\general\')

for d = 121:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
    
    sync_time_jmm();

    wcv_to_spktim();
    
    hip_wc_lfp_spk_shift_jmm(sess_data(d).gains);
    
end

    
    