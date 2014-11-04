function [] = process_wc_data(dir_name)

cd(dir_name)

% sync_time_jmm;
% clear all
wcv_to_spktim;
clear all
hip_wc_lfp_spk_shift_combined();
