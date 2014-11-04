clear all
close all

addpath(genpath('C:\Code\'))
addpath(genpath('C:\WC_Germany\'))

cd C:/WC_Germany/persistent_downs/
load ./new_pdown_dir
true_heka_fs = 1.99984e+004;

uset = 36:length(new_pdown_dir);
% uset = 31
% uset(12) = [];
for d = uset
    d
    cd(new_pdown_dir{d})
    pwd   
    
    load ./used_data wcv
    load ./sync_times
    load ./heka_data
    
    dc_times = (1:length(dc_data))/true_heka_fs;
    wcv_int = wcv;
    [dc_offset(d),dc_maxcorr(d)] = align_dc_ac_sigs_initial_v2(dc_data,dc_times,wcv_int);
    dc_time = dc_times + dc_offset(d);
    ac_time = (1:length(wcv_int))/2016;
    if min(diff(dc_time)) <= 0
        error('Alignment Problem!')
    end
    save aligned_heka dc_data dc_time ac_time
    
    clear dc_data 
end