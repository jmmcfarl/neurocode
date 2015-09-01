%inputs:
% raw_data
% raw_Fs
clear all
close all

addpath('G:\Code\smoothing\software')
addpath('G:\Code\FullBNT-1.0.4\KPMstats\')
addpath('G:\Code\FullBNT-1.0.4\netlab3.3')
addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010

raw_Fs = 2016;
Fs_desired = raw_Fs/8;

for d = 1:length(sess_data)
    
    cd(sess_data(d).directory)
    load used_data lf8 wcv_minus_spike
    
    [np_state_seq,t_axis,Fs] = get_np_state_sequence(wcv_minus_spike,raw_Fs);
    save np_state_seq
   
    
end

