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
Fs_desired = 50;
Fs_new = 2016/8;

l23 = 1:n_23_pyr_par;
l56 = n_23_pyr_par + n_int_par+1:n_23_pyr_par + n_int_par+n_56_pyr_par;
used_cells = [l23 l56];

for d = 1:length(used_cells)
  d  
    cd(sess_data(used_cells(d)).directory)
    load used_data lf8 wcv_minus_spike
    load hsmm_state_seq8_lf_pert15
    
    [lf_features,t_axis,Fs] = get_lf_features(lf8,raw_Fs,Fs_desired);
    [lf_featuresw,t_axis,Fs] = get_lf_features(wcv_minus_spike,raw_Fs,Fs_desired);
    
%     [gamma] = get_hsmm_gamma(hmm8,lf_features);
    gamma = hmm8.gamma;

%     [scal_features,wfreqs,t_axis_new,Fs_hf] = get_scalogram_features(lf8,raw_Fs,Fs_new);
    [hf_features,t_axis_new,Fs_hf] = get_hf_features(lf8,raw_Fs,Fs_new);
    
    [comb_hsmm_state_seq,no_peak_ups,no_peak_downs] = pert_optimize_hsmm_transitions...
        (hf_features,hmm8,hsmm_state_seq8,gamma,Fs,Fs_new,'fixed',[1 0]);
    Fs_comb = Fs_hf;
    
%     plot(t_axis,lf_features), hold on
%     plot(t_axis_new,comb_hsmm_state_seq,'r','linewidth',2)
%     plot(t_axis,hsmm_state_seq8-2,'k','linewidth',2)
%     plot(t_axis_new,hf_features,'g')
% plot(t_axis_new,hsmm_bbstate_seq-1,'c','linewidth',2)
% pause
    
    save hsmm_state_seq8_comb_sm4 comb_hsmm* Fs_comb no_peak*
    
end

