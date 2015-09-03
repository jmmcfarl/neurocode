clear all
close all

load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hsmm_state_detection\')
addpath('F:\WC_Germany\persistent_revised\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

cd F:\WC_Germany\persistent_9_27_2010\
load ./pa_corresponding_lfp_state_data
load ./pa_heka_UDS_data
load ./spike_rate_data

for d = 1:length(sess_data)
  
    plot(mp_upstate_heka_meanamp_sig{d},mp_downstate_heka_meanamp_sig{d},'.')
    pause
    clf
    
    
end

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));

