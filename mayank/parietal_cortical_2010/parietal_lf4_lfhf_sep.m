%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'G';

addpath('G:\WC_Germany\hmm_state_detect\\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010')
% addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\persistent_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load desynch_times_individual

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
Fs = 50.4;

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);

min_state_dur = 1;
max_state_dur = round(Fs*30);
dur_range = (1:max_state_dur)/Fs;

dur_hist_bincenters = (min_state_dur:max_state_dur)/Fs;

hcf_vals = [0.5 1 2 4 10];
smooth_vals = [0.025 0.05 0.1 0.2 0.3];

short_threshold = 0.2;

n = length(sess_data);

%% For LF4
for d = 1:length(thom_el)
    d
    direct = sess_data(thom_el(d)).directory;
    direct(1) = 'G';
    cd(direct)
    load used_data lf4
    [hmm_bbstate_seq4_lfhf,hmm4_lfhf] = parietal_get_hmm_uds_state_seq_lfhf(lf4,raw_Fs,desynch_times_lf4{thom_el(d)});
    meandiff1 = [];
    meandiff2 = [];
    for i = 1:hmm4_lfhf.Nsegs
        meandiff1 = [meandiff1; hmm4_lfhf.state(2).meanfun{i}-hmm4_lfhf.state(1).meanfun{i}];
        meandiff2 = [meandiff2; hmm4_lfhf.state(2).meanfun2{i}-hmm4_lfhf.state(1).meanfun2{i}];
    end
    meandiff = [meandiff1 meandiff2];
    covar1 = hmm4_lfhf.state(1).var;
    covar2 = hmm4_lfhf.state(2).var;
    kl1 = gauss_kl_div(meandiff',covar1,covar2);
    kl2 = gauss_kl_div(-meandiff',covar2,covar1);
    lf4_kl_lfhf(d) = kl1+kl2;
    
    save hmm_lf4_lfhf hmm*
    clear hmm*
end
cd G:\WC_Germany\parietal_cortical_2010
save hmm_lf4_lfhf_kl lf4_kl_lfhf
