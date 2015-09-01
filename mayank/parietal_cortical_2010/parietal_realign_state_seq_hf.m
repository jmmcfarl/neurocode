%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'G';

addpath(strcat(drive_letter,':\Code\smoothing\software'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
% addpath(strcat(drive_letter,':\WC_Germany\new_stellate_analysis\'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))

cd(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
load parietal_cortical_2010
load desynch_times_mp_lf8_2_24_09
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
high_lcf = 15;
high_hcf = 40;
hf_smooth = 0.03;
up_max_perthf = [0.4 0];
down_max_perthf =[0.4 0.4];
%get rid of interneurons
% interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
% sess_data(interneurons) = [];

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);

% thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
% sess_data = sess_data(thomas_el);

for d = 1:length(sess_data)
    d
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load used_data lf8 wcv_minus_spike lf4

    %     lf5_hf = get_hf_features(lf5,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    lf8_hf = get_hf_features(lf8,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    wcv_hf = get_hf_features(wcv_minus_spike,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    if sess_data(d).thom_elec
        lf4_hf = get_hf_features(lf4,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    end

%     load  hsmm_state_seq_seg_lf_pert15 
%     [hsmm_hfstate_seq] = pert_optimize_hsmm_transitions_seg...
%         (wcv_hf,hmm,hsmm_state_seq,hmm.gamma,hmm.Fs,Fsd,max_perthf);
%     save hsmm_state_seq_seg_lf_pert15 hsmm* hmm* Fs* fract_low_post*
%     clear hsmm* hmm* fract*

    load  hsmm_state_seq8_seg_lf_pert15 
    [hsmm_hfstate_seq8] = pert_optimize_hsmm_transitions_seg...
        (lf8_hf,hmm8,hsmm_state_seq8,hmm8.gamma,hmm8.Fs,Fsd,up_max_perthf,down_max_perthf);
    save hsmm_state_seq8_seg_lf_pert15 hsmm* hmm* Fs* fract_low_post*
    clear hsmm* hmm* fract*

        load  hsmm_state_seq_seg_lf_pert15 
    [hsmm_hfstate_seq] = pert_optimize_hsmm_transitions_seg...
        (wcv_hf,hmm,hsmm_state_seq,hmm.gamma,hmm.Fs,Fsd,up_max_perthf,down_max_perthf);
    save hsmm_state_seq_seg_lf_pert15 hsmm* hmm* Fs* fract_low_post*
    clear hsmm* hmm* fract*

    if sess_data(d).thom_elec
        load  hsmm_state_seq4_seg_lf_pert15
        [hsmm_hfstate_seq4] = pert_optimize_hsmm_transitions_seg...
            (lf4_hf,hmm4,hsmm_state_seq4,hmm4.gamma,hmm4.Fs,Fsd,up_max_perthf,down_max_perthf);
        save hsmm_state_seq_seg_lf_pert15 hsmm* hmm* Fs* fract_low_post*
        clear hsmm* hmm* fract*
    end
    
end

