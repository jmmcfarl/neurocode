%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'G';

addpath(strcat(drive_letter,':\Code\smoothing\software'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
addpath(strcat(drive_letter,':\Code\maphmmbox\'))

cd(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
load parietal_cortical_2010
load desynch_times_individual
raw_Fs = 2016;

% get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

n = length(sess_data);
for d = 36:n
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load used_data lf8 wcv_minus_spike lf4

    %% for lf8
    [hsmm_bbstate_seq8_lfhf,hsmm_state_seq8_lfhf,hmm_bbstate_seq8_lfhf,hmm8_lfhf,Fs_bb,Fs] = ...
        parietal_get_hsmm_uds_state_seq_lfhf(lf8,raw_Fs,desynch_times_lf8{d},strcat(sess_data(d).name,'_lf8'));
    save hsmm_state_seq8_seg_lfhf hsmm* hmm* Fs* 
    clear hsmm* hmm* fract*
% 
%         %% for MP
%     [hsmm_bbstate_seq_lfhf,hsmm_state_seq_lfhf,hmm_bbstate_seq_lfhf,hmm_lfhf,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq_lfhf(wcv_minus_spike,raw_Fs,desynch_times_mp{d},strcat(sess_data(d).name,'_mp'));
%     save hsmm_state_seq_seg_lfhf hsmm* hmm* Fs* 
%     clear hsmm* hmm* fract*
    
    if sess_data(d).thom_elec
        %% for LF4
        [hsmm_bbstate_seq4_lfhf,hsmm_state_seq4_lfhf,hmm_bbstate_seq4_lfhf,hmm4_lfhf,Fs_bb,Fs] = ...
            parietal_get_hsmm_uds_state_seq_lfhf(lf4,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_lf4'));
        save hsmm_state_seq4_seg_lfhf hsmm* hmm* Fs*
        clear hsmm* hmm* fract*
    end
    
end

cd G:\WC_Germany\parietal_cortical_2010
