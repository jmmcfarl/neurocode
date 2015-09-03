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
addpath(strcat(drive_letter,':\WC_Germany\persistent_9_27_2010\'))
addpath(strcat(drive_letter,':\Code\maphmmbox\'))

raw_Fs = 2016;

load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\WC_Germany\parietal_cortical_2010\')
used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

n = length(sess_data);
for d = 1:n
% d=14;
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load ./used_data wcv_minus_spike lf8
    load ./desynch_times_lf8
    
%% compute MP state sequences
    tic
    [hsmm_bbstate_seq,hsmm_state_seq,hmm_bbstate_seq,hsmm,hmm,Fs_bb,Fs] = ...
        pa_get_hsmm_uds_state_seq_new(wcv_minus_spike,raw_Fs,desynch_times_lf8,sess_data(d).name);
    hsmm.run_dur = toc;
    bad_models = 0;
    for i = 1:hsmm.Nsegs
        curT = length(hsmm.UDS_segs(i,1):hsmm.UDS_segs(i,2));
        fract_low_post{i}(d) = sum(hsmm.posteriors{i}(:) < 0.05)/curT;
       bad_models = bad_models + length(hsmm.rob_model_inds{i});
    end
    fract_bad_model(d) = bad_models/hsmm.Fs/hsmm.uds_dur;
    save pa_hsmm_state_seq_new3 hsmm* hmm* Fs* fract_low_post*
    clear hsmm* hmm* fract*
% load pa_hsmm_state_seq_new2
% hsmm_bbstate_seq = thresh_state_smooth_seg(hsmm_bbstate_seq,252,300,300); %smooth out excessively short state durations to improve model fits
% save pa_hsmm_state_seq_sm300

%% compute LF8 state sequences
        disp('computing LF8 state sequences')
        tic
        [hsmm_bbstate_seq8,hsmm_state_seq8,hmm_bbstate_seq8,hsmm8,hmm8,Fs_bb,Fs] = ...
            pa_get_hsmm_uds_state_seq_new(lf8,raw_Fs,desynch_times_lf8,strcat(sess_data(d).name,'_lf8'));
        hsmm8.run_dur8 = toc;
        bad_models = 0;
        for i = 1:hsmm8.Nsegs
            curT = length(hsmm8.UDS_segs(i,1):hsmm8.UDS_segs(i,2));
            bad_models = bad_models + length(hsmm8.rob_model_inds{i});
            fract_low_post8{i}(d) = sum(hsmm8.posteriors{i}(:) < 0.05)/curT;
        end
        save pa_hsmm_state_seq8_new2 hsmm* hmm* Fs* fract_low_post*
        clear hsmm* hmm* fract*
% load pa_hsmm_state_seq8_new2
% hsmm_bbstate_seq8 = thresh_state_smooth_seg(hsmm_bbstate_seq8,252,300,300); %smooth out excessively short state durations to improve model fits
% save pa_hsmm_state_seq8_sm300

end
cd G:\WC_Germany\persistent_9_27_2010
