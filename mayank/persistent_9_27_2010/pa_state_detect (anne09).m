%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'F';

addpath(strcat(drive_letter,':\Code\smoothing\software'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
addpath(strcat(drive_letter,':\WC_Germany\persistent_9_27_2010\'))
addpath(strcat(drive_letter,':\Code\maphmmbox\'))

raw_Fs = 2016;

load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\WC_Germany\parietal_cortical_2010\')
used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

n = length(sess_data);
for d = 7:10
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load ./used_data wcv_minus_spike lf8
    load ./desynch_times_lf8
    
%% compute MP state sequences
    tic
    [hsmm_bbstate_seq,hsmm_state_seq,hmm_bbstate_seq,hmm,Fs_bb,Fs] = ...
        pa_get_hsmm_uds_state_seq(wcv_minus_spike,raw_Fs,desynch_times_lf8,sess_data(d).name);
    hmm.run_dur = toc;
    for i = 1:hmm.Nsegs
        curT = length(hmm.UDS_segs(i,1):hmm.UDS_segs(i,2));
        fract_low_post{i}(d) = sum(hmm.posteriors{i}(:) < 0.05)/curT;
    end
    save pa_hsmm_state_seq hsmm* hmm* Fs* fract_low_post*
    clear hsmm* hmm* fract*
    

%% compute LF8 state sequences
        disp('computing LF8 state sequences')
        tic
        [hsmm_bbstate_seq8,hsmm_state_seq8,hmm_bbstate_seq8,hmm8,Fs_bb,Fs] = ...
            pa_get_hsmm_uds_state_seq(lf8,raw_Fs,desynch_times_lf8,strcat(sess_data(d).name,'_lf8'));
        hmm8.run_dur8 = toc;
        for i = 1:hmm8.Nsegs
            curT = length(hmm8.UDS_segs(i,1):hmm8.UDS_segs(i,2));
            fract_low_post8{i}(d) = sum(hmm8.posteriors{i}(:) < 0.05)/curT;
        end
        save pa_hsmm_state_seq8 hsmm* hmm* Fs* fract_low_post*
        clear hsmm* hmm* fract*

end

cd F:\WC_Germany\persistent_9_27_2010
