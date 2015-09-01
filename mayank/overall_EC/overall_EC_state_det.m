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
addpath(strcat(drive_letter,':\WC_Germany\overall_EC\'))

cd(strcat(drive_letter,':\WC_Germany\overall_EC'))
load overall_EC_dir
raw_Fs = 2016;

n = length(sess_data);
for d = 36
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

    load used_data lf8 wcv_minus_spike
    load desynch_times_lf8
%     load used_data lf2 lf5 
%     load desynch_times_lf2r
%     lf2 = lf2/sess_data(d).gains(2);
%     lf5 = lf5/sess_data(d).gains(5);
%     lf2_r = lf2-lf5;
%     if max(isnan(lf2_r)) > 0
%         bad_lf2r(d) = 1;
%     else
%         bad_lf2r(d) = 0;
%     end

load used_data lf3 lf8
% load desynch_times_lf3
% 
%     total_time = length(lf2)/raw_Fs;
%     desynch_time = sum(desynch_times_lf2r(:,2)-desynch_times_lf2r(:,1));
%     uds_time = total_time - desynch_time;
    
    total_time = length(lf3)/raw_Fs;
    desynch_time = sum(desynch_times_lf8(:,2)-desynch_times_lf8(:,1));
    uds_time = total_time - desynch_time;
    
%% compute MP state sequences
% if ~exist('./ec_hmm_state_seq.mat','file')
%     disp('computing MP state sequences')
%     [hmm_bbstate_seq,hmm,Fs_bb,Fs] = ...
%         overall_ec_get_hsmm_uds_state_seq(wcv_minus_spike,raw_Fs,desynch_times_lf8,strcat(s_name,'_MP'));
%     save ec_hmm_state_seq hmm* Fs* 
%     clear hmm* fract*
% end
%% compute LF8 state sequence
% if ~exist('./ec_hmm_state_seq8.mat','file')
    disp('computing LF8 state sequences')
    [hmm_bbstate_seq8,hmm8,Fs_bb,Fs] = ...
        overall_ec_get_hsmm_uds_state_seq(lf8,raw_Fs,desynch_times_lf8,strcat(s_name,'_LF8'));
%     save ec_hmm_state_seq8 hmm* Fs* 
    clear hmm* fract*
% end
%% compute LF2r state sequence
% if bad_lf2r(d) == 0 && uds_time > 180 
%     if ~exist('./ec_hmm_state_seq2r.mat','file')
%     disp('computing LF2r state sequences')
%     [hmm_bbstate_seq2r,hmm2r,Fs_bb,Fs] = ...
%         overall_ec_get_hsmm_uds_state_seq(lf2_r,raw_Fs,desynch_times_lf2r,strcat(s_name,'_LF2r'));
%     save ec_hmm_state_seq2r hmm* Fs* 
%     clear hmm* fract*
%     end
% end

%% compute LF3 state sequence
% if uds_time > 250
%     disp('computing LF3 state sequences')
% %     [hmm_bbstate_seq3,hmm3,Fs_bb,Fs] = ...
% %         overall_ec_get_hsmm_uds_state_seq(lf3,raw_Fs,desynch_times_lf8,strcat(s_name,'_LF3'));
% [bb_hmm_state_seq,hmm_state_seq,hmm,Fs] = parietal_get_hmm_uds_state_seq_hf...
%     (lf3,raw_Fs,desynch_times_lf8,'_LF3');    
% % save ec_hmm_state_seq3 hmm* Fs* 
%     clear hmm* fract*
% end
% 
end

cd G:\WC_Germany\overall_EC
