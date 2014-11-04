clear all
close all
drive_letter = 'C';

addpath(strcat(drive_letter,':\Code\smoothing\software'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
addpath(strcat(drive_letter,':\WC_Germany\persistent_9_27_2010\'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_uds_code'))
addpath(strcat(drive_letter,':\Code\maphmmbox\'))
addpath('C:\WC_Germany\hsmm_uds_code\')

cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
combined_dir = combined_dir(uset);

raw_Fs = 2016;

%%
for d = 1:length(combined_dir)
% for d = [56]
    cd(combined_dir{d})
    pwd
    
    last_slash = find(combined_dir{d} == '\',1,'last');
    f_name = combined_dir{d}(last_slash+1:end);
    
    load ./used_data wcv_minus_spike lf8 lf7 lf6
    if ctx_lfp(d) == 7
        lf8 = lf7;
    end
    if ctx_lfp(d) == 6
        lf8 = lf6;
    end
%     [desynch_times_lf7,desynch_inds,P_lf7,f,t] = locate_desynch_times_individual_v2(lf7);
     [desynch_times_ctx,desynch_inds,P_ctx,f,t] = locate_desynch_times_individual_v2(lf8);
   
    %% compute MP state sequences
    tic
    [hsmm_bbstate_seq,hsmm_state_seq,hmm_bbstate_seq,hsmm,hmm,Fs_bb,Fs] = ...
        pa_get_hsmm_uds_state_seq_new(wcv_minus_spike,raw_Fs,desynch_times_ctx,['mp_' f_name]);
    save pa_hsmm_state_seq_combined_fin_nd hsmm* hmm* Fs*
%     save pa_hsmm_state_seq_combined_fin_newdes hsmm* hmm* Fs*
    clear hsmm* hmm* fract*
    
    %% compute LF8 state sequences
%     disp('computing LF8 state sequences')
%     tic
%     [hsmm_bbstate_seq8,hsmm_state_seq8,hmm_bbstate_seq8,hsmm8,hmm8,Fs_bb,Fs] = ...
%         pa_get_hsmm_uds_state_seq_new(lf8,raw_Fs,desynch_times_lf8,['lf8_' f_name]);
%     save pa_hsmm_state_seq8_combined hsmm* hmm* Fs*
%     clear hsmm* hmm* fract*
    
    %% compute LF7 state sequences
    disp('computing CTX state sequences')
    tic
    [hsmm_bbstate_seq7,hsmm_state_seq7,hmm_bbstate_seq7,hsmm7,hmm7,Fs_bb,Fs] = ...
        pa_get_hsmm_uds_state_seq_new(lf8,raw_Fs,desynch_times_ctx,['lf7_' f_name]);
    save pa_hsmm_state_seq7_combined_fin_nd hsmm* hmm* Fs*
%     save pa_hsmm_state_seq7_combined_fin_newdes hsmm* hmm* Fs*
    clear hsmm* hmm* fract*
    
    %% compute LF3 state sequences
    %     disp('computing LF3 state sequences')
    %     tic
    %     [hsmm_bbstate_seq3,hsmm_state_seq3,hmm_bbstate_seq3,hsmm3,hmm3,Fs_bb,Fs] = ...
    %         pa_get_hsmm_uds_state_seq_new(lf2-lf5,raw_Fs,desynch_times_lf8,['lf3_' f_name]);
    %     save pa_hsmm_state_seq3_combined hsmm* hmm* Fs*
    %     clear hsmm* hmm* fract*
    
end