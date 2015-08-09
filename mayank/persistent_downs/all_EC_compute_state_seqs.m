function [] = all_EC_compute_state_seqs(f_name)

% addpath(genpath(strcat(drive_letter,':\Code\Chronux')))
% addpath(strcat(drive_letter,':\Code\smoothing\software'))
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
% addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
% addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
% addpath(strcat(drive_letter,':\WC_Germany\persistent_9_27_2010\'))
% addpath(strcat(drive_letter,':\WC_Germany\sven_thomas_combined\'))
% addpath(strcat(drive_letter,':\WC_Germany\hsmm_uds_code'))
% addpath(strcat(drive_letter,':\Code\maphmmbox\'))
% addpath('C:\WC_Germany\hsmm_uds_code\')
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))
addpath('~/James_scripts/mayank/hsmm_state_detection/');
addpath('~/James_scripts/mayank/parietal_cortical_2010/');
addpath('~/James_scripts/mayank/persistent_9_27_2010/');
addpath('~/James_scripts/mayank/hsmm_uds_code/');
% addpath(strcat(drive_letter,':\Code\smoothing\software'))
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
% addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
% addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
% addpath(strcat(drive_letter,':\WC_Germany\persistent_9_27_2010\'))
% addpath(strcat(drive_letter,':\WC_Germany\sven_thomas_combined\'))
% addpath(strcat(drive_letter,':\WC_Germany\hsmm_uds_code'))
% addpath(strcat(drive_letter,':\Code\maphmmbox\'))
% addpath('C:\WC_Germany\hsmm_uds_code\')


raw_Fs = 2016;

load ./used_data wcv_minus_spike lf7 %LF7 is cortical LFP, wcv_minus_spike is the de-spiked MP signal
[desynch_times_ctx,desynch_inds,P_ctx,f,t] = locate_desynch_times_individual_v2(lf7); %get indices of desynchronized epochs (where cortical LFP does not have UDS)

%times of UDS epochs
synch_times(:,1) = [0; desynch_times_ctx(:,2)];
synch_times(:,2) = [desynch_times_ctx(:,1); length(lf7)/raw_Fs];
synch_durs = synch_times(:,2)-synch_times(:,1);

%% compute MP state sequences
if max(synch_durs) > 60
    %     [hmm_bbstate_seq,hmm,Fs_bb,Fs] = all_hmm_uds_state_seq(wcv_minus_spike,raw_Fs,desynch_times_ctx,['mp_' f_name]);
    [hsmm_bbstate_seq,hsmm_state_seq,hmm_bbstate_seq,hsmm,hmm,Fs_bb,Fs] = ...
        pa_get_hsmm_uds_state_seq_new(wcv_minus_spike,raw_Fs,desynch_times_ctx,['mp_' f_name]);
else
    hmm_bbstate_seq = [];
    hsmm_bbstate_seq = [];
    hsmm = [];
    hmm = [];
    Fs = [];
end
% save all_combined_mp_uds hmm* Fs*
save pa_hsmm_state_seq_combined_fin_nd hmm* hsmm* Fs*
clear hmm* hsmm* fract*

%% compute LF7 state sequences
disp('computing CTX state sequences')
if max(synch_durs) > 60
%     [hmm_bbstate_seq7,hmm7,Fs_bb,Fs] = all_hmm_uds_state_seq(lf7,raw_Fs,desynch_times_ctx,['lf7_' f_name]);
    [hsmm_bbstate_seq7,hsmm_state_seq7,hmm_bbstate_seq7,hsmm7,hmm7,Fs_bb,Fs] = ...
        pa_get_hsmm_uds_state_seq_new(lf7,raw_Fs,desynch_times_ctx,['lf7_' f_name]);
else
    hsmm_bbstate_seq7 = [];
    hmm_bbstate_seq7 = [];
    hsmm7 = [];
    hmm7 = [];
end
% save all_combined_lf7_uds hmm* Fs*
save pa_hsmm_state_seq7_combined_fin_nd hmm* hsmm* Fs*
clear hmm* hsmm* fract*

%     end
