function [fm_state_seq,fhmm,Fs] = parietal_get_fixthresh_uds_state_seq_new(raw_data,raw_Fs,desynch_times,thresh_type)

drive_letter = 'G';
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))

dsf = 8;
Fs = raw_Fs/dsf;
hmmFs = raw_Fs/40;
niqf = raw_Fs/2;
lf_lcf = 0.05; %low cut-off for the LF amp filter (Default 0.05)
lf_hcf = 2; %high cut-off for the LF amp filter
num_states = 2; %number of hidden states (must be 2 in current implementation)
min_seg_dur = 60; %minimum duration of a UDS segment

%% initializations
[obs_lf,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs,[lf_lcf lf_hcf]);

%% compute UDS segments
UDS_segs = hsmm_uds_get_uds_segments(desynch_times,Fs,length(obs_lf),min_seg_dur);

%% initialize hsmm
emiss = obs_lf;
if thresh_type == 1
    [fm_state_seq,fhmm] = hsmm_uds_get_fixthresh_uds_state_seq(emiss,Fs,UDS_segs);
elseif thresh_type == 2
    [fm_state_seq,fhmm] = hsmm_uds_get_fixthresh_uds_state_seqz(emiss,Fs,UDS_segs);
elseif thresh_type == 3
    [fm_state_seq,fhmm] = hsmm_uds_get_fixthresh_uds_state_seq_mix(emiss,Fs,UDS_segs);
else
    error('Invalid threshold type')
end

fhmm.UDS_segs = UDS_segs;

