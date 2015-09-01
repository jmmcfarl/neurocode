function [sm_state_seq,state_seq,Fs,fhmm] = parietal_get_fixthresh_uds_state_seq(raw_data,raw_Fs,desynch_times,f_names)

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
N = round(length(obs_lf)/(Fs/hmmFs));
desynch_indices = round(desynch_times*hmmFs);
temp_inds = zeros(N,1);
for i = 1:size(desynch_indices,1)
    temp_inds(desynch_indices(i,1):desynch_indices(i,2)) = 1;
end
UDS_start = find(temp_inds(1:end-1)==1 & temp_inds(2:end)==0);
UDS_stop = find(temp_inds(1:end-1)==0 & temp_inds(2:end)==1);
if temp_inds(1)==0
    UDS_start = [1; UDS_start];
end
if temp_inds(end)==0
    UDS_stop = [UDS_stop; N];
end
UDS_segs = [UDS_start UDS_stop];

UDS_seg_durs = (UDS_segs(:,2)-UDS_segs(:,1))/hmmFs;
too_short = find(UDS_seg_durs < min_seg_dur);
UDS_segs(too_short,:) = [];

%% initialize hsmm
emiss = obs_lf;
[fhmm] = initialize_uds_fixthresh_seg(emiss,Fs,UDS_segs);

[state_seq] = get_fixedthresh_state_seq(fhmm,emiss);

%for estimating a 'smoothed' state sequence, eliminating all states which
%have duration less than some threshold minimum.  Note, this is an
%ambiguous procedure
[sm_state_seq] = fixedthresh_state_seq_smooth(fhmm,emiss,state_seq); 
