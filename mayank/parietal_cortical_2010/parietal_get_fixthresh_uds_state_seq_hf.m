function [sm_state_seq,state_seq,Fs,fhmm] = parietal_get_fixthresh_uds_state_seq_hf(raw_data,raw_Fs,desynch_times,f_names)

drive_letter = 'F';
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))

dsf = 8;
Fs = raw_Fs/dsf;
hmmFs = raw_Fs/40;
niqf = raw_Fs/2;
hf_lcf = 20; %low cut-off for the low-freq filter
hf_hcf = 80; %high cut-off for the low-freq filter
hf_smooth = 0.15;
num_states = 2; %number of hidden states (must be 2 in current implementation)
min_seg_dur = 60;

%% initializations
[obs_hf,t_axis,Fs] = get_hf_features(raw_data,raw_Fs,Fs,[hf_lcf hf_hcf],hf_smooth);

%% compute UDS segments
N = round(length(obs_hf)/(Fs/hmmFs));
desynch_indices = round(desynch_times*hmmFs);
% if size(desynch_indices,1) ~= 0
%     if size(desynch_indices,1)==1
%         if desynch_indices(1,2) < N
%             UDS_segs = [desynch_indices(1,2) N];
%         else
%             UDS_segs = [];
%         end
%     elseif size(desynch_indices,1) > 1
%         UDS_segs = [];
%         for i = 1:size(desynch_indices,1)-1
%             UDS_segs = [UDS_segs; desynch_indices(i,2) desynch_indices(i+1,1)];
%         end
%     end
%     if desynch_indices(1,1) ~= 1
%         UDS_segs = [1 desynch_indices(1,1); UDS_segs];
%     end
% else
%     UDS_segs = [1 N];
% end
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
emiss = obs_hf;
[fhmm] = initialize_uds_fixthresh_seg(emiss,Fs,UDS_segs);

[state_seq] = get_fixedthresh_state_seq(fhmm,emiss);
[sm_state_seq] = fixedthresh_state_seq_smooth(fhmm,emiss,state_seq);
