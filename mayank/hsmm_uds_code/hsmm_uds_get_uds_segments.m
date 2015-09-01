function UDS_segs = hsmm_uds_get_uds_segments(desynch_times,Fs,T,min_seg_dur)

% get the sample indices corresponding to the start and stops of each UDS
% segemnt
%
% Input: 
%       desynch_times:  Nx2 matrix containing the stop and start times (in
%       s) of each desynchronized epoch (epoch without UDS)
%       Fs:  sample frequency of the observation sequence to be used for
%       UDS inference
%       T: length of the observation sequence
%       min_seg_dur:  minimum duration of a UDS segment used for analysis
%       (in seconds). Default = 60s
%              
% Output:
%        UDS_segs: Nx2 matrix containing the start and stop indices (wrt to the observation sequence) 
%         of each segment containing UDS

if nargin < 4
    min_seg_dur = 60; %set default min_seg_dur at 60s
end

desynch_indices = round(desynch_times*Fs);
temp_inds = zeros(T,1);
for i = 1:size(desynch_indices,1)
    temp_inds(desynch_indices(i,1):desynch_indices(i,2)) = 1;
end
UDS_start = find(temp_inds(1:end-1)==1 & temp_inds(2:end)==0);
UDS_stop = find(temp_inds(1:end-1)==0 & temp_inds(2:end)==1);
if temp_inds(1)==0
    UDS_start = [1; UDS_start];
end
if temp_inds(end)==0
    UDS_stop = [UDS_stop; T];
end
UDS_segs = [UDS_start UDS_stop]; 

%get rid of UDS segments which are too short for analysis
UDS_seg_durs = (UDS_segs(:,2)-UDS_segs(:,1))/Fs; %segment durations in seconds
too_short = UDS_seg_durs < min_seg_dur;
UDS_segs(too_short,:) = [];

