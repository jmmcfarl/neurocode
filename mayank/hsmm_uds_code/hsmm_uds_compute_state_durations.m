function [state_durations] = hsmm_uds_compute_state_durations(state_seq,Fs)

% compute state durations given a cell array of state sequences
%
% Input: 
%       state_seq: cell array containing a 2-state sequence
%       Fs: sample frequency of the state-seq array
%              
% Output:
%       state_durations: 2-element cell array containing vectors of state
%       durations (in seconds) for the down (1) and up (2) states

state_durations{1} = []; %for down states
state_durations{2} = []; %for up states
for i = 1:length(state_seq)
    state_seq{i} = state_seq{i}(:);
    up_trans = find(state_seq{i}(2:end)==2 & state_seq{i}(1:end-1)==1);
    down_trans = find(state_seq{i}(2:end)==1 & state_seq{i}(1:end-1)==2);

    %start on an up-transition and end on a down transitions
    down_trans(down_trans < up_trans(1)) = [];
    up_trans(up_trans > down_trans(end)) = [];

    %the following insures that there are the same number of up and down
    %states, last down state has undefined duration 
    state_durations{1} = [state_durations{1}; (up_trans(2:end)-down_trans(1:end-1))/Fs; nan];
    state_durations{2} = [state_durations{2}; (down_trans-up_trans)/Fs];
end