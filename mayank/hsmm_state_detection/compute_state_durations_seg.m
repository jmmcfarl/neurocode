function [state_durations] = compute_state_durations_seg(state_seq,Fs)

state_durations{1} = [];
state_durations{2} = [];
for i = 1:length(state_seq)
    state_seq{i} = state_seq{i}(:);
    up_trans = find(state_seq{i}(2:end)==2 & state_seq{i}(1:end-1)==1);
    down_trans = find(state_seq{i}(2:end)==1 & state_seq{i}(1:end-1)==2);

    if ~isempty(up_trans) && ~isempty(down_trans)
    
    %start on an up-transition and end on a down transitions
    down_trans(down_trans < up_trans(1)) = [];
    up_trans(up_trans > down_trans(end)) = [];

    %the following insures that there are the same number of up and down
    %states, last down state has undefined duration 
    state_durations{1} = [state_durations{1}; (up_trans(2:end)-down_trans(1:end-1))/Fs; nan];
    state_durations{2} = [state_durations{2}; (down_trans-up_trans)/Fs];
    end
end