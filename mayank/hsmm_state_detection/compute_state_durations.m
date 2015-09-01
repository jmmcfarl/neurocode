function [state_durations] = compute_state_durations(state_seq,Fs)

up_trans = find(state_seq(2:end)==2 & state_seq(1:end-1)==1);
down_trans = find(state_seq(2:end)==1 & state_seq(1:end-1)==2);

down_trans(down_trans < up_trans(1)) = [];
up_trans(up_trans > down_trans(end)) = [];

state_durations{1} = (up_trans(2:end)-down_trans(1:end-1))/Fs;
state_durations{2} = (down_trans-up_trans)/Fs;