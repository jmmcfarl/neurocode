function [up_lags,down_lags,up_times1,up_times2,down_times1,down_times2] = compute_synch_translags(state_seq1,state_seq2,desynch_times,Fs)
addpath('G:\Code\LIBRA_27aug09')

up_trans1 = find(state_seq1(1:end-1) == 1 & state_seq1(2:end)==2);
down_trans1 = find(state_seq1(1:end-1)==2 & state_seq1(2:end)==1);

up_trans2 = find(state_seq2(1:end-1) == 1 & state_seq2(2:end)==2);
down_trans2 = find(state_seq2(1:end-1)==2 & state_seq2(2:end)==1);

up_trans1(up_trans1 > down_trans1(end)) = [];
down_trans1(down_trans1 < up_trans1(1)) = [];

up_trans2(up_trans2 > down_trans2(end)) = [];
down_trans2(down_trans2 < up_trans2(1)) = [];

desynch_ups1 = [];
desynch_ups2 = [];
for i = 1:size(desynch_times,2)
    desynch_ups1 = [desynch_ups1; ...
        find(up_trans1/Fs > desynch_times(1,i) & up_trans1/Fs < desynch_times(2,i))];
    desynch_ups2 = [desynch_ups2; ...
        find(up_trans2/Fs > desynch_times(1,i) & up_trans2/Fs < desynch_times(2,i))];
end

desynch_downs1 = [];
desynch_downs2 = [];
for i = 1:size(desynch_times,2)
    desynch_downs1 = [desynch_downs1; ...
        find(down_trans1/Fs > desynch_times(1,i) & down_trans1/Fs < desynch_times(2,i))];
    desynch_downs2 = [desynch_downs2; ...
        find(down_trans2/Fs > desynch_times(1,i) & down_trans2/Fs < desynch_times(2,i))];
end

up_trans1(desynch_ups1) = [];
down_trans1(desynch_downs1) = [];
up_trans2(desynch_ups2) = [];
down_trans2(desynch_downs2) = [];

up_times1 = up_trans1/Fs;
up_times2 = up_trans2/Fs;
down_times1 = down_trans1/Fs;
down_times2 = down_trans2/Fs;


%% compute relative timing
up_lags = zeros(size(up_trans1));
for i = 1:length(up_trans1)
    [dummy,near_up] = min(abs(up_trans1(i)-up_trans2));
    up_lags(i) = up_trans1(i)-up_trans2(near_up);
end
up_lags = up_lags/Fs;

down_lags = zeros(size(down_trans1));
for i = 1:length(down_trans1)
    [dummy,near_down] = min(abs(down_trans1(i)-down_trans2));
    down_lags(i) = down_trans1(i)-down_trans2(near_down);
end
down_lags = down_lags/Fs;



