function [up_state_ind,down_state_ind] = get_used_state_trans(state_seq,Fs,desynch_times,edge_buffer)

up_trans = find(state_seq(1:end-1)==1 & state_seq(2:end) == 2);
down_trans = find(state_seq(1:end-1)==2 & state_seq(2:end) == 1);

%get rid of transitions near the start and end
up_trans(up_trans < edge_buffer) = [];
down_trans(down_trans < edge_buffer) = [];
up_trans(up_trans > length(state_seq)-edge_buffer) = [];
down_trans(down_trans > length(state_seq)-edge_buffer) = [];

up_trans(up_trans > down_trans(end)) = [];
down_trans(down_trans < up_trans(1)) = [];

bad_trans = [];
for i = 1:size(desynch_times,1)
    bad_trans = [bad_trans find(up_trans/Fs > desynch_times(i,1) & up_trans/Fs < desynch_times(i,2))];
end
up_trans(bad_trans) = [];
down_trans(bad_trans) = [];

up_start = up_trans;
up_stop = down_trans;

down_start = down_trans(1:end-1);
down_stop = up_trans(2:end);

up_state_ind = [up_start(:) up_stop(:)];
down_state_ind = [down_start(:) down_stop(:)];

