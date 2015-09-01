function [state_dutycycles] = compute_state_dutycycles_seg(state_seq)

state_dutycycles = [];
for i = 1:length(state_seq)
    state_seq{i} = state_seq{i}(:);
    up_trans = find(state_seq{i}(2:end)==2 & state_seq{i}(1:end-1)==1);
    down_trans = find(state_seq{i}(2:end)==1 & state_seq{i}(1:end-1)==2);

    down_trans(down_trans < up_trans(1)) = [];
    up_trans(up_trans > down_trans(end)) = [];

    down_durs = up_trans(2:end)-down_trans(1:end-1);
    up_durs = down_trans-up_trans;
    
    cur_dutycycles = nan(size(up_durs));
    cur_dutycycles(2:end) = up_durs(2:end)./(up_durs(2:end)+down_durs);
    
    state_dutycycles = [state_dutycycles; cur_dutycycles];
end