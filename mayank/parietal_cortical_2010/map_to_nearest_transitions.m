function [new_state_seq] = map_to_nearest_transitions(state_seq1,state_seq2)

up_trans1 = find(state_seq1(1:end-1) == 1 & state_seq1(2:end) == 2);
up_trans2 = find(state_seq2(1:end-1) == 1 & state_seq2(2:end) == 2);
down_trans1 = find(state_seq1(1:end-1) == 2 & state_seq1(2:end) == 1);
down_trans2 = find(state_seq2(1:end-1) == 2 & state_seq2(2:end) == 1);

new_up_trans2 = up_trans1;
new_down_trans2 = down_trans1;
for i = 1:length(up_trans1)
    [dummy,nearest_up2] = min(abs(up_trans1(i)-up_trans2));
    prev_down = down_trans1(find(down_trans1 < up_trans1(i),1,'last'));
    if isempty(prev_down)
        prev_down = 1;
    end
    next_down = down_trans1(find(down_trans1 > up_trans1(i),1,'first'));
    if isempty(next_down)
        next_down = length(state_seq1);
    end
    if up_trans2(nearest_up2) > prev_down && up_trans2(nearest_up2) < next_down
        new_up_trans2(i) = up_trans2(nearest_up2);
    end
end
for i = 1:length(down_trans1)
    [dummy,nearest_down2] = min(abs(down_trans1(i)-down_trans2));
    prev_up = up_trans1(find(up_trans1 < down_trans1(i),1,'last'));
    if isempty(prev_up)
        prev_up = 1;
    end
    next_up = up_trans1(find(up_trans1 > down_trans1(i),1,'first'));
    if isempty(next_up)
        next_up = length(state_seq1);
    end
    if down_trans2(nearest_down2) > prev_up && down_trans2(nearest_down2) < next_up
        new_down_trans2(i) = down_trans2(nearest_down2);
    end
end

%% reconstruct optimized state sequence
new_state_seq = ones(size(state_seq1));
for i = 1:length(new_up_trans2)
    next_down = new_down_trans2(find(new_down_trans2 > new_up_trans2(i),1,'first'));
    if ~isempty(next_down)
        new_state_seq(new_up_trans2(i):next_down) = 2;
    else
        new_state_seq(new_up_trans2(i):end) = 2;
    end
end
if new_down_trans2(1) < new_up_trans2(1)
    new_state_seq(1:new_down_trans2(1)) = 2;
end