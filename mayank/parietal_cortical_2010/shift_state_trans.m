function [corrected_state_seq] = shift_state_trans(state_seq,up_shift,down_shift)

up_trans = find(state_seq(1:end-1) == 1 & state_seq(2:end) == 2);
down_trans = find(state_seq(1:end-1) == 2 & state_seq(2:end) == 1);

for i = 1:length(up_trans)
    cur_up = up_trans(i);
    prev_down = down_trans(find(down_trans < cur_up,1,'last'));
    if isempty(prev_down)
        prev_down = 1;
    end
    if cur_up > up_shift && cur_up - prev_down > up_shift
        up_trans(i) = cur_up - up_shift;
    end
end

for i = 1:length(down_trans)
    cur_down = down_trans(i);
    prev_up = up_trans(find(up_trans < cur_down,1,'last'));
    if isempty(prev_up)
        prev_up = 1;
    end
    if cur_down > down_shift && cur_down - prev_up > down_shift
        down_trans(i) = cur_down - down_shift;
    end
end

corrected_state_seq = ones(length(state_seq),1);
for i = 1:length(up_trans)
    next_down = down_trans(find(down_trans > up_trans(i),1,'first'));
    if ~isempty(next_down)
        corrected_state_seq(up_trans(i):next_down) = 2;
    else
        corrected_state_seq(up_trans(i):end) = 2;
    end
end
if down_trans(1) < up_trans(1)
    corrected_state_seq(1:down_trans(1)) = 2;
end
