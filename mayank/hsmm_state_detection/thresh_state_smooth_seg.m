function [smoothed_seq] = thresh_state_smooth_seg(state_seq,Fs,up_thresh,down_thresh)

%state duration thresholds in ms

% reject_dur = [];
for i = 1:length(state_seq)

    smoothed_seq{i} = state_seq{i}(:);

    state_changes = [1;find(diff(smoothed_seq{i}) ~=0)];
    sojourn_times = diff(state_changes);
    sojourn_times = sojourn_times/Fs*1000; %in ms
    state_changes(end) = [];

    [dummy,sojourn_order] = sort(sojourn_times);
    %sort from shortest to longest states with up states first for proposals
    % sojourn_order_up2 = sojourn_order(state_seq(state_changes(sojourn_order)+1)==3);
    sojourn_order_up = sojourn_order(state_seq{i}(state_changes(sojourn_order)+1)==2);
    sojourn_order_down = sojourn_order(state_seq{i}(state_changes(sojourn_order)+1)==1);
    % sojourn_order = [sojourn_order_up2 sojourn_order_up1 sojourn_order_down];

    % sojourn_order(sojourn_order==1|sojourn_order==length(sojourn_order)) = [];

    while ~isempty(sojourn_order_up)
        cur_state = sojourn_order_up(1);
        if cur_state > 1 && cur_state < length(sojourn_times)
            cur_start = state_changes(cur_state)+1;
            cur_stop = state_changes(cur_state+1);
            prev_state = smoothed_seq{i}(cur_start-1);
            next_state = smoothed_seq{i}(cur_stop+1);
            if prev_state == next_state
                other_state = prev_state;
                if sojourn_times(cur_state) < up_thresh
%                     reject_dur = [reject_dur sojourn_times(cur_state)];
                    smoothed_seq{i}(cur_start:cur_stop) = other_state;
                    state_changes([cur_state cur_state+1]) = [];
                    sojourn_times(cur_state-1) = sojourn_times(cur_state-1)...
                        +sojourn_times(cur_state)+sojourn_times(cur_state+1);
                    sojourn_times([cur_state cur_state+1]) = [];
                    sojourn_order_up(sojourn_order_up==cur_state+1) = [];
                    sojourn_order_up(sojourn_order_up >= cur_state) = ...
                        sojourn_order_up(sojourn_order_up >= cur_state)-2;
                end
            end
        end
        sojourn_order_up(1) = [];
    end

    while ~isempty(sojourn_order_down)
        cur_state = sojourn_order_down(1);
        if cur_state > 1 && cur_state < length(sojourn_times)
            cur_start = state_changes(cur_state)+1;
            cur_stop = state_changes(cur_state+1);
            prev_state = smoothed_seq{i}(cur_start-1);
            next_state = smoothed_seq{i}(cur_stop+1);
            if prev_state == next_state
                other_state = prev_state;
                if sojourn_times(cur_state) < down_thresh
%                     reject_dur = [reject_dur sojourn_times(cur_state)];
                    smoothed_seq{i}(cur_start:cur_stop) = other_state;
                    state_changes([cur_state cur_state+1]) = [];
                    sojourn_times(cur_state-1) = sojourn_times(cur_state-1)...
                        +sojourn_times(cur_state)+sojourn_times(cur_state+1);
                    sojourn_times([cur_state cur_state+1]) = [];
                    sojourn_order_down(sojourn_order_down==cur_state+1) = [];
                    sojourn_order_down(sojourn_order_down >= cur_state) = ...
                        sojourn_order_down(sojourn_order_down >= cur_state)-2;
                end
            end
        end
        sojourn_order_down(1) = [];
    end
end
% fprintf('rejected %d of %d states\n',length(reject_dur),length(orig_dur));
