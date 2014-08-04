function [blink_data,in_blink,tot_disp_f] = get_blinks(rEyeXY,lEyeXY,eye_ts)

%%
blink_thresh = 2;
% blink_sep = 0.5;
blink_bound_thresh = 0.6;

thresh_eyespeed = 5;

eye_dt = median(diff(eye_ts));
eye_fs = 1/eye_dt;
lcf = 0.2;
hcf = 40;
[b,a] = butter(2,[lcf hcf]/(eye_fs/2));

%%
h_disp = rEyeXY(:,1) - lEyeXY(:,1);
v_disp = rEyeXY(:,2) - lEyeXY(:,2);
h_disp_f = filtfilt(b,a,h_disp);
v_disp_f = filtfilt(b,a,v_disp);
tot_disp_f = sqrt(h_disp_f.^2+v_disp_f.^2);
% tot_disp = sqrt(h_disp.^2+v_disp.^2);

blink_data = struct('start_times',[],'stop_times',[]);
in_blink = zeros(size(eye_ts));
if max(tot_disp_f) > blink_thresh

    proposed_blink_starts = 1+find(tot_disp_f(1:end-1) < blink_bound_thresh & tot_disp_f(2:end) >= blink_bound_thresh);
    proposed_blink_stops = 1+find(tot_disp_f(1:end-1) >= blink_bound_thresh & tot_disp_f(2:end) < blink_bound_thresh);
    if tot_disp_f(1) >= blink_bound_thresh
        proposed_blink_starts = [1; proposed_blink_starts];
    end
    if tot_disp_f(end) >= blink_bound_thresh
        proposed_blink_stops = [proposed_blink_stops; length(tot_disp_f)];
    end
    fail_blinks = [];
    for i = 1:length(proposed_blink_starts)
        if max(tot_disp_f(proposed_blink_starts(i):proposed_blink_stops(i))) < blink_thresh
            fail_blinks = [fail_blinks i];
        end
    end
    proposed_blink_starts(fail_blinks) = [];
    proposed_blink_stops(fail_blinks) = [];
    
    n_blinks = length(proposed_blink_starts);
    for i = 1:n_blinks
        in_blink(proposed_blink_starts(i):proposed_blink_stops(i)) = 1;
    end
    
    blink_start_times = eye_ts(proposed_blink_starts);
    blink_stop_times = eye_ts(proposed_blink_stops);
    
    if n_blinks > 0
        blink_data = struct('start_times',mat2cell(blink_start_times(:),ones(n_blinks,1)),...
            'stop_times',mat2cell(blink_stop_times(:),ones(n_blinks,1)));
    end
else
    n_blinks = 0;
end

% fprintf('%d blinks located\n',n_blinks);
