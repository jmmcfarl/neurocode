function [blink_data,in_blink,tot_disp_f] = get_blinks_v2(rEyeXY,lEyeXY,eye_ts)

%%
blink_thresh = 2;

thresh_eyespeed = 3;
min_interblink_interval = 0.5;
min_blink_dur = 0.05;
eye_dt = median(diff(eye_ts));
eye_fs = 1/eye_dt;
lcf = 0.2;
hcf = 40;
[b,a] = butter(2,[lcf hcf]/(eye_fs/2));

leye_vel = [0 0; diff(lEyeXY)]/eye_dt;
reye_vel = [0 0; diff(rEyeXY)]/eye_dt;
leye_speed = sqrt(leye_vel(:,1).^2+leye_vel(:,2).^2);
reye_speed = sqrt(reye_vel(:,1).^2+reye_vel(:,2).^2);
eye_speed = max(leye_speed,reye_speed);

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
    
    proposed_blink_starts = 1+find(tot_disp_f(1:end-1) < blink_thresh & tot_disp_f(2:end) >= blink_thresh);
    
    blink_start_inds = proposed_blink_starts;
    blink_stop_inds = nan(size(proposed_blink_starts));
    for i = 1:length(proposed_blink_starts)
        temp = find(eye_speed(1:proposed_blink_starts(i)) < thresh_eyespeed,1,'last');
        if ~isempty(temp)
            blink_start_inds(i) = temp;
        end
        temp = find(eye_speed(proposed_blink_starts(i)+1:end) < thresh_eyespeed,1,'first');
        if i < length(proposed_blink_starts)
            if temp > proposed_blink_starts(i+1)-proposed_blink_starts(i)
                temp = nan;
            end
        end
        if ~isempty(temp)
            blink_stop_inds(i) = proposed_blink_starts(i)+temp-1;
        end
    end
    bad_sacs = find(isnan(blink_stop_inds));
    blink_start_inds(bad_sacs) = [];
    blink_stop_inds(bad_sacs) = [];
            
    isi = diff(blink_start_inds);
    too_short_int = find(isi < min_interblink_interval);
    blink_start_inds(too_short_int + 1) = [];
    blink_stop_inds(too_short_int) = [];
  
    blink_durs = (blink_stop_inds-blink_start_inds)*eye_dt;
    too_short = find(blink_durs < min_blink_dur);
    blink_start_inds(too_short) = [];
    blink_stop_inds(too_short) = [];
    
    
    n_blinks = length(blink_start_inds);
    for i = 1:n_blinks
        in_blink(blink_start_inds(i):blink_stop_inds(i)) = 1;
    end
    
    blink_start_times = eye_ts(blink_start_inds);
    blink_stop_times = eye_ts(blink_stop_inds);
    
    if n_blinks > 0
        blink_data = struct('start_times',mat2cell(blink_start_times(:),ones(n_blinks,1)),...
            'stop_times',mat2cell(blink_stop_times(:),ones(n_blinks,1)));
    end
else
    n_blinks = 0;
end

% fprintf('%d blinks located\n',n_blinks);
