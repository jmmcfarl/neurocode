function is_blink = detect_blinks_old(all_eye_ts,all_eye_vals,saccades,et_params)

blink_smooth = 0.005;
sac_dur_thresh = 0.1; 
max_blink_amp = Inf; %2 (set to Inf for old version)

all_eye_vel = [zeros(1,4); diff(all_eye_vals)];
if et_params.use_coils(1) == 0
    all_eye_vel(:,[1 2]) = nan;
elseif et_params.use_coils(2) == 0
    all_eye_vel(:,[3 4]) = nan;
end
all_eye_speed = sqrt(nanmax(all_eye_vel.^2,[],2));

all_eye_speed = jmm_smooth_1d_cor(all_eye_speed,round(blink_smooth*et_params.eye_fs));

% saccade_start_times = [saccades(:).start_time];
saccade_start_times = [saccades(:).peak_time];
sac_start_inds = round(interp1(all_eye_ts,1:length(all_eye_ts),saccade_start_times));

thresh_speed = 3*nanmedian(all_eye_speed);
thresh_cups = find(all_eye_speed(1:end-1) < thresh_speed & all_eye_speed(2:end) >= thresh_speed) + 1;
thresh_cdowns = find(all_eye_speed(1:end-1) >= thresh_speed & all_eye_speed(2:end) < thresh_speed) + 1;
if thresh_cdowns(1) < thresh_cups(1)
    thresh_cdowns(1) = [];
end
if thresh_cups(end) > thresh_cdowns(end)
    thresh_cups(end) = [];
end

n_sacs = length(saccades);
sac_dur = nan(n_sacs,1);
for ii = 1:n_sacs
    if all_eye_speed(sac_start_inds(ii)) > thresh_speed
        prev_thresh_cross = find(thresh_cups <= sac_start_inds(ii),1,'last');
        next_thresh_cross = find(thresh_cdowns >= sac_start_inds(ii),1);
        if ~isempty(prev_thresh_cross) && ~isempty(next_thresh_cross)
            sac_dur(ii) = (thresh_cdowns(prev_thresh_cross)-thresh_cups(next_thresh_cross));
        end
    end
end
sac_dur = sac_dur/et_params.eye_fs;
sac_amp = [saccades(:).amplitude];

is_blink = sac_dur > sac_dur_thresh & sac_amp' < max_blink_amp;