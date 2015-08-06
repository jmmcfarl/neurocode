function is_blink = detect_blinks(all_eye_ts,all_eye_vals,saccades,et_params)
% is_blink = detect_blinks(all_eye_ts,all_eye_vals,saccades,et_params)
% INPUTS: 
%     all_eye_ts: time samples of eye position data
%     all_eye_vals: raw eye position data
%     saccades: struct array of saccade info
%     et_params; param struct
% OUTPUTS: 
%     is_blink: bool array indicating whether each sacccade gets counted as a blink

blink_smooth = 0.005; %smoothing sigma (s) for estimating eye speed
sac_dur_thresh = 0.1; %maximum saccade duration
max_blink_amp = 2; % maximum amplitude to be classified as a blink, rather than a slow sac (set to Inf for old version)

all_eye_vel = [zeros(1,4); diff(all_eye_vals)];
if et_params.use_coils(1) == 0
    all_eye_vel(:,[1 2]) = nan;
elseif et_params.use_coils(2) == 0
    all_eye_vel(:,[3 4]) = nan;
end
all_eye_speed = sqrt(nanmax(all_eye_vel.^2,[],2)); %get maximum instantaneous speed over any usable eye velocity dimension
all_eye_speed = jmm_smooth_1d_cor(all_eye_speed,round(blink_smooth*et_params.eye_fs)); %smooth speed estimate

% saccade_start_times = [saccades(:).start_time];
saccade_start_times = [saccades(:).peak_time];
sac_start_inds = round(interp1(all_eye_ts,1:length(all_eye_ts),saccade_start_times));

thresh_speed = 3*nanmedian(all_eye_speed); %set threshold eye speed
thresh_cups = find(all_eye_speed(1:end-1) < thresh_speed & all_eye_speed(2:end) >= thresh_speed) + 1; %upward threshold crossings
thresh_cdowns = find(all_eye_speed(1:end-1) >= thresh_speed & all_eye_speed(2:end) < thresh_speed) + 1; %downward threshold crossings
if thresh_cdowns(1) < thresh_cups(1)
    thresh_cdowns(1) = [];
end
if thresh_cups(end) > thresh_cdowns(end)
    thresh_cups(end) = [];
end

n_sacs = length(saccades);
sac_dur = nan(n_sacs,1);
for ii = 1:n_sacs
    if all_eye_speed(sac_start_inds(ii)) > thresh_speed %if eye speed is above threhsold at start of sac
        prev_thresh_cross = find(thresh_cups <= sac_start_inds(ii),1,'last'); %find previous speed thresh cross
        next_thresh_cross = find(thresh_cdowns >= sac_start_inds(ii),1); %find next speed down-cross
        if ~isempty(prev_thresh_cross) && ~isempty(next_thresh_cross)
            sac_dur(ii) = (thresh_cdowns(prev_thresh_cross)-thresh_cups(next_thresh_cross)); %get saccade duration
        end
    end
end
sac_dur = sac_dur/et_params.eye_fs;
sac_amp = [saccades(:).amplitude];

%blinks are defined as having duration > sac_dur_thresh and amplitude <
%max_blink_amp
is_blink = sac_dur > sac_dur_thresh & sac_amp' < max_blink_amp;