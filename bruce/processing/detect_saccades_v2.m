function [saccades,et_params] = detect_saccades_v2(all_eye_vals,orig_eye_vals,all_eye_speed,all_eye_ts,et_params)
% [saccades,et_params] = detect_saccades_v2(all_eye_vals,orig_eye_vals,all_eye_speed,all_eye_ts,et_params)
% INPUTS: 
%     all_eye_vals: corrected eye position data (orig sampling)
%     orig_eye_vals: raw (uncorrected) eye position data
%     all_eye_speed: instantaneous eye speed
%     all_eye_ts: time axis for eye data
%     et_params: struct of params
% OUTPUTS:
%     saccades: struct array with saccade info
%     et_params: new param struct


%% SACCADE DETECTION
fprintf('Detecting saccades\n');

%parameters
sac_thresh = 10; %threshold eye speed
peri_thresh = 3; %threshold eye speed for defining saccade boundary inds
min_isi = 0.05; max_isi = Inf; %min/max inter-saccade intervals

%find local maxima of eye speed signal exceeding sac_thresh
peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);

%find times when speed signal crossed above and below the peri-saccade
%threshold
thresh_cross_up = 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh);
thresh_cross_down = 1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh);
sac_start_inds = nan(size(saccade_inds));
sac_stop_inds = nan(size(saccade_inds));
for ii = 1:length(saccade_inds)
    next_tc = find(thresh_cross_down > saccade_inds(ii),1,'first');
    if ~isempty(next_tc)
        sac_stop_inds(ii) = thresh_cross_down(next_tc);
    end
    prev_tc = find(thresh_cross_up < saccade_inds(ii),1,'last');
    if ~isempty(prev_tc)
        sac_start_inds(ii) = thresh_cross_up(prev_tc);
    end
end

%get rid of double-peaks
isis = [Inf; diff(sac_start_inds)]/et_params.eye_fs;
bad_isis = (isis < min_isi | isis > max_isi);
bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];

next_isi = [isis(2:end); Inf];
saccade_times = all_eye_ts(saccade_inds); %saccade peak times
sac_start_times = all_eye_ts(sac_start_inds); %saccade start times
sac_stop_times = all_eye_ts(sac_stop_inds); %saccade end times
sac_durs = sac_stop_times - sac_start_times; %saccade durations 
sac_peakvel = all_eye_speed(saccade_inds); %peak eye vels

sac_dbuff = round(0.0*et_params.eye_fs); %for computing saccade amplitudes, option to include a buffer to account for transients
pre_inds = sac_start_inds - sac_dbuff;
pre_inds(pre_inds < 1) = 1;
sac_pre_pos = all_eye_vals(pre_inds,:); %pre-sac corrected position
sac_pre_pos_raw = orig_eye_vals(pre_inds,:); %pre-sac raw position
post_inds = sac_stop_inds + sac_dbuff;
post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
sac_post_pos = all_eye_vals(post_inds,:); %post sac corrected position
sac_post_pos_raw = orig_eye_vals(post_inds,:); %post-sac raw position

if et_params.use_coils(1) && ~et_params.use_coils(2) %if using left coil only, just use left-position data
    sac_pre_pos = sac_pre_pos(:,1:2);
    sac_post_pos = sac_post_pos(:,1:2);
    sac_pre_pos_raw = sac_pre_pos_raw(:,1:2);
    sac_post_pos_raw = sac_post_pos_raw(:,1:2);
elseif ~et_params.use_coils(1) && et_params.use_coils(2) %if using right coil only, ditto
    sac_pre_pos = sac_pre_pos(:,3:4);
    sac_post_pos = sac_post_pos(:,3:4);
    sac_pre_pos_raw = sac_pre_pos_raw(:,3:4);
    sac_post_pos_raw = sac_post_pos_raw(:,3:4);
elseif et_params.use_coils(1) && et_params.use_coils(2) %if using both coils, take average positions
    sac_pre_pos = 0.5*sac_pre_pos(:,1:2) + 0.5*sac_pre_pos(:,3:4);
    sac_post_pos = 0.5*sac_post_pos(:,1:2) + 0.5*sac_post_pos(:,3:4);
    sac_pre_pos_raw = 0.5*sac_pre_pos_raw(:,1:2) + 0.5*sac_pre_pos_raw(:,3:4);
    sac_post_pos_raw = 0.5*sac_post_pos_raw(:,1:2) + 0.5*sac_post_pos_raw(:,3:4);
end

sac_delta_pos = sac_post_pos - sac_pre_pos; %change in position (corrected) 
sac_delta_pos_raw = sac_post_pos_raw - sac_pre_pos_raw; %change in position (raw) 
sac_amps = sqrt(sum(sac_delta_pos.^2,2)); %sacadde amplitude 
%minus sign on x-coordinate gives these directions in SCREEN coordinates!
sac_dirs = atan2(sac_delta_pos_raw(:,2),-sac_delta_pos_raw(:,1));

temp = ones(length(saccade_times),1);
saccades = struct('peak_time',mat2cell(saccade_times,temp),'start_time',mat2cell(sac_start_times,temp),...
    'stop_time',mat2cell(sac_stop_times,temp),'isi',mat2cell(isis,temp),'next_isi',mat2cell(next_isi,temp),...
    'duration',mat2cell(sac_durs,temp),'amplitude',mat2cell(sac_amps,temp),'direction',mat2cell(sac_dirs,temp),...
    'pre_pos',mat2cell(sac_pre_pos,temp),'post_pos',mat2cell(sac_post_pos,temp),'peakvel',mat2cell(sac_peakvel,temp));

et_params.sac_thresh = sac_thresh;
et_params.peri_thresh = peri_thresh;
et_params.min_isi = min_isi;
et_params.max_isi = max_isi;
et_params.sac_dbuff = sac_dbuff;
