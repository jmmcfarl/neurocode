function [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data(eye_structure,ts_range)

eye_dt = eye_structure.Header.CRrates(1);
eye_track_starts = [eye_structure.Trials(:).Start]/1e4;
eye_track_ends = [eye_structure.Trials(:).End]/1e4;

first_needed_trial = find(eye_track_starts <= ts_range(1),1,'last');
last_needed_trial = find(eye_track_ends >= ts_range(2),1,'first');
if isempty(first_needed_trial) | isempty(last_needed_trial)
    error('Dont have all eye-tracking data needed')
end

segs_to_get = first_needed_trial:last_needed_trial;
eye_vals = [];
eye_ts = [];
for i = segs_to_get
    cur_rh = eye_structure.Trials(i).Eyevals.rh;
    cur_lh = eye_structure.Trials(i).Eyevals.lh;
    cur_rv = eye_structure.Trials(i).Eyevals.rv;
    cur_lv = eye_structure.Trials(i).Eyevals.lv;
    cur_T1 = (eye_track_starts(i)+eye_dt/2):eye_dt:(eye_track_ends(i)-eye_dt/2);

    eye_data_lens(i,1) = length(cur_lh);
    eye_data_lens(i,2) = length(cur_lv);
    eye_data_lens(i,3) = length(cur_rh);
    eye_data_lens(i,4) = length(cur_rv);
    
    cv = length(cur_lv);
    ch = length(cur_lh);
    if cv > ch
        cur_lv(ch+1:end) = [];
    end
    if cv < ch
        cur_lv = [cur_lv; cur_rv(cv+1:end)];
    end
%     if length(cur_rv) > length(cur_rh)
%         cur_rv(length(cur_rh)+1:end) = [];
%     end
%     if length(cur_lv) < length(cur_rv)
%         cur_lv = [cur_lv; cur_rv(length(cur_lv)+1:end)];
%     elseif length(cur_lv) > length(cur_rv)
%         cur_lv(length(cur_rv)+1:end) = [];
%     end
    

    
    cur_T = (eye_track_starts(i)+eye_dt/2):eye_dt:(eye_track_starts(i)+eye_dt/2+eye_dt*(length(cur_rv)-1));
    if length(cur_T1) ~= length(cur_T)
disp('Warning: Trial parsing issue'); 
    end
    
    eye_ts = [eye_ts; cur_T(:)];
    eye_vals = [eye_vals; cur_lh(:) cur_lv(:) cur_rh(:) cur_rv(:)];
end

dts = [0; diff(eye_ts)];
cur_back_pts = find(dts < 0);
use_samples = ones(size(eye_ts));
for i = 1:length(cur_back_pts)
    next_use = find(eye_ts(cur_back_pts(i):end) > eye_ts(cur_back_pts(i)-1),1,'first');
    use_samples(cur_back_pts(i):cur_back_pts(i)+next_use) = 0;
end
use_samples = logical(use_samples);

eye_ts_interp = eye_ts(1):eye_dt:eye_ts(end);
eye_vals_interp = interp1(eye_ts(use_samples),eye_vals(use_samples,:),eye_ts_interp);
eye_insample = zeros(size(eye_ts_interp));
for i = segs_to_get
    eye_insample(eye_ts_interp > eye_track_starts(i) & eye_ts_interp < eye_track_ends(i)) = 1;
end

