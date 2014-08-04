function [fixation_data,in_fixation] = parse_fixations(in_blink,in_sac,eye_ts)

in_fixation = in_blink == 0 & in_sac == 0;

fixation_starts = find(in_fixation(1:end-1) == 0 & in_fixation(2:end) == 1);
fixation_stops = find(in_fixation(1:end-1) == 1 & in_fixation(2:end) == 0);
if in_fixation(1) == 1
    fixation_starts = [1 fixation_starts];
end
if in_fixation(end) == 1
    fixation_stops = [fixation_stops length(in_fixation)];
end

n_fixations = length(fixation_starts);
fixation_start_times = eye_ts(fixation_starts);
fixation_stop_times = eye_ts(fixation_stops);
if n_fixations == 0
    fixation_data = struct('start_times',[],'stop_times',[],'start_inds',[],'stop_inds',[]);
else
    fixation_data = struct('start_times',mat2cell(fixation_start_times(:),ones(n_fixations,1)),...
        'stop_times',mat2cell(fixation_stop_times(:),ones(n_fixations,1)),...
        'start_inds',mat2cell(fixation_starts(:),ones(n_fixations,1)),...
                'stop_inds',mat2cell(fixation_stops(:),ones(n_fixations,1)));
end

% fprintf('%d fixations located\n',n_fixations);