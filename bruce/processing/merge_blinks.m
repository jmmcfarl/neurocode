function [saccades,is_blink] = merge_blinks(saccades,is_blink)

min_ibi = 0.05; %minimum inter-blink interval
sac_start_times = [saccades(:).start_time];
sac_stop_times = [saccades(:).stop_time];

double_blink = find(is_blink(1:end-1) & is_blink(2:end));
inter_blink_interval = sac_start_times(double_blink+1)-sac_stop_times(double_blink);
double_blink = double_blink(inter_blink_interval <= min_ibi);

for ii = 1:length(double_blink)
    new_stop_time = saccades(double_blink(ii)+1).stop_time;
    new_duration = new_stop_time - saccades(double_blink(ii)).start_time;
    saccades(double_blink(ii)).stop_time = new_stop_time;
    saccades(double_blink(ii)).duration = new_duration;
end
saccades(double_blink+1) = [];
is_blink(double_blink+1) = [];