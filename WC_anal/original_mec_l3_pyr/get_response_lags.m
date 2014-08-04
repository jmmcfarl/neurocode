function [lags] = get_response_lags(trig_avgs,response_avgs,lag_times)

trig_min = min(trig_avgs,[],2);
[trig_max,trig_max_pt] = max(trig_avgs,[],2);
trig_half = (trig_min+trig_max)/2;

resp_min = min(response_avgs,[],2);
[resp_max,resp_max_pt] = max(response_avgs,[],2);
resp_half = (resp_min+resp_max)/2;

num_cells = size(trig_avgs,1);

trig_half_cross = zeros(1,num_cells);
resp_half_cross = zeros(1,num_cells);
for i = 1:num_cells
   
    temp = find(trig_avgs(i,1:trig_max_pt) < trig_half(i),1,'last');
    if isempty(temp)
        trig_half_cross(i) = nan;
        resp_half_cross(i) = nan;
    else
    trig_half_cross(i) = lag_times(temp);
    
    temp = find(response_avgs(i,1:resp_max_pt) < resp_half(i),1,'last');
    resp_half_cross(i) = lag_times(temp);
    end
end

lags = resp_half_cross-trig_half_cross;