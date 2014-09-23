function [null_peaks,null_valleys] = get_null_peaks(n_sacs,nboot,Robs,backlag,forwardlag,search_range,dt)

NT = length(Robs);
lags = -backlag:forwardlag;
all_tavgs = nan(nboot,length(lags));
for ii = 1:nboot
    randinds = randi(NT,n_sacs,1);
    all_tavgs(ii,:) = get_event_trig_avg_v3(Robs,randinds,backlag,forwardlag);
end

[null_peaks] = get_tavg_peaks((all_tavgs-1),lags*dt,search_range);
[null_valleys] = get_tavg_peaks(-(all_tavgs-1),lags*dt,search_range);

