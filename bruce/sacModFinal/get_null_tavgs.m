function [null_tavgs] = get_null_tavgs(n_sacs,nboot,Robs,backlag,forwardlag)

NT = length(Robs);
lags = -backlag:forwardlag;
null_tavgs = nan(nboot,length(lags));
for ii = 1:nboot
    randinds = randi(NT,n_sacs,1);
    null_tavgs(ii,:) = get_event_trig_avg_v3(Robs,randinds,backlag,forwardlag);
end


