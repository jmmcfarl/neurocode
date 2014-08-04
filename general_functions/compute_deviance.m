function deviance = compute_deviance(obs_Y,pred_Y)

obs_Y = obs_Y(:); pred_Y = pred_Y(:);
obs_LL = sum(obs_Y.*log(pred_Y)-pred_Y);
sat_LL = obs_Y.*log(obs_Y)-obs_Y;
sat_LL(obs_Y==0) = 0;
sat_LL = sum(sat_LL);

deviance = -2*(obs_LL - sat_LL);