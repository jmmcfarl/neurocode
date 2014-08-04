function [res,dev] = compute_poisreg_residuals(obs,pred)

zset = obs == 0;
deviance = 2*(obs.*log(obs./pred) - obs + pred);
deviance(zset) = 2*pred(zset);
dev = sum(deviance);
res = sign(obs-pred).*sqrt(deviance);