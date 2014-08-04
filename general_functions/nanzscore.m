function sig = nanzscore(sig)

sig = bsxfun(@minus,sig,nanmean(sig));
sig = bsxfun(@rdivide,sig,nanstd(sig));