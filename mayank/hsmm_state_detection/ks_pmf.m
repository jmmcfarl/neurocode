function [ks] = ks_pmf(pmf1,pmf2)

cmf1 = cumsum(squeeze(pmf1));
cmf2 = cumsum(squeeze(pmf2));

ks = max(abs(cmf1-cmf2));