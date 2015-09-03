function [mu,var] = lognorm_mlfit_prob(xvals,pmf)

zvals = find(xvals == 0);
pmf(zvals) = [];
xvals(zvals) = [];

mu = sum(pmf.*log(xvals));
var = sum(pmf.*((log(xvals)-mu).^2));

