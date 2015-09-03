function [mu,lambda] = inverse_gauss_mlfit_pmf(pmf,x)

mu = sum(pmf.*x);
lambda = (sum(pmf.*(1./x - 1/mu)))^(-1);

