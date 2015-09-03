function [exp_diff] = invgauss_exp_diff(x,exp_np,dur_range)

mu = x(1);
lambda = x(2);
stat1 = mu;
stat2 = sum(inverse_gaussian_pmf(dur_range,mu,lambda)./dur_range);
exp_diff = abs(stat1 - exp_np(1) + stat2 - exp_np(2));