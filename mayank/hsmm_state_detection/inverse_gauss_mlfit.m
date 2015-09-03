function [mu,lambda] = inverse_gauss_mlfit(data)

mu = nanmean(data);
lambda = (nanmean(1./data - 1/mu))^(-1);

