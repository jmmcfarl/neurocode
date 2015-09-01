function [mu,lambda] = inverse_gauss_mlfit(data)

%compute maximum likelihood parameters mu and lambda of an inverse gaussian
%distribution fit to the sequence of durations in vector data.

mu = nanmean(data);
lambda = (nanmean(1./data - 1/mu))^(-1);

