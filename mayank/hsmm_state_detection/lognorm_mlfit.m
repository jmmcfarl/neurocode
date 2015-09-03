function [mu,var] = lognorm_mlfit(data)

mu = mean(log(data));
var = mean((log(data)-mu).^2);

