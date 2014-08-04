function phat = inv_gauss_fit(data)

mu = mean(data);

inv_lambda = mean(1./data - 1./mu);
lambda = 1/inv_lambda;

phat = [mu lambda];