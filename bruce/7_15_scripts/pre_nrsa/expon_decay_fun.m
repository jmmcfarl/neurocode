function y = expon_decay_fun(beta,x)

y = 1-beta(1)*(1-exp(-x/beta(2)));