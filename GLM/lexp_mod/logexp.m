function f = logexp(x,beta)

f = 1/beta*log(1+exp(beta*x)) - 1/beta*log(1+exp(beta*min(x)));