function z = logexp(x,K)

alpha = K(1);
beta = K(2);
z = log(1+exp((x-alpha)/beta));