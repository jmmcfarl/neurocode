function [exp_diff] = gamma_exp_diff(x,exp_np)

alpha = x(1);
beta = x(2);

exp_diff = abs(alpha/beta - exp_np(1) + psi(alpha) - log(beta) - exp_np(2));