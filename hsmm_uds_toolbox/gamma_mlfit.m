function [alpha,beta] = gamma_mlfit(data)

%compute maximum likelihood parameters alpha and beta of a gamma
%distribution fit to the sequence of durations in vector data.

N = length(data);
tol = 1e-4;
delta = 1;

%first estimate alpha numerically
s=log(1/N*nansum(data))-1/N*nansum(log(data));
alpha_cur = (3-s+sqrt((s-3)^2+24*s))/(12*s);
alpha_next = alpha_cur;
while delta > tol
    alpha_cur = alpha_next;
    alpha_next = alpha_cur - (log(alpha_cur)-psi(alpha_cur)-s)/(1/alpha_cur-psi(1,alpha_cur));
    delta = abs((alpha_cur-alpha_next)/alpha_cur);
end

alpha = alpha_next;

%now compute beta given ML alpha
beta = (1/(alpha*N)*nansum(data))^(-1);