function [alpha,beta] = gamma_mlfit_prob(pmf,xvals)

zvals = find(xvals == 0);
pmf(zvals) = [];
xvals(zvals) = [];

tol = 1e-4;
delta = 1;

%first estimate sufficient stats
st1 = sum(pmf.*xvals);
st2 = sum(pmf.*log(xvals));

s=log(st1)-st2;
alpha_cur = (3-s+sqrt((s-3)^2+24*s))/(12*s);
alpha_next = alpha_cur;
while delta > tol
    alpha_cur = alpha_next;
    alpha_next = alpha_cur - (log(alpha_cur)-psi(alpha_cur)-s)/(1/alpha_cur-psi(1,alpha_cur));
    delta = abs((alpha_cur-alpha_next)/alpha_cur);
end

alpha = alpha_next;
beta = (1/(alpha)*st1)^(-1);