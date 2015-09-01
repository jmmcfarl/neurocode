function [pmf] = gamma_pmf(x,alpha,beta)

pmf = exp(-beta*x).*x.^(alpha-1);
pmf = pmf/nansum(pmf);

