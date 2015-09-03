function [pmf] = gamma_pmf(x,alpha,beta)

%compute discrete gamma distribution at points in vector x using parameters
%alpha and beta

pmf = exp(-beta*x).*x.^(alpha-1);
pmf = pmf/nansum(pmf);

