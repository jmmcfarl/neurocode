function [pmf] = exponential_pmf(x,lambda)

pmf = exp(-lambda*x);
pmf = pmf/sum(pmf);

