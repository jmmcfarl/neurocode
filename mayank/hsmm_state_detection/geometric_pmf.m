function [pmf] = geometric_pmf(x,p)

pmf = (1-p).^(x-1).*p;
pmf = pmf/sum(pmf);

