function [pmf] = geometric_pmf(x,p)

%compute geometric probability distribution over points in x using
%parameter p

pmf = (1-p).^(x-1).*p;
pmf = pmf/sum(pmf);

