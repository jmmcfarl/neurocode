function [pmf] = inverse_gauss_pmf(x,mu,lambda)

pmf = (lambda./(x.^3)).^(1/2).*exp(-lambda*(x-mu).^2./(2*mu^2*x));
pmf = pmf/nansum(pmf);
