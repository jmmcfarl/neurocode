function [pmf] = inverse_gaussian_pmf(x,mu,lambda)

%compute inverse gaussian discrete probability distribution at points in
%vector x using parameters mu and lambda

pmf = (lambda./(x.^3)).^(1/2).*exp(-lambda*(x-mu).^2./(2*mu^2*x));
pmf = pmf/nansum(pmf);
