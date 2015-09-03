function [pmf] = log_norm_pmf(x,mu,var)

pmf = 1./(x*sqrt(var)).*exp(-((log(x)-mu).^2/(2*var)));
pmf = pmf/nansum(pmf);
pmf(isnan(pmf)) = 0;