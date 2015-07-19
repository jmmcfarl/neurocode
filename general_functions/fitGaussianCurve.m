function [fit_params,fit_z] = fitGaussianCurve(xi,zi)
%[fit_params,fit_z] = fitGaussianCurve(xi,zi)
%finds least squares gaussian curve fit of the data zi on domain xi,
%fit_params contains the parameters: [mu SD A c] where A is the amplitude
%of the gaussian and c is a constant offset. fit_z is the gaussian function
%evaluated on xi.

sz = length(zi);
boundx = [min(xi),max(xi)];
rgx = diff(boundx);
minsigma = rgx/sz/5; %smallest allowed sigma
maxsigma = rgx*.3; %largest allowed sigma

%% INITIALIZE USING MEAN ESTIMATE
%can get a more robust initial estimate of mu by thresholding zi to zero at
%~15th percentile
thresh_cut = prctile(zi,15); 
zi_thresh = zi - thresh_cut;
zi_thresh(zi_thresh < 0) = 0;
mean_est = sum(xi.*zi_thresh)/sum(zi_thresh); %initial estimate of the mean, treating zi as a distribution

%now get an initial estimate of the SD
[max_val,max_loc] = max(zi);
sp = find(zi > max_val/2,1); %first crossing above half-max
ep = find(zi(max_loc:end) < max_val/2,1); %first crossing under half-max
if ~isempty(sp) & ~isempty(ep) %if both exist, make initial SD estimate FWHM/2
    fwhm = ep+max_loc-1-sp; 
    std_est = fwhm/2;
else %otherwise estimate the SD treating zi as a distribution
    xi_diff = (xi - mean_est).^2; 
    std_est = sqrt(sum(xi_diff.*zi_thresh)/sum(zi_thresh));
end

%% First find least squares fit using the initial mean estimate 
x0 = [mean_est std_est max(zi_thresh) thresh_cut];
lb = [min(xi) minsigma 0 min(zi)];
ub = [max(xi) maxsigma Inf max(zi)];
opts.Display = 'off';
[fit_params1,rnorm1] = lsqcurvefit(@gaussfun,x0,xi,zi,lb,ub,opts);

%% NOW find best fit using max(zi) as the initial estimate of mu
x0 = [max_loc std_est max(zi_thresh) thresh_cut];
[fit_params2,rnorm2] = lsqcurvefit(@gaussfun,x0,xi,zi,lb,ub,opts);

%% find which one is the better fit
if rnorm2 < rnorm1
    fit_params = fit_params2;
else
    fit_params = fit_params1;
end
fit_z = gaussfun(fit_params,xi);

end

function G = gaussfun(x,xdata)

G = exp(-(xdata - x(1)).^2/(2*x(2)^2))*x(3) + x(4);

end
