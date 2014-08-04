function [fit_params,fit_z] = fitGaussianCurve(xi,zi)


sz = length(zi);
boundx = [min(xi),max(xi)];
rgx = diff(boundx);
minsigma = rgx/sz/5;
maxsigma = rgx*.3;

thresh_cut = prctile(zi,15);
zi_thresh = zi - thresh_cut;
zi_thresh(zi_thresh < 0) = 0;

mean_est = sum(xi.*zi_thresh)/sum(zi_thresh);

[max_val,max_loc] = max(zi);
sp = find(zi > max_val/2,1);
ep = find(zi(max_loc:end) < max_val/2,1);
if ~isempty(sp) & ~isempty(ep)
    fwhm = ep+max_loc-1-sp;
    std_est = fwhm/2;
else
    xi_diff = (xi - mean_est).^2;
    std_est = sqrt(sum(xi_diff.*zi_thresh)/sum(zi_thresh));
end

x0 = [mean_est std_est max(zi_thresh) thresh_cut];

lb = [min(xi) minsigma 0 min(zi)];
ub = [max(xi) maxsigma Inf max(zi)];
opts.Display = 'off';
fit_params = lsqcurvefit(@gaussfun,x0,xi,zi,lb,ub,opts);
fit_z = gaussfun(fit_params,xi);

end

function G = gaussfun(x,xdata)

G = exp(-(xdata - x(1)).^2/(2*x(2)^2))*x(3) + x(4);

end
