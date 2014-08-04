function Q = dnf_arkal_ljungbox(E, ploton) 
% lb_test = dnf_arkal_ljungbox(En, ploton)
% 
% Performs hypothesis test for all lags (usually we consider up to 25)
%  
% The autocorrelation of the residual can be used to test for whiteness.
% If the lags of the autocorr as a whole are bounded by chisquared(K - p)
%  at either 10% or 5% of the distribution, then the model is said to be a
%  good fit, otherwise there will be doubt of the fit.
%
% En is the residual, the output from dnf_arkal.m
% ploton = 0 or 1, plot the results of the test
%
% Box and Jenkins page 314.
% Copyright (c) 2009, David P. Nguyen <dpnguyen@mit.edu>

Klim = 100;
N = numel(E);
N = 200;
ci_bound = 1.96/sqrt(N);

if ~exist('ploton', 'var')
  ploton = 0;
end

% compute the autocorr  
r = xcorr(E(:))./N;

% autocorr is symmetric, so just take the second half
L = floor(numel(r)/2) + 1;
r = r(L:end);
r = r(:)';

% this is ljung-box-pierce statistic
Q = N.*(N+2).*cumsum(1./(N-(1:Klim)).*(r(1:Klim).^2));


if ploton == 1
  subplot(2,1,1); cla; hold off;
  plot(r(1:Klim), 'k', 'linewidth', 2); hold on;
  plot([0 Klim], [ci_bound ci_bound], 'k', 'linewidth', 3);
  set(gca,'box', 'on', 'linewidth', 2, 'fontsize', 14, 'fontweight', 'bold');
  xlim([0 100]);
  ylim([-0.5 0.5]);
  xlabel('Lags in Autcorrelation');
  ylabel('Autocorrelation');
  title('Visual Examination using 95% CI');
  
  subplot(2,1,2); cla; hold off;
  upperlim = chi2inv(0.95, 1:Klim - 2);
  lowerlim = chi2inv(0.90, 1:Klim - 2);
  plot(upperlim, 'k', 'linewidth', 3); hold on;
  %plot(lowerlim, 'k', 'linewidth', 3); 
  plot(Q, 'r--', 'linewidth', 3);
  set(gca,'box', 'on', 'linewidth', 2, 'fontsize', 14, 'fontweight', 'bold');
  axis tight;
  xlim([0 100]);
  xlabel('Lags in Autocorrelation');
  ylabel('\chi ^2 Degrees of Freedom');
  title('Ljung-Box-Pierce Test');
  legend('upper limit', 'data');
end

