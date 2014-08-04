function ic = dnf_arkal_aic(En, p)
% ic = dnf_arkal_aic(En, p)
%
% En = Residual vector from output of Kalman Filter
% p = model order (p = 2, default)
%
% Output: Akaike, Final Prediction Error, Rissanen, Schwarz, 
%         Hannan-Quinn
%
% Ref: Ivanov (2005), Nonlinear Dynamics.
%
% Copyright (c) 2009, David P. Nguyen <dpnguyen@mit.edu>
  
  ic.aic = [];
  ic.fpe = [];
  ic.ris = [];
  ic.sic = [];
  ic.hqc = [];
  ic.mse = [];
  
  En = En(En ~= nan);
  En = En(isfinite(En(:)) == 1);
  En = En(:);
      
  N = numel(En);  % number of observations
  V = var(En);
    
  ic.p = p;
  ic.aic = log(V) + 2*L*L*p/N;
  ic.fpe = (((N+L.*p+1)./(N-L.*p-1)).^L).*V;
  ic.ris = N*log(V) + L*L*p*log(N);
  ic.hqc = log(V) + 2*log(log(N))*p*L*L/N;
  ic.sic = log(V) + log(N)*p*L*L/N;
  ic.mse = mean(En.^2);
  ic.labels = {'AKAIKE', 'FPE', 'RISSANEN', 'HANNAN-QUINN','SCHWARTZ', ...
               'MEAN2ERROR'};
