function [Fn, opt] = dnf_arkal_instfreq(Y, samplingrate, varargin)
% [Fn, opt] = dnf_arkal_instfreq(Y, SAMPLINGRATE, varargin)
%
% Kalman Filter for Univariate AR Models
% This function is especially for instantaneous frequency estimation
% for bandlimited signals.
%
% Assumptions:
%  Amplitude demodulation will occur such that at any point in time,
%  oscillations will have a maximum value of 1.
%
% Inputs: 
% Y = input signal (1-D, column oriented), should be bandpass filtered
% SAMPLINGRATE = sampling rate of Y
%  
%   'sigma2v' = defaults to 0.07; observation error var
%   'sigma2w' = defaults to 0.007; state error var
%   'ARorder'   = AR ORDER defaults to 2
%   'unbiased' = remove amplitude bias from tracking, default = 1
%
% Outpus:
%  Estimated AR coefficients and residual E
%  y(n) = A(1)y(n-1) + A(2)y(n-2) + ... + A(P)y(n-P)
%
%
% Copyright (c) 2009, David P. Nguyen <dpnguyen@mit.edu>


if ~exist('samplingrate')
  opt.samplingrate = 2*pi;
else
  opt.samplingrate = samplingrate;
end

%% tvs2e -> unbiased; tvlambda -> smoothness
opt.rescalefactor = [];
opt.ARorder = 2;
opt.Sigma2W = 0.007; % Variance of Kalman state transition process
opt.Sigma2V = 0.07; % Observation Noise Variance
opt.unbiased = 1; % do amplitude demodulation first
opt.boundError = 1; % bound the error so it never blows up
opt.ksmooth = 1; % always ksmooth

opt = parsevarargin(varargin, opt);

% either normalize or amplitude demodulate
Y = Y(:);
if opt.unbiased == 1
  Y = Y./abs(hilbert(Y));
else
  my = abs(max(Y));  
  Y = Y./my;
end


N = length(Y);  % number of samples
M = 1;  % number of channels

maxY = max(abs(Y));
threshE = maxY(:).*2;  

P = opt.ARorder;

%% INITIALIZE THE ALGORITHM
starttime = cputime;

An = zeros(P,1);  % AR coefficient Matrix
B = eye(P,P);
Ptt = eye(P,P);
Pt1step = Ptt;

if opt.ksmooth == 1
  Sm = zeros(P,P,N);
end

Sigma2W = opt.Sigma2W(1);                    % this is variance, not std
Sigma2W = Sigma2W.*eye(P, P); 		     % state transition noise
Sigma2V = opt.Sigma2V(1);                    % this is variance, not std


%warning('off', 'MATLAB:divideByZero');
%warning off;

%% Initalize using the yule walker equations
Xhn = zeros(P, N);  % estimated state (AR coefficients)
atmp = aryule(Y(:), opt.ARorder);
atmp = atmp(:);
Xh = -atmp(2:end);
Xhn(:,P) = Xh;

En = zeros(N,1);        % error
Test = 0;
E = 0;

% push the kalman gain and variances to steady state
n = 1+P;
C = Y(n-1:-1:n-P)';
for nn=1:25
  K = (Pt1step*C')./(C*Pt1step*C' + Sigma2V);
  Pt1step = B*Ptt*B' + Sigma2W;
  Ptt = Pt1step - K*C*Pt1step;
end

%% BIG LOOP
for n = 1+P:N
  
  % Data stacked from t = n-p, to t = 1
  Ys = Y(n-1:-1:n-P);
  
  % Observation matrix
  C = Ys';
 
  % Kalman Gain
  K = (Pt1step*C')./(C*Pt1step*C' + Sigma2V);

  % one-step state error covariance
  Pt1step = B*Ptt*B' + Sigma2W;
   
  % Smoothing matrix
  if opt.ksmooth == 1
    Sm(:,:,n) = Ptt*B'*inv(B*Pt1step*B' + Sigma2W);
  end

  % full state error covariance
  Ptt = Pt1step - K*C*Pt1step;
 
  % innovations
  E = (Y(n) - C*Xh);
  
  % bound the errors to prevent transient instabilities
  % from messing up the whole simulation
  if opt.boundError == 1
    E = min(threshE, E);
    E = max(-threshE, E);
  end
  
  % state estimate correction
  Xh = B*Xh + K*E;
     
  idbad = find(isfinite(Xh) == 0);
  if length(idbad) > 0
    Xh(idbad) = 0;
  end
  
  % save the update
  Xhn(:,n) = Xh;
  En(n) = E;
end

%% SMOOTH THE ESTIMATES
if 1 == opt.ksmooth
  Xhnsm = zeros(size(Xhn));
  Xhnsm(:,end) = Xhn(:,end);
  Xhnsm(:,end-1) = Xhn(:,end);
  for n = N-1:-1:1+P
    sum(abs(Xhnsm(:,n+1) - Xhn(:,n+1)));
    Xhnsm(:,n) = Xhn(:,n) + Sm(:,:,n)*(Xhnsm(:,n+1) - Xhn(:,n)); 
  end
  Xhn = Xhnsm;
  
  % RECALCULATE THE ERROR AFTER SMOOTHING
  for n = 1+P:N
    En(n) = Y(n) - Y(n-1:-1:n-P)'*Xhn(:,n);
  end
end

% RESHAPE X(n) to get A(n)
An = single(Xhn);

% FILL IN THE INITIAL ESTIMATES
for k = 1:P+1
  An(:,k) = An(:,P+2);
  En(k) = En(P+2);
end

opt.En = En;

if 1 == 1

  % CONVERT THE AR COEFFICIENTS TO FREQUENCY
  tmproots = zeros(P,N);
  
  for nn = 1:N
    tmp = roots(fliplr([1, -An(:,nn)']));
    tmproots(:,nn) = [tmp; zeros(P-length(tmp),1)];
  end
  Fn = angle(tmproots);
  Fn = sort(Fn, 1);
  Fn = Fn(2,:);
  Fn = Fn.*(1/pi).*(samplingrate/2);

end

