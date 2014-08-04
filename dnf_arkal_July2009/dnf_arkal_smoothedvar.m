function [An, En, Vn, opt] = dnf_arkal_smoothedvar(Y, samplingrate, varargin)
% [An, En, SAn, opt] = dnf_arkal_smoothedvar(Y, SAMPLINGRATE, varargin)
%
% Kalman Filter for Univariate AR Models
%
%
% The main difference between this file and others is the variance of the
% state estimate, conditioned on all the data, is also computed (SAn)
%
%
% Assumptions:
%  Amplitude normalization will automatically occur inside the function,
%  i.e. the amplitude of the entire signal will be adjusted to have a
%  maximum value of 1.
%
%  If the 'unbiased' option below is left on, amplitude demodulation will
%  also occur.
%
% Inputs: 
% Y = input signal (1-D, column oriented), should be bandpass filtered
% SAMPLINGRATE = sampling rate of Y
%  
%   'sigma2v' = defaults to 0.07; observation error var
%   'sigma2w' = defaults to 0.007; state error var
%   'ARorder'   = AR ORDER defaults to 2
%   'unbiased' = remove amplitude bias from tracking, default = 1
%   'boundError' = bound the residual to prevent instability, default 1
%
% Outputs:
%  (1) Estimated AR coefficients and (2) residual E
%  y(n) = A(1)y(n-1) + A(2)y(n-2) + ... + A(P)y(n-P)
%
%  If the output is unstable try lowering the values of sigma2v and
%  sigma2w in steps of 0.05.  In addition, unbiased=1 will give more
%  unstabilities over time.
%
%
% 
% Copyright (c) 2009, David P. Nguyen <dpnguyen@mit.edu>


if ~exist('samplingrate', 'var')
  opt.samplingrate = 2*pi;
else
  opt.samplingrate = samplingrate;
end

%% tvs2e -> unbiased; tvlambda -> smoothness
opt.rescalefactor = [];
opt.ARorder = 2;
opt.Sigma2W = 0.007; % Variance of Kalman state transition process
opt.Sigma2V = 0.07; % Observation Noise Variance
opt.unbiased = 1; % compute time-varying variance
opt.boundError = 1; % bound the error so it never blows up
opt.ksmooth = 1; % always ksmooth

opt = parsevarargin(varargin, opt);

% always normalize or amplitude demodulate
Y = Y(:);
my = abs(max(Y));  
Y = Y./my;
opt.rescalefactor = my;

if opt.unbiased == 1
  Y = Y./abs(hilbert(Y));
end

N = length(Y);  % number of samples

maxY = max(abs(Y));
threshE = maxY(:);  

P = opt.ARorder;

%% INITIALIZE THE ALGORITHM

B = eye(P,P);
Ptt = repmat(eye(P,P), [1 1 N]);  % state variance
Pt1step = Ptt;     % 1step state variance
Vn = Ptt;          % fixed interval smoothed state variance

if opt.ksmooth == 1
  Sm = zeros(P,P,N);
end

Sigma2W = opt.Sigma2W(1);                    % this is variance, not std
Sigma2W = Sigma2W.*eye(P, P);                % state transition noise
Sigma2V = opt.Sigma2V(1);                    % this is variance, not std

warning('off', 'MATLAB:divideByZero');

%% Initalize using the yule walker equations
Xhn = zeros(P, N);  % estimated state (AR coefficients)
atmp = aryule(Y(:), opt.ARorder);
atmp = atmp(:);
Xh = -atmp(2:end);
Xhn(:,P) = Xh;

En = zeros(N,1);        % error

% push the kalman gain and variances to steady state
n = 1+P;
C = Y(n-1:-1:n-P)';
for nn=1:25
  K = (Pt1step(:,:,n)*C')./(C*Pt1step(:,:,n)*C' + Sigma2V);
  Pt1step(:,:,n) = B*Ptt(:,:,n)*B' + Sigma2W;
  Ptt(:,:,n) = Pt1step(:,:,n) - K*C*Pt1step(:,:,n);
end


%% BIG LOOP
for n = 1+P:N

  if(length(opt.Sigma2V) == N)
    Sigma2V = opt.Sigma2V(n);
  end
  
  % Data stacked from t = n-p, to t = 1
  Ys = Y(n-1:-1:n-P);
  
  % Observation matrix
  C = Ys';
     
  % compute kalman gain
  K = (Pt1step(:,:,n)*C')./(C*Pt1step(:,:,n)*C' + Sigma2V);

  % one-step state error covariance
  Pt1step(:,:,n+1) = B*Ptt(:,:,n)*B' + Sigma2W;
   
  % Smoothing matrix
  if opt.ksmooth == 1
    Sm(:,:,n) = Ptt(:,:,n)*B'*inv(B*Pt1step(:,:,n+1)*B' + Sigma2W);
  end

  % full state error covariance
  Ptt(:,:,n+1) = Pt1step(:,:,n+1) - K*C*Pt1step(:,:,n+1);
 
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
  
  % limit the amplitude of the pole so the spectrum
  % does not get too flat.
    
  idbad = find(isfinite(Xh) == 0);
  if isempty(idbad) ~= 1 
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
    %sum(abs(Xhnsm(:,n+1) - Xhn(:,n+1)));
    Xhnsm(:,n) = Xhn(:,n) + Sm(:,:,n)*(Xhnsm(:,n+1) - Xhn(:,n));
    
    %compute the smoothed state variance as well
    Vn(:,:,n) = Ptt(:,:,n) + Sm(:,:,n)*(Vn(:,:,n+1) - Pt1step(:,:,n+1))*Sm(:,:,n)';
  end
  Xhn = Xhnsm;

  % RECALCULATE THE ERROR AFTER SMOOTHING
  for n = 1+P:N
    En(n) = Y(n) - Y(n-1:-1:n-P)'*Xhn(:,n);
  end
    
end

%% RESHAPE X(n) to get A(n)
An = single(Xhn);
Vn = single(Vn);

%% FILL IN THE INITIAL ESTIMATES
for k = 1:P+1
  An(:,k) = An(:,P+2);
  En(k) = En(P+2);
  Vn(:,:,k) = Vn(:,:,P+2);
end

