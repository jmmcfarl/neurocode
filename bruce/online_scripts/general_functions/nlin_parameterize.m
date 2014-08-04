function [theta,beta,alpha,nlin] = nlin_parameterize( gd, lm, gpar )
%
% Usage: [theta,beta,alpha,nlin] = nlin_parameterize( gd, lm, <gpar> )
%
% Fit f(g) = alpha/beta * log( 1 + exp(beta*(g+theta)) );

if nargin < 3
  GMAX = 6;
  GINC = 0.05;
else
  GMAX = gpar(1);
  GINC = gpar(2);
end

ga = -GMAX:GINC:GMAX;
gd = make_row(gd);
lm = make_row(lm);
% Initial conditions
[dum alpha0 theta0 x0 x1] = nlin_totfit( gd, lm );

if (alpha0 <= 0) | (theta0 < length(ga)/2)
  nlin0 = calc_nlin2(gd,lm);
  %theta0 = min(find(nlin0 > 0));
  theta0 = floor(0.5*(x0+x1));
  [alpha0 b] = linfit(theta0:length(ga),nlin0(theta0:end),gd(theta0:end));
end  

beta0 = 4;

% plot(alpha0/beta0/0.05*log(1+exp(beta0*((1:241)-theta0)*0.05)),'c')

% Minimize error 
initial = [theta0 beta0*0.05 alpha0];
opts = optimset('GradObj','on','LargeScale','off','Display','off','MaxIter',2000,'MaxFunEvals',10000);

[model chi2] = fminunc( @(K) min_nlin(K,gd,lm), initial,opts );

theta = (model(1)-floor(length(ga)/2))*GINC;
beta = model(2)/GINC;
alpha = model(3)*beta;

nlin = alpha/beta*log(1+exp(beta*(ga-theta)));