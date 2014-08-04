function [LL, LLgrad] = gabor_weights_LLinternal(params, Robs, X)
%
% Usage: [LL, LLgrad] = RG_LLinternal_elog( params, Robs, X, Skb, model, targets, lamrange, fprimes )
%

[NT,klen] = size(X);
b = params(end);
k = params(1:end-1);
g = X*k;

too_large = find(g > 100);
expg = exp(g + b)';
r = log(1+expg);
r(too_large) = g(too_large);
r(r < 1e-20) = 1e-20; %minimum predicted rate

LL = sum(Robs.*log(r) - r);

residual = (Robs./r - 1) .* expg ./ (1+expg);
residual(too_large) = (Robs(too_large)./r(too_large) - 1);

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);

% Calculate output of derivative module
for n = 1:length(k)
    LLgrad(n) = residual * X(:,n);
end
%%****************************

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

