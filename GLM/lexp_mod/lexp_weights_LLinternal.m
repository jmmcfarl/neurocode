function [LL, LLgrad] = lexp_weights_LLinternal(params, Robs, X, model)
%
% Usage: [LL, LLgrad] = RG_LLinternal_elog( params, Robs, X, Skb, model, targets, lamrange, fprimes )
%

Nmods = length(model.mods);
stimlen = size(X,1);
SDIM = model.mods(1).SDIM;
klen = length(model.mods(1).k);
flen = klen/SDIM;

b = params(end);
k = params(1:end-1);
g = X*k;

g(g > 100) = 100;
expg = exp(g + b)';
if strcmp(model.spk_nl,'logexp')
    r = log(1+expg);
elseif strcmp(model.spk_nl,'exp')
    r = expg;
else
    error('invalid spk nl');
end
r(r < 1e-20) = 1e-20; %minimum predicted rate

LL = sum(Robs(flen:end) .* log(r(flen:end)) - r(flen:end));

if strcmp(model.spk_nl,'logexp')
    residual = (Robs./r - 1) .* expg ./ (1+expg);
else
    residual = Robs - r;
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);

% Calculate output of derivative module
for n = 1:Nmods
    LLgrad(n) = residual(flen:end) * X(flen:end,n);
end
%%****************************

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

