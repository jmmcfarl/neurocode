function [LL, LLgrad] = FULLBF_temponly_LLinternal(params, Robs, T, model)
%
% Usage: [LL, LLgrad] = RG_LLinternal_elog( params, Robs, X, Skb, model, targets, lamrange, fprimes )
%

Nmods = length(model.mods);
[stimlen,wklen] = size(T);

hlen = 1; %no psc term
flen = 1;

b = params(end);
g = zeros(stimlen,1);
g = g + T*params(1:end-1);

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

% calculate gradient for temporal params
LLgrad(1:end-1) = residual*T;

%%****************************

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

