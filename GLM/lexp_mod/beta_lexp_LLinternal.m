function [LL] = beta_lexp_LLinternal( params, Robs, g_mat, model)

% function [LL, LLgrad] = beta_lexp_LLinternal( params, Robs, g_mat, model)


Nmods = length(model.mods);
[NT,NSTC_dims] = size(g_mat);
hlen = 1; %no psc term

b = params(end);

g = zeros(NT,1);
for n = 1:Nmods
    % Calculate convolution with internal receptive field
    fgint = 1/params(n)*log(1+exp(params(n)*g_mat(:,n)));
    
    %multiply by model weight
    g = g + fgint*model.mods(n).w; %convoles the current modules f(g) with its PSC term
end

kx = g+b;
kx(kx > 50) = 50;
expg = exp(kx)';
if strcmp(model.spk_nl,'logexp')
    r = log(1+expg);
elseif strcmp(model.spk_nl,'exp')
    r = expg;
end
r(r < 1e-20) = 1e-20;

LL = sum(Robs .* log(r) - r);

residual = (Robs./r - 1) .* expg ./ (1+expg);

LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);

% Calculate output of derivative module
for n = 1:Nmods    
    temp = 1/params(n)*g_mat(:,n).*exp(params(n)*g_mat(:,n))./(1+exp(params(n)*g_mat(:,n)));
    temp = temp - 1/params(n)^2*log(1+exp(params(n)*g_mat(:,n)));
    LLgrad(n) = residual*temp;
end

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

%disp(sprintf( '%f\t', [LL params] ))

