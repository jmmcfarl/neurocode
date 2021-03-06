function [LL, LLgrad] = fitGNM_filter_scales_internal(params, Robs, X, model, fprimes)

Nmods = length(model.mods);
[stimlen,klen] = size(X);
flen = model.stim_params.flen;
fsdim = model.stim_params.fsdim;
sdim = model.stim_params.sdim;

b = params(end);
g = zeros(stimlen,1);
kerns = zeros(Nmods,klen);
for n = 1:Nmods
    kerns(n,:) = model.mods(n).k;
    % Calculate convolution with internal receptive field
    gint{n} = X * kerns(n,:)'*params(n);
    
    if strcmp(model.mods(n).nltype,'lexp')
        min_x = min(model.mods(1).nlx);
        bgint = (gint{n}-model.mods(n).lexp_theta)*model.mods(n).lexp_beta;
        too_large = bgint > 50;
        fgint = 1/model.mods(n).lexp_beta*log(1+exp(bgint)) - 1/model.mods(n).lexp_beta*log(1+exp(model.mods(n).lexp_beta*min_x));
        fgint(too_large) = 1/model.mods(n).lexp_beta*bgint(too_large) - 1/model.mods(n).lexp_beta*min_x;
    elseif strcmp(model.mods(n).nltype,'uncon')
        fgint = nlin_proc_stim(gint{n},model.mods(n).nly,model.mods(n).nlx);
    elseif strcmp(model.mods(n).nltype,'quad')
        fgint = gint{n}.^2;
    elseif strcmp(model.mods(n).nltype,'lin')
        fgint = gint{n};
    elseif strcmp(model.mods(n).nltype,'threshlin')
        fgint = gint{n};
        fgint(fgint < 0) = 0;
    elseif strcmp(model.mods(n).nltype,'rquad')
        fgint = gint{n}.^2;
        fgint(gint{n} < 0) = 0;
    else
        error('Invalid internal NL');
    end
    
    %multiply by weight and add to generating function
    g = g + fgint*model.mods(n).w;
end
g = g + b;
if max(g) > 100
    g(g > 100) = 100;%max g value to prevent overflow
    disp('Warning, overflow in g')
end

if strcmp(model.spk_nl,'logexp')
    expg = exp(model.spk_beta*(g-model.spk_theta))';
    r = model.spk_alpha*log(1+expg);
elseif strcmp(model.spk_nl,'exp')
    expg = exp(g)';
    r = expg;
else
    error('invalid spk nl');
end
if min(r) < 1e-20
    r(r < 1e-20) = 1e-20; %minimum predicted rate
    disp('Warning, underflow in predicted r')
end

LL = sum(Robs.* log(r) - r);

if strcmp(model.spk_nl,'logexp')
    residual = model.spk_alpha*model.spk_beta*(Robs./r - 1) .* expg ./ (1+expg);
elseif strcmp(model.spk_nl,'exp')
    residual = Robs - r;
else
    error('Unsupported spiking NL')
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);

% Calculate output of derivative module
chunk_size = 1000; %maximum chunk size for processing
for n = 1:Nmods
    
    if strcmp(model.mods(n).nltype,'lexp')
        %internal gen fun processed by f'
        bgint = (gint{n}-model.mods(n).lexp_theta)*model.mods(n).lexp_beta;
        g = exp(bgint)./(1+exp(bgint));
        too_large = bgint > 50;
        g(too_large) = 1;
    elseif strcmp(model.mods(n).nltype,'uncon')
        g = fprime(gint{n},fprimes{n}, model.mods(n).nlx);
    elseif strcmp(model.mods(n).nltype,'quad')
        g = 2*gint{n};
    elseif strcmp(model.mods(n).nltype,'lin')
        g = ones(size(gint{n}));
    elseif strcmp(model.mods(n).nltype,'threshlin')
        g = ones(size(gint{n}));
        g(gint{n} < 0) = 0;
    elseif strcmp(model.mods(n).nltype,'rquad')
        g = 2*gint{n};
        g(gint{n} < 0) = 0;
    end
    
    temp = X*kerns(n,:)'.*g*model.mods(n).w;
    LLgrad(n) = residual*temp;
end


Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

