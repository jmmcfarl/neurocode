function scale_grad = compute_gnm_scalegrad(lambda_scale,gnm,X,Robs,imod,g_rest,cur_fprime,lpen)

cur_mod = gnm.mods(imod);
cur_k = cur_mod.k;
kx = X*cur_k;
if strcmp(cur_mod.nltype,'lexp')
    min_x = min(kx);
    cur_g = 1/cur_mod.lexp_beta*log(1+exp(cur_mod.lexp_beta*(kx-cur_mod.lexp_theta))) - ...
        1/cur_mod.lexp_beta*log(1+exp(cur_mod.lexp_beta*(min_x-cur_mod.lexp_theta)));
elseif strcmp(cur_mod.nltype,'quad')
    cur_g = kx.^2;
elseif strcmp(cur_mod.nltype,'lin')
    cur_g = kx;
elseif strcmp(cur_mod.nltype,'threshlin')
    cur_g = kx;
    cur_g(cur_g < 0) = 0;
elseif strcmp(cur_mod.nltype,'rquad')
    cur_g = kx.^2;
    cur_g(cur_g < 0) = 0;
elseif strcmp(cur_mod.nltype,'uncon')
    cur_g = nlin_proc_stim(kx,cur_mod.nly,cur_mod.nlx); %pass internal output through NL
else
    error('Invalid nl type');
end
g_tot = g_rest + cur_mod.w*cur_g;

if strcmp(gnm.spk_nl,'logexp')
    internal = gnm.spk_beta*g_tot;
    too_big = find(internal > 60);
    expg = exp(internal);
    r = gnm.spk_alpha*log(1+expg);
    r(too_big) = gnm.spk_alpha*internal(too_big);
elseif strcmp(gnm.spk_nl,'exp')
    if max(g_tot) > 100
        g_tot(g_tot > 100) = 100;%max g value to prevent overflow
        disp('Warning, overflow in g')
    end
    expg = exp(g_tot)';
    r = expg;
else
    error('invalid spk nl');
end
if min(r) < 1e-20
    r(r < 1e-20) = 1e-20; %minimum predicted rate
    disp('Warning, underflow in predicted r')
end

if strcmp(gnm.spk_nl,'logexp')
    residual = gnm.spk_alpha*gnm.spk_beta*(Robs./r - 1).* expg ./ (1+expg);
    residual(too_big) = gnm.spk_alpha*gnm.spk_beta*(Robs(too_big)./r(too_big)-1);
elseif strcmp(gnm.spk_nl,'exp')
    residual = Robs - r;
else
    error('Unsupported spiking NL')
end


if strcmp(cur_mod.nltype,'lexp')
    %internal gen fun processed by f'
    bgint = (kx-cur_mod.lexp_theta)*cur_mod.lexp_beta;
    g = exp(bgint)./(1+exp(bgint));
    too_large = bgint > 50;
    g(too_large) = 1;
elseif strcmp(cur_mod.nltype,'uncon')
    g = fprime(kx,cur_fprime, cur_mod.nlx);
elseif strcmp(cur_mod.nltype,'quad')
    g = 2*kx;
elseif strcmp(cur_mod.nltype,'lin')
    g = ones(size(kx));
elseif strcmp(cur_mod.nltype,'threshlin')
    g = ones(size(kx));
    g(kx < 0) = 0;
elseif strcmp(cur_mod.nltype,'rquad')
    g = 2*kx;
    g(kx < 0) = 0;
end

temp = kx.*g*cur_mod.w;
ll_grad = sum(residual.*temp);

pen_grad = lambda_scale*(lpen.l1x(imod) + ...
  2*lpen.lapl(imod) + 2*lpen.l2x(imod) + 2*lpen.laplXT(imod));

scale_grad = abs(ll_grad - pen_grad);
