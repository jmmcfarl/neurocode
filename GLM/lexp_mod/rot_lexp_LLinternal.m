function [LL, LLgrad] = rot_lexp_LLinternal( params, Robs, Skb, model,fprimes)


Nmods = length(model.mods);
[NT,NSTC_dims] = size(Skb);
hlen = 1; %no psc term

b = params(end);

g = zeros(NT,1);
min_x = min(model.mods(1).nlx);

for n = 1:Nmods
    % Calculate convolution with internal receptive field
    gint{n} = Skb * params((n-1)*NSTC_dims + (1:NSTC_dims));
%     gint{n}(gint{n}>100) = 100;

    if strcmp(model.mods(n).nltype,'lexp')
    bgint = gint{n}*model.mods(n).beta;
    too_large = find(bgint > 50);
    fgint = 1/model.mods(n).beta*log(1+exp(bgint)) - 1/model.mods(n).beta*log(1+exp(model.mods(n).beta*min_x));
    fgint(too_large) = 1/model.mods(n).beta*bgint(too_large) - 1/model.mods(n).beta*min_x;
    elseif strcmp(model.mods(n).nltype,'lin')
        fgint = gint{n};
    elseif strcmp(model.mods(n).nltype,'quad')
        fgint = gint{n}.^2;
    elseif strcmp(model.mods(n).nltype,'threshlin')
        fgint = gint{n};
        fgint(fgint < 0) = 0;
    elseif strcmp(model.mods(n).nltype,'rquad')
        fgint = gint{n}.^2;
        fgint(gint{n} < 0) = 0;
    elseif strcmp(model.mods(n).nltype,'uncon')
        %pass through current internal NL
        fgint = nlin_proc_stim(gint{n},model.mods(n).nly,model.mods(n).nlx);
    else
        error('unsupported NL')
    end
    %     fgint = 1/model.mods(n).beta*log(1+exp(bgint));
    %     fgint(too_large) = 1/model.mods(n).beta*bgint(too_large);
    %multiply by model weight
    g = g + fgint*model.mods(n).w; %convoles the current modules f(g) with its PSC term
end

kx = g+b;
kx(kx > 100) = 100;
expg = exp(kx)';
if strcmp(model.spk_nl,'logexp')
    r = log(1+expg);
elseif strcmp(model.spk_nl,'exp')
    r = expg;
end
r(r < 1e-20) = 1e-20;

LL = sum(Robs .* log(r) - r);

if strcmp(model.spk_nl,'logexp')
    residual = (Robs./r - 1) .* expg ./ (1+expg);
elseif strcmp(model.spk_nl,'exp')
    residual = Robs-r;
end

LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);

% Calculate output of derivative module
for n = 1:Nmods
    if strcmp(model.mods(n).nltype,'lexp')
    bgint = gint{n}*model.mods(n).beta;
    g = exp(bgint)./(1+exp(bgint));
    too_large = find(bgint > 50);
    g(too_large) = 1;
    elseif strcmp(model.mods(n).nltype,'lin')
        g = ones(size(gint{n}));
    elseif strcmp(model.mods(n).nltype,'quad')
        g = 2*gint{n};
    elseif strcmp(model.mods(n).nltype,'threshlin')
        g = ones(size(gint{n}));
        g(gint{n} < 0) = 0;
    elseif strcmp(model.mods(n).nltype,'rquad')
        g = 2*gint{n};
        g(gint{n} < 0) = 0;
    elseif strcmp(model.mods(n).nltype,'uncon')
        g = fprime(gint{n}, fprimes{n}, model.mods(n).nlx);
    end
    temp = Skb.*repmat(g,1,NSTC_dims)*model.mods(n).w;
    LLgrad(((n-1)*NSTC_dims+1):n*NSTC_dims) = (residual * temp)';
end

% loc_pens = arrayfun(@(x) x.locLambda,model.mods);
% if max(loc_pens) > 0
%     if strcmp(model.image_type,'2d')
%         error('loc penalty not working for 2d yet')
%     end
%     kern_l = length(model.mods(1).k);
%     sdim = model.mods(1).SDIM;
%     kern_t = kern_l/sdim;
%     loc_penalty = zeros(Nmods,1);
%     [smooth_x,smooth_y] = meshgrid(-3:3,-3:3);
%     smooth_sigma = 1;
%     smooth_kern = exp(-(smooth_x.^2+smooth_y.^2)/(2*smooth_sigma^2));
%     [X,T] = meshgrid(1:sdim,1:kern_t);
%     for n = 1:Nmods
%         cur_cfs = params((n-1)*NSTC_dims + (1:NSTC_dims));
%         cur_kmat = reshape(model.STCbasis*cur_cfs,kern_t,sdim);
%         cur_smooth_kmat = conv2(cur_kmat.^2,smooth_kern,'same');
%         [~,peakloc] = max(cur_smooth_kmat(:));
%         [peakt,peakx] = ind2sub(size(cur_smooth_kmat),peakloc);
%         distmat = exp((X-peakx).^2/(2*model.mods(n).locSigmaX^2) + (T-peakt).^2/(2*model.mods(n).locSigmaT^2));
%         l2_pen{n} = model.mods(n).locLambda*distmat(:);       
%         loc_penalty(n) = model.mods(n).locLambda*sum(cur_kmat(:).^2.*l2_pen{n});
%     end
%     
%     LL = LL - sum(loc_penalty);
% end
% 
% eps = 1e-8;
% LLgrad_pen = zeros(size(LLgrad));
% if max(loc_pens) > 0
%     for n = 1:Nmods
%         cur_cfs = params((n-1)*NSTC_dims + (1:NSTC_dims));
%         for j = 1:length(cur_cfs)
%             cur_cfs_eps = cur_cfs;
%             cur_cfs_eps(j) = cur_cfs_eps(j) + eps;
%             cur_kmat = reshape(model.STCbasis*cur_cfs_eps,kern_t,sdim);
%             loc_pen_eps = model.mods(n).locLambda*sum(cur_kmat(:).^2.*l2_pen{n});
%             LLgrad_pen((n-1)*NSTC_dims+j) = (loc_pen_eps-loc_penalty(n))/eps;
%         end
%     end
% end

%%Add penalty terms from localization
kern_l = length(model.mods(1).k);
sdim = model.mods(1).SDIM;
kern_t = kern_l/sdim;
loc_penalty = zeros(Nmods,1);
for n = 1:Nmods
    cur_cfs = params((n-1)*NSTC_dims + (1:NSTC_dims));
    loc_penalty(n) = model.mods(n).locLambda*kernel_std(model.STCbasis,cur_cfs,sdim);
end

LL = LL - sum(loc_penalty);


eps = 1e-8;
LLgrad_pen = zeros(size(LLgrad));
for n = 1:Nmods
    cur_cfs = params((n-1)*NSTC_dims + (1:NSTC_dims));   
    %     for numerical gradient calculation
    for j = 1:model.STCdim
        cur_cfs_eps = cur_cfs;
        cur_cfs_eps(j) = cur_cfs_eps(j) + eps;
        ent_dx = kernel_std(model.STCbasis,cur_cfs_eps,sdim);
        ent = kernel_std(model.STCbasis,cur_cfs,sdim);
        LLgrad_pen((n-1)*NSTC_dims+j) = model.mods(n).locLambda*(ent_dx-ent)/eps;
    end
    
end

LLgrad = LLgrad - LLgrad_pen;

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

%disp(sprintf( '%f\t', [LL params] ))

