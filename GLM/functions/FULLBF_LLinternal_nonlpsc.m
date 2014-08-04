function [LL, LLgrad] = FULLBF_LLinternal_nonlpsc(params, Robs, X, model, lamrange,targets)
%
% Usage: [LL, LLgrad] = RG_LLinternal_elog( params, Robs, X, Skb, model, targets, lamrange, fprimes )
%

if isempty(lamrange)
    els = [];  r1s = [];
else
    els = lamrange(:,1);
    if sum(els) > 0
        r1s = lamrange(:,2);
        r2s = lamrange(:,3);
    else
        els = [];  r1s = [];
    end
end

Nmods = length(model.mods);
[NT,klen] = size(X);
hlen = 1; %no psc term

b = params(end);

all_mods = 1:Nmods;
fixed_mods = setdiff(all_mods,targets);

g = zeros(NT,1);
for n = 1:length(targets)
    % Calculate convolution with internal receptive field
    gint{n} = X * params((n-1)*klen + (1:klen));   
    fgint = gint{n};
    fgint(fgint < 0) = 0;
    g = g + fgint*model.mods(targets(n)).w; %convoles the current modules f(g) with its PSC term
end
for n = 1:length(fixed_mods)
    % Calculate convolution with internal receptive field
    fgint = X * model.mods(fixed_mods(n)).k;   
    fgint(fgint < 0) = 0;   
    g = g + fgint*model.mods(fixed_mods(n)).w; %convoles the current modules f(g) with its PSC term
end

expg = exp(g + b)';
r = log(1+expg);
r(r < 1e-20) = 1e-20;

LL = sum(Robs .* log(r) - r);

residual = (Robs./r - 1) .* expg ./ (1+expg);

LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);

% Calculate output of derivative module
for n = 1:length(targets)
    g = zeros(size(gint{targets(n)}));
    g(gint{n} >= 0) = 1;
%     g = fprime(gint{n}, fprimes{n}, model.mods(n).nlx);

    temp = X.*repmat(g,1,klen)*model.mods(targets(n)).w;
    LLgrad(((n-1)*klen+1):n*klen) = (residual * temp)';
%     for m = 1:NSTC_dims
%         temp = X(:,m).*g*model.mods(n).w;
%         LLgrad((n-1)*NSTC_dims + m) = residual * temp(hlen:end);
%     end
end

%%Add penalty terms from localization
kern_l = length(model.mods(1).k);
sdim = model.mods(1).SDIM;
kern_t = kern_l/sdim;
loc_penalty = zeros(length(targets),1);
for n = 1:length(targets)
    cur_kern = params((n-1)*klen + (1:klen));
    kern_mat{n} = reshape(cur_kern,kern_t,sdim);
    %     spatial_dist = var(reshape(cur_ks,kern_t,sdim));
    %     spatial_dist = spatial_dist/sum(spatial_dist);
    %     x_ax = 1:sdim;
    %     com = x_ax*spatial_dist';
    %     kstd = (x_ax-com).^2*spatial_dist';
    %     loc_penalty(n) = model.mods(n).locLambda*kstd;
    %         smooth_penalty(n) = model.dist_pen*cur_kern'*model.prior_precisionmat*cur_kern;
    
    % lapl = del2(kern_mat{n});
    % smooth_penalty(n) = model.dist_pen*sum(lapl(:).^2);

    %**THIS WORKS
    space_derivs = (kern_mat{n}(:,2:end) - kern_mat{n}(:,1:end-1));
    time_derivs = (kern_mat{n}(2:end,:) - kern_mat{n}(1:end-1,:));
    smooth_penalty(n) = model.dist_pen*(sum(space_derivs(:).^2) + sum(time_derivs(:).^2));
    %**
    
    %          smooth_penalty(n) = model.dist_pen*cur_kern'*model.prior_precisionmat*cur_kern;
    
    %     temp = cur_kern' / model.R;
    %     smooth_penalty(n) = sum(temp.^2);
    %     Y = model.V*cur_kern;
    %     smooth_penalty(n) = Y'*model.Dinv*Y;
end

LL = LL - sum(smooth_penalty);


%%STILL PROBLEMS HERE!!
eps = 1e-8;
zerovec = zeros(1,kern_l)';
LLgrad_pen = zeros(size(LLgrad));
for n = 1:length(targets)
    cur_cfs = params((n-1)*kern_l + (1:kern_l));
    
    %     for numerical gradient calculation
    %     for j = 1:kern_l
    % %         cur_cfs_eps = cur_cfs;
    % %         cur_cfs_eps(j) = cur_cfs_eps(j) + eps;
    % %         y
    % %         ent_dx = kernel_std(model.STCbasis,cur_cfs_eps,sdim);
    % %         ent = kernel_std(model.STCbasis,cur_cfs,sdim);
    % %         LLgrad_pen((n-1)*NSTC_dims+j) = model.mods(n).locLambda*(ent_dx-ent)/eps;
    %     end
    
%     lapl = del2(kern_mat{n});
%     LLgrad_pen((n-1)*kern_l + (1:kern_l)) = -2*model.dist_pen*lapl(:);

    
[gradmat,sp_gradmat,ti_gradmat] = deal(zeros(size(kern_mat{n})));
sp_gradmat(:,2:end) = kern_mat{n}(:,2:end) - kern_mat{n}(:,1:end-1);
sp_gradmat(:,1:end-1) = sp_gradmat(:,1:end-1) + kern_mat{n}(:,1:end-1) - kern_mat{n}(:,2:end);
ti_gradmat(2:end,:) = kern_mat{n}(2:end,:) - kern_mat{n}(1:end-1,:);
ti_gradmat(1:end-1,:) = ti_gradmat(1:end-1,:) + kern_mat{n}(1:end-1,:) - kern_mat{n}(2:end,:);
gradmat = sp_gradmat + ti_gradmat;  

    LLgrad_pen((n-1)*kern_l + (1:kern_l)) = 2*model.dist_pen*gradmat(:);

%     LLgrad_pen((n-1)*kern_l + (1:kern_l)) = 2*model.dist_pen*cur_kern'*model.prior_precisionmat;

%     Y = model.V*cur_cfs;
%     LLgrad_pen((n-1)*kern_l + (1:kern_l)) = 2*Y'*model.Dinv;
% 
end

LLgrad = LLgrad - LLgrad_pen;

% % Add penalty terms from slope
% for i = 1:length(els)
%   range = r1s(i):r2s(i);
%   chunk = params(range);
%   LL = LL - els(i) * sum((chunk(2:end) - chunk(1:end-1)).^2);
%   LLgrad(range(2:end)) = LLgrad(range(2:end)) - 2*els(i)*(chunk(2:end) - chunk(1:end-1));
%   LLgrad(range(1:end-1)) = LLgrad(range(1:end-1)) - 2*els(i)*(chunk(1:end-1) - chunk(2:end));
% end

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

%disp(sprintf( '%f\t', [LL params] ))

