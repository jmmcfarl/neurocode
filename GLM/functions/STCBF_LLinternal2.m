function [LL, LLgrad] = STCBF_LLinternal2( params, Robs, Skb, model, lamrange, fprimes, pscs)
%
% Usage: [LL, LLgrad] = RG_LLinternal_elog( params, Robs, X, Skb, model, targets, lamrange, fprimes, pscs )
%

pos_dims = find(model.dim_signs == 1);
neg_dims = find(model.dim_signs == -1);

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

[NT,NSTC_dims] = size(Skb);
hlen = length(model.mods(1).h);

b = params(end);

g = zeros(NT,1);
for n = 1:Nmods
    % Calculate convolution with internal receptive field
    gint{n} = Skb * params((n-1)*NSTC_dims + (1:NSTC_dims));
   
    %pass through current internal NL
    fgint = nlin_proc_stim(gint{n},model.mods(n).nly,model.mods(n).nlx);

    %convolve with current PSC term
    g = g + g_convolve_jmm(fgint,pscs{n}); %convoles the current modules f(g) with its PSC term
end
g(1:hlen-1) = [];

expg = exp(g + b)';
r = log(1+expg);
r(r < 1e-20) = 1e-20;

LL = sum(Robs .* log(r) - r);

residual = (Robs./r - 1) .* expg ./ (1+expg);

LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);

% Calculate output of derivative module
for n = 1:Nmods
    g = fprime(gint{n}, fprimes{n}, model.mods(n).nlx);
    if model.mod_signs(n) == 1 %if this is a positive module
        for m = 1:length(pos_dims)
           temp = g_convolve_jmm(Skb(:,pos_dims(m)).*g, pscs{n});
           LLgrad((n-1)*NSTC_dims + pos_dims(m)) = residual * temp(hlen:end);
        end
    elseif model.mod_signs(n) == -1 %if this is a negative module
        for m = 1:length(neg_dims)
            temp = g_convolve_jmm(Skb(:,neg_dims(m)).*g, pscs{n});
            LLgrad((n-1)*NSTC_dims + neg_dims(m)) = residual * temp(hlen:end);
        end
    else
        for m = 1:NSTC_dims
            temp = g_convolve_jmm(Skb(:,m).*g, pscs{n});
            LLgrad((n-1)*NSTC_dims + m) = residual * temp(hlen:end);
        end
    end
end

%%Add penalty terms from localization
kern_l = length(model.mods(1).k);
sdim = model.mods(1).SDIM;
kern_t = kern_l/sdim;
loc_penalty = zeros(Nmods,1);
for n = 1:Nmods
    cur_cfs = params((n-1)*NSTC_dims + (1:NSTC_dims));
    
    %     loc_penalty(n) = model.mods(n).locLambda*kernel_entropy(model.STCbasis,cur_cfs,sdim);
        loc_penalty(n) = model.mods(n).locLambda*kernel_std(model.STCbasis,cur_cfs,sdim);
    
%     %for prior covmat penalty on k's
%     cur_kern = model.STCbasis*cur_cfs;
%     loc_penalty(n) = cur_kern'*model.mods(n).prior_precisionmat*cur_kern;
end

LL = LL - sum(loc_penalty);


%%STILL PROBLEMS HERE!!
eps = 1e-8;
zerovec = zeros(1,model.STCdim)';
LLgrad_pen = zeros(size(LLgrad));
for n = 1:Nmods
    cur_cfs = params((n-1)*NSTC_dims + (1:NSTC_dims));
    
    %     for numerical gradient calculation
    if model.mod_signs(n) == 1 %if positive module
        for j = 1:length(pos_dims)
            cur_cfs_eps = cur_cfs;
            cur_cfs_eps(pos_dims(j)) = cur_cfs_eps(pos_dims(j)) + eps;
            ent_dx = kernel_std(model.STCbasis,cur_cfs_eps,sdim);
            ent = kernel_std(model.STCbasis,cur_cfs,sdim);
            LLgrad_pen((n-1)*NSTC_dims+pos_dims(j)) = model.mods(n).locLambda*(ent_dx-ent)/eps;
        end
    elseif model.mod_signs(n) == -1 %if negative module
        for j = 1:length(neg_dims)
            cur_cfs_eps = cur_cfs;
            cur_cfs_eps(neg_dims(j)) = cur_cfs_eps(neg_dims(j)) + eps;
            ent_dx = kernel_std(model.STCbasis,cur_cfs_eps,sdim);
            ent = kernel_std(model.STCbasis,cur_cfs,sdim);
            LLgrad_pen((n-1)*NSTC_dims+neg_dims(j)) = model.mods(n).locLambda*(ent_dx-ent)/eps;
        end
    else
        for j = 1:model.STCdim
            cur_cfs_eps = cur_cfs;
            cur_cfs_eps(j) = cur_cfs_eps(j) + eps;
            ent_dx = kernel_std(model.STCbasis,cur_cfs_eps,sdim);
            ent = kernel_std(model.STCbasis,cur_cfs,sdim);
            LLgrad_pen((n-1)*NSTC_dims+j) = model.mods(n).locLambda*(ent_dx-ent)/eps;
        end
        
    end
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

