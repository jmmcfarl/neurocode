function [LL, best_const] = STCBF_LLinternal_nopsc_constmin(params, Robs, Skb, model)
%
% Usage: [LL, LLgrad] = RG_LLinternal_elog( params, Robs, X, Skb, model, targets, lamrange, fprimes )
%

Nmods = length(model.mods);
[NT,NSTC_dims] = size(Skb);
hlen = 1; %no psc term

kern_l = length(model.mods(1).k);
sdim = model.mods(1).SDIM;
kern_t = kern_l/sdim;
loc_penalty = zeros(Nmods,1);

nl_penalty = zeros(Nmods,1);

g = zeros(NT,1);
for n = 1:Nmods
    cur_xmod_inds = (n-1)*NSTC_dims + (1:NSTC_dims);
    cur_nlmod_inds = (n-1)*11 + (1:11);
    
    % Calculate convolution with internal receptive field
    gint = Skb * params{1}(cur_xmod_inds);
    
    %     pass through current internal NL
    fgint = nlin_proc_stim(gint,params{2}(cur_nlmod_inds),model.mods(n).nlx);
    
    %multiply by model weight
    g = g + fgint*model.mods(n).w; %convoles the current modules f(g) with its PSC term
    
    loc_penalty(n) = model.mods(n).locLambda*kernel_std(model.STCbasis,params{1}(cur_xmod_inds),sdim);

    chunk = params{2}(cur_nlmod_inds);
    nl_penalty(n) = model.mods(n).lnl2*sum((2*chunk(2:end-1) - chunk(1:end-2) - chunk(3:end)).^2);
    
end

opts.Display = 'off';
opts.MaxFunEvals = 20;
opts.MaxIter = 100;
opts.progTol = 1e-4;
opts.optTol = 1e-3;
[best_const LL eflag] = minFunc( @(K) helper_LL(K,g,Robs), 0, opts );


LL = LL + sum(loc_penalty) + sum(nl_penalty);

Nspks = sum(Robs);
LL = LL/Nspks;

end

function [LL,LLgrad] = helper_LL(const,g,Robs)

expg = exp(g+const)';
r = log(1+expg);
LL = sum(Robs.*log(r)-r);
residual =(Robs./r-1).*expg./(1+expg);
LLgrad = sum(residual);
LL = -LL;
LLgrad = -LLgrad;

end