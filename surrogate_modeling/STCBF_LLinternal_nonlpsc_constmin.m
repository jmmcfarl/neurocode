function [LL, best_const] = STCBF_LLinternal_nonlpsc_constmin(params, Robs, Skb, model)
%
% Usage: [LL, LLgrad] = RG_LLinternal_elog( params, Robs, X, Skb, model, targets, lamrange, fprimes )
%
Nmods = length(model.mods);
[NT,NSTC_dims] = size(Skb);
hlen = 1; %no psc term

g = zeros(NT,1);
for n = 1:Nmods
    % Calculate convolution with internal receptive field
    gint{n} = Skb * params((n-1)*NSTC_dims + (1:NSTC_dims));
    
    %pass through current internal NL
    %         fgint = nlin_proc_stim(gint{n},model.mods(n).nly,model.mods(n).nlx);
    %using threshold nonlinearity
    fgint = gint{n};
    fgint(fgint < 0) = 0;
    
    %multiply by model weight
    g = g + fgint*model.mods(n).w; %convoles the current modules f(g) with its PSC term
end

opts.Display = 'off';
opts.MaxFunEvals = 20;
opts.MaxIter = 100;
opts.progTol = 1e-4;
opts.optTol = 1e-3;
[best_const LL eflag] = minFunc( @(K) helper_LL(K,g,Robs), 0, opts );

%%Add penalty terms from localization
kern_l = length(model.mods(1).k);
sdim = model.mods(1).SDIM;
kern_t = kern_l/sdim;
loc_penalty = zeros(Nmods,1);
for n = 1:Nmods
    cur_cfs = params((n-1)*NSTC_dims + (1:NSTC_dims));
    loc_penalty(n) = model.mods(n).locLambda*kernel_std(model.STCbasis,cur_cfs,sdim);
end

LL = LL + sum(loc_penalty);

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