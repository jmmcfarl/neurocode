function LL_adj = adjustLL_STCBF(model,LL,nspks)


%%Add penalty terms from localization
Nmods = length(model.mods);
kern_l = length(model.mods(1).k);
sdim = model.mods(1).SDIM;
loc_penalty = zeros(Nmods,1);
for n = 1:Nmods
%     loc_penalty(n) = model.mods(n).locLambda*kernel_entropy(model.STCbasis,model.mods(n).STCcf,sdim);    
        loc_penalty(n) = model.mods(n).locLambda*kernel_std(model.STCbasis,model.mods(n).STCcf,sdim);

    %for prior covmat penalty on k's
%     cur_kern = model.STCbasis*model.mods(n).STCcf;
%     loc_penalty(n) = cur_kern'*model.mods(n).prior_precisionmat*cur_kern;
end

LL_adj = LL - sum(loc_penalty)/nspks;
