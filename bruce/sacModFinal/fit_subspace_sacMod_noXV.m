function [all_subMods] = fit_subspace_sacMod_noXV(stim_mod,Robs,Xmat,sac_Xmat,fit_inds,poss_d2T,post_sacMod,off_d2T,basemod_predrate)

NT = length(Robs);
silent = 1;
n_slags = size(sac_Xmat,2);
[~,~,~,~,filt_outs] = NMMmodel_eval(stim_mod,Robs,Xmat);

%%
cur_stim_params(1) = NMMcreate_stim_params(length(stim_mod.mods));
mod_signs = [1 1 1];
NL_types = {'lin','quad','quad'};
base_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types);
base_mod = NMMfit_filters(base_mod,Robs,filt_outs,[],fit_inds);

%%

X{1} = sac_Xmat;
X{2} = reshape(bsxfun(@times,sac_Xmat,reshape(filt_outs,NT,1,[])),NT,[]);

cur_stim_params(1) = NMMcreate_stim_params(n_slags);
cur_stim_params(2) = NMMcreate_stim_params([n_slags length(stim_mod.mods)]);

%initialize model with same number (and type) of subunits, but
%where the filter coefs can mix within the subspace
% mod_signs = [1 stim_mod.mods(:).sign];
% Xtargs = [1 2*ones(1,length(stim_mod.mods))];
% NL_types = cat(2,{'lin'},{stim_mod.mods(:).NLtype});
mod_signs = [1 1 1 1];
Xtargs = [1 2 2 2];
NL_types = {'lin','lin','quad','quad'};

clear all_sub_mods
for jj = 1:length(poss_d2T)
    modSeq_d2T = poss_d2T(jj);
    fprintf('subspace d2T %d of %d\n',jj,length(poss_d2T));
    
    %use opt d2T lambda from gain/offset model for the offset
    %kernel
    reg_params = NMMcreate_reg_params('lambda_d2T',[off_d2T repmat(modSeq_d2T,1,length(mod_signs)-1)]','boundary_conds',repmat([Inf Inf Inf],length(mod_signs),1));
    init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
    init_mod.mods(1).reg_params.boundary_conds = [Inf Inf Inf]; %loose boundary on offset term
    
    %set initial filters to have no mixing 
    for ii = 1:length(base_mod.mods)
        init_filt = repmat(base_mod.mods(ii).filtK',n_slags,1);
        init_mod.mods(ii+1).filtK = init_filt(:);
    end
    init_mod.spk_NL_params = stim_mod.spk_NL_params;
%     init_mod.mods(1).filtK = post_sacMod.mods(2).filtK; %set initial offset filter
    subspace_mod = NMMfit_filters(init_mod,Robs,X,[],fit_inds,silent);
    subspace_mod = NMMfit_logexp_spkNL(subspace_mod,Robs,X,[],fit_inds);
    
    [LL,~,pred_rate,~,~,~,nullLL] = NMMmodel_eval(subspace_mod,Robs(fit_inds),get_Xcell_tInds(X,fit_inds));
    subspace_mod.ovInfo = mean(pred_rate/mean(pred_rate).*log2(pred_rate/mean(pred_rate)));
    subspace_mod.ovLLimp = (LL-nullLL)/log(2);
    
    [sac_offset,sac_gain,sac_info,sac_LLimp] = deal(nan(n_slags,1));
    for ss = 1:n_slags
        temp = find(sac_Xmat(fit_inds,ss) == 1);
        cur_Robs = Robs(fit_inds(temp));
        
        rr = regress(pred_rate(temp),[ones(length(temp),1) basemod_predrate(fit_inds(temp))]);
        sac_offset(ss) = rr(1);
        sac_gain(ss) = rr(2);
        
        sac_info(ss) = nanmean(pred_rate(temp).*log2(pred_rate(temp)/mean(pred_rate(temp))))/mean(pred_rate(temp));
        
        %compute nullLL for data at this latency
        cur_nullLL = nansum(cur_Robs.*log2(mean(cur_Robs)) - mean(cur_Robs));
        cur_Nspks = sum(cur_Robs);
        
        sac_LLimp(ss) = (nansum(cur_Robs.*log2(pred_rate(temp)) - pred_rate(temp)) - cur_nullLL)/cur_Nspks;
    end
    
    subspace_mod.sac_offset = sac_offset;
    subspace_mod.sac_gain = sac_gain;
    subspace_mod.sac_modinfo = sac_info;
    subspace_mod.sac_LLimp = sac_LLimp;
    
    all_subMods{jj} = subspace_mod;
end

