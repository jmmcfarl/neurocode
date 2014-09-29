function [subspace_mod,predrate] = fit_subspace_sacMod(stim_mod,Robs,Xmat,sac_Xmat,tr_inds,xv_inds,poss_d2T,post_sacMod)

NT = length(Robs);
silent = 1;
n_slags = size(sac_Xmat,2);
[~,~,~,~,filt_outs] = NMMmodel_eval(stim_mod,Robs,Xmat);
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

offset_d2T = post_sacMod.opt_d2T_off;

%cycle over a range of possible d2T reg lambdas
subspace_xvLL = nan(length(poss_d2T),1);
null_prate = mean(Robs(tr_inds));
null_xvLL = sum(Robs(xv_inds).*log(ones(size(xv_inds))*null_prate) - ones(size(xv_inds))*null_prate)/sum(Robs(xv_inds));

clear all_sub_mods
for jj = 1:length(poss_d2T)
    modSeq_d2T = poss_d2T(jj);
    fprintf('Xval %d of %d\n',jj,length(poss_d2T));
    %use opt d2T lambda from gain/offset model for the offset
    %kernel
    reg_params = NMMcreate_reg_params('lambda_d2T',[offset_d2T repmat(modSeq_d2T,1,length(mod_signs)-1)]','boundary_conds',repmat([Inf Inf Inf],length(mod_signs),1));
    init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
    init_mod.mods(1).reg_params.boundary_conds = [0 Inf Inf]; %0 boundary on offset term
    
%     %set initial filters to have no mixing 
%     for ii = 1:length(stim_mod.mods)
%         init_filt = zeros(n_slags,length(stim_mod.mods));
%         init_filt(:,ii) = 1;
%         init_mod.mods(ii+1).filtK = init_filt(:);
%     end
    init_mod.spk_NL_params = stim_mod.spk_NL_params;
    init_mod.mods(1).filtK = post_sacMod.mods(2).filtK; %set initial offset filter
    subspace_mod = NMMfit_filters(init_mod,Robs,X,[],tr_inds,silent);
    subspace_xvLL(jj) = NMMmodel_eval(subspace_mod,Robs(xv_inds),get_Xcell_tInds(X,xv_inds));
    all_sub_mods(jj) = subspace_mod;
end
%select best lambda using xvLL
[~,optloc] = max(subspace_xvLL);
subspace_optd2T = poss_d2T(optloc);

%%
all_inds = union(tr_inds,xv_inds);

%now initialize model with optimal d2T lambda and fit to all
%used data
reg_params = NMMcreate_reg_params('lambda_d2T',[offset_d2T repmat(subspace_optd2T,1,length(mod_signs)-1)]','boundary_conds',repmat([Inf Inf Inf],length(mod_signs),1));
init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
init_mod.spk_NL_params = stim_mod.spk_NL_params;
init_mod.mods(1).reg_params.boundary_conds = [0 Inf Inf];
% for ii = 1:length(init_mod.mods)-1
%     init_filt = zeros(n_slags,length(stim_mod.mods));
%     init_filt(:,ii) = 1;
%     init_mod.mods(ii+1).filtK = init_filt(:);
% end
init_mod.mods(1).filtK = post_sacMod.mods(2).filtK;
init_mod.spk_NL_params(1) = stim_mod.spk_NL_params(1);
subspace_mod = NMMfit_filters(init_mod,Robs,X,[],all_inds,silent);

%evaluate overall info on anysac inds
[LL,~,predrate,~,~,~,nullLL] = NMMmodel_eval(subspace_mod,Robs(all_inds),get_Xcell_tInds(X,all_inds));
subspace_mod.ovInfo = mean(predrate/mean(predrate).*log2(predrate/mean(predrate)));
subspace_mod.ovLLimp = (LL-nullLL)/log(2);

subspace_mod.opt_L2 = subspace_optd2T;
subspace_mod.xvLLimp = (subspace_xvLL - null_xvLL)/log(2);

[~,~,predrate] = NMMmodel_eval(subspace_mod,Robs,X);
