function [sacMod,pred_rate] = sacMod_scan_doubleregularization(base_mod,Robs,Xsac,gain_sigs,tr_inds,xv_inds,poss_d2T_off,poss_d2T_gain)
%[sacMod,xvLLs,opt_d2T,opt_L2] = sacMod_scan_regularization(base_mod,Robs,Xsac,gain_sigs,tr_inds,xv_inds,poss_d2T,poss_L2)
%Fits a (post-filtering) saccade gain/offset model. Can use a number of
%independent gain filters specified by the number of gain signals in
%"gain_sigs". Also uses xvalLL to scan a range of d2T and L2 reg
%parameters.

n_slags = size(Xsac,2);
n_gains = size(gain_sigs,2);

tr_stim{1} = gain_sigs; 
tr_stim{2} = Xsac; %saccade timing indicator matrix
tr_stim{3} = reshape(bsxfun(@times,Xsac,reshape(gain_sigs,[],1,n_gains)), size(gain_sigs,1),[]);

sac_stim_params(1) = NMMcreate_stim_params(n_gains);
sac_stim_params(2) = NMMcreate_stim_params(n_slags);
sac_stim_params(3) = NMMcreate_stim_params([n_slags n_gains]);

mod_signs = [1 1 1];
Xtargets = [1 2 3];
NL_types = {'lin','lin','lin'};
sac_reg_params = NMMcreate_reg_params('boundary_conds',repmat([0 0 0],length(mod_signs),1));

%%

xvLLs = nan(length(poss_d2T_off),length(poss_d2T_gain));
null_prate = mean(Robs(tr_inds));
null_xvLL = sum(Robs(xv_inds).*log(ones(size(xv_inds))*null_prate) - ones(size(xv_inds))*null_prate)/sum(Robs(xv_inds));
if length(poss_d2T_gain) > 1 || length(poss_d2T_off) > 1
    
    for jj = 1:length(poss_d2T_off)
        for ii = 1:length(poss_d2T_gain)
            cur_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
            cur_mod.mods(1).filtK(:) = 1; %initialize base gains to 1
            cur_mod.spk_NL_params = base_mod.spk_NL_params;
            cur_mod = NMMadjust_regularization(cur_mod,[2],'lambda_d2T',poss_d2T_off(jj));
            cur_mod = NMMadjust_regularization(cur_mod,[3],'lambda_d2T',poss_d2T_gain(ii));
            cur_mod = NMMfit_filters(cur_mod,Robs,tr_stim,[],tr_inds,1,[],[],[2 3]);
            xvLLs(jj,ii) = NMMmodel_eval(cur_mod,Robs(xv_inds),get_Xcell_tInds(tr_stim,xv_inds));
        end
    end
    
    [~,optloc] = max(xvLLs(:));
    [optloc_x,optloc_y] = ind2sub([length(poss_d2T_off) length(poss_d2T_gain)],optloc);
    opt_d2T_off = poss_d2T_off(optloc_x);
    opt_d2T_gain = poss_d2T_gain(optloc_y);
else
    opt_d2T_off = poss_d2T_off;
    opt_d2T_gain = poss_d2T_gain;
end

%%
all_inds = union(tr_inds,xv_inds);
sacMod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
sacMod.mods(1).filtK(:) = 1; %set base gains to 1
sacMod.spk_NL_params = base_mod.spk_NL_params;
sacMod = NMMadjust_regularization(sacMod,[2],'lambda_d2T',opt_d2T_off);
sacMod = NMMadjust_regularization(sacMod,[3],'lambda_d2T',opt_d2T_gain);
sacMod = NMMfit_filters(sacMod,Robs,tr_stim,[],all_inds,1,[],[],[2 3]);
sacMod = NMMfit_logexp_spkNL(sacMod,Robs,tr_stim,[],all_inds);

%%
[LL,~,pred_rate,~,~,~,nullLL] = NMMmodel_eval(sacMod,Robs(all_inds),get_Xcell_tInds(tr_stim,all_inds));
sacMod.ovInfo = mean(pred_rate/mean(pred_rate).*log2(pred_rate/mean(pred_rate)));
sacMod.ovLLimp = (LL-nullLL)/log(2);
sacMod.nullLL = nullLL;

[~,~,pred_rate] = NMMmodel_eval(sacMod,Robs,tr_stim);

%%
sacMod.opt_d2T_off = opt_d2T_off;
sacMod.opt_d2T_gain = opt_d2T_gain;
sacMod.poss_d2T_off = poss_d2T_off;
sacMod.poss_d2T_gain = poss_d2T_gain;
sacMod.lambda_xvLLImp = (xvLLs - null_xvLL)/log(2);
