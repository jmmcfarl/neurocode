function [sacMod] = fit_sacMod_basic(base_mod,Robs,Xsac,gain_sigs,fit_inds,d2T_off,d2T_gain,off_only)
%[allSacMods,pred_rate] = sacMod_scan_regularization(base_mod,Robs,Xsac,gain_sigs,poss_d2T,poss_L2)
%Fits a (post-filtering) saccade gain/offset model. Can use a number of
%independent gain filters specified by the number of gain signals in
%"gain_sigs". Also uses xvalLL to scan a range of d2T and L2 reg
%parameters.

if nargin < 8 
    off_only = false;
end

n_slags = size(Xsac,2);
n_gains = size(gain_sigs,2);

sac_stim_params(1) = NMMcreate_stim_params(n_gains);
sac_stim_params(2) = NMMcreate_stim_params(n_slags);

tr_stim{1} = gain_sigs;
tr_stim{2} = Xsac; %saccade timing indicator matrix
if ~off_only
    tr_stim{3} = reshape(bsxfun(@times,Xsac,reshape(gain_sigs,[],1,n_gains)), size(gain_sigs,1),[]);
    sac_stim_params(3) = NMMcreate_stim_params([n_slags n_gains]);
end

if ~off_only
    mod_signs = [1 1 1];
    Xtargets = [1 2 3];
    NL_types = {'lin','lin','lin'};
else
    mod_signs = [1 1];
    Xtargets = [1 2];
    NL_types = {'lin','lin'};
end
silent = 1;
sac_reg_params = NMMcreate_reg_params('boundary_conds',repmat([Inf 0 0],length(mod_signs),1));

%%

cur_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
cur_mod.mods(1).filtK(:) = 1; %initialize base gains to 1
cur_mod.spk_NL_params = base_mod.spk_NL_params;
cur_mod = NMMadjust_regularization(cur_mod,2,'lambda_d2T',d2T_off); %temporal smoothness reg for offset filter
if ~off_only
    cur_mod = NMMadjust_regularization(cur_mod,3,'lambda_d2T',d2T_gain); %temporal smoothness reg for gain filter
    cur_mod = NMMfit_filters(cur_mod,Robs,tr_stim,[],fit_inds,silent,[],[],[2:3]);
else
    cur_mod = NMMfit_filters(cur_mod,Robs,tr_stim,[],fit_inds,silent,[],[],[2]);
end
sacMod = NMMfit_logexp_spkNL(cur_mod,Robs,tr_stim,[],fit_inds);

[LL,nullLL,pred_rate] = NMMeval_model(sacMod,Robs,tr_stim,[],fit_inds);
sacMod.ovInfo = mean(pred_rate/mean(pred_rate).*log2(pred_rate/mean(pred_rate)));
sacMod.ovLLimp = (LL-nullLL)/log(2);
sacMod.nullLL = nullLL;


