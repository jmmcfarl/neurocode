function [allSacMods] = sacMod_scan_doubleregularization_noXV(base_mod,Robs,Xsac,gain_sigs,fit_inds,basemod_pred_rate,...
    poss_d2T_off,poss_d2T_gain,poss_L2_gain,off_only)
%[allSacMods,pred_rate] = sacMod_scan_regularization(base_mod,Robs,Xsac,gain_sigs,poss_d2T,poss_L2)
%Fits a (post-filtering) saccade gain/offset model. Can use a number of
%independent gain filters specified by the number of gain signals in
%"gain_sigs". Also uses xvalLL to scan a range of d2T and L2 reg
%parameters.

if nargin < 8 || isempty(poss_d2T_gain)
    poss_d2T_gain = 0;
end
if nargin < 9 || isempty(poss_L2_gain)
    poss_L2_gain = 0;
end
if nargin < 10 
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
sac_reg_params = NMMcreate_reg_params('boundary_conds',repmat([Inf 0 0],length(mod_signs),1));

%%

    for jj = 1:length(poss_d2T_off)
        for ii = 1:length(poss_d2T_gain)
            for kk = 1:length(poss_L2_gain)
%                 [jj ii kk]
                cur_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
                cur_mod.mods(1).filtK(:) = 1; %initialize base gains to 1
                cur_mod.spk_NL_params = base_mod.spk_NL_params;
                cur_mod = NMMadjust_regularization(cur_mod,[2],'lambda_d2T',poss_d2T_off(jj));
                if ~off_only
                    cur_mod = NMMadjust_regularization(cur_mod,[3],'lambda_d2T',poss_d2T_gain(ii),'lambda_L2',poss_L2_gain(kk));
                    cur_mod = NMMfit_filters(cur_mod,Robs,tr_stim,[],fit_inds,1,[],[],[2 3]);
                else
                    cur_mod = NMMfit_filters(cur_mod,Robs,tr_stim,[],fit_inds,1);
                end
                cur_mod = NMMfit_logexp_spkNL(cur_mod,Robs,tr_stim,[],fit_inds);
                
                [LL,~,pred_rate,~,~,~,nullLL] = NMMmodel_eval(cur_mod,Robs(fit_inds),get_Xcell_tInds(tr_stim,fit_inds));
                cur_mod.ovInfo = mean(pred_rate/mean(pred_rate).*log2(pred_rate/mean(pred_rate)));
                cur_mod.ovLLimp = (LL-nullLL)/log(2);
                cur_mod.nullLL = nullLL;
                
                [sac_offset,sac_gain,sac_info,sac_LLimp] = deal(nan(n_slags,1));
                for ss = 1:n_slags
                    temp = find(Xsac(fit_inds,ss) == 1);
                    cur_Robs = Robs(fit_inds(temp));
                    
                    %this code uses just the models for calculations. It
                    %takes longer but factors out any autocorrelation in
                    %saccade timing
%                     poss_Xsac = zeros(length(fit_inds),size(Xsac,2));
%                     poss_Xsac(:,ss) = 1;
%                     cur_tr_stim{1} = gain_sigs(fit_inds,:);
%                     cur_tr_stim{2} = poss_Xsac; %saccade timing indicator matrix
%                     if ~off_only
%                         cur_tr_stim{3} = reshape(bsxfun(@times,poss_Xsac,reshape(gain_sigs(fit_inds,:),[],1,n_gains)), length(fit_inds),[]);
%                     end
%                     [~,~,pred_rate] = NMMmodel_eval(cur_mod,[],cur_tr_stim);
%                      rr = regress(pred_rate,[ones(length(pred_rate),1) basemod_pred_rate]);
%                     sac_info(ss) = nanmean(pred_rate.*log2(pred_rate/mean(pred_rate)))/mean(pred_rate);
 
                    rr = regress(pred_rate(temp),[ones(length(temp),1) basemod_pred_rate(temp)]);
                   sac_offset(ss) = rr(1);
                    sac_gain(ss) = rr(2);
%                     rr = polyfit(basemod_pred_rate(temp),pred_rate(temp),3);
%                    sac_offset(ss) = rr(4);
%                     sac_gain(ss) = rr(3);
%                     sac_2(ss) = rr(2);
%                     sac_3(ss) = rr(1);
                    
                    sac_info(ss) = nanmean(pred_rate(temp).*log2(pred_rate(temp)/mean(pred_rate(temp))))/mean(pred_rate(temp));
                    
                    %compute nullLL for data at this latency
                    cur_nullLL = nansum(cur_Robs.*log2(mean(cur_Robs)) - mean(cur_Robs));
                    cur_Nspks = sum(cur_Robs);
                    
                    sac_LLimp(ss) = (nansum(cur_Robs.*log2(pred_rate(temp)) - pred_rate(temp)) - cur_nullLL)/cur_Nspks;
                end
                
                cur_mod.sac_offset = sac_offset;
                cur_mod.sac_gain = sac_gain;
                cur_mod.sac_modinfo = sac_info;
                cur_mod.sac_LLimp = sac_LLimp;
                
                allSacMods{jj,ii,kk} = cur_mod;
            end
        end
    end


