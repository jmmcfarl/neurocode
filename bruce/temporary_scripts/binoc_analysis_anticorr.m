close all
clear all
clc
load('/home/james/Downloads/cell8.mat')

%%
% flen = params_stim.dims(1);
% nPix = params_stim.dims(2)/2;
% 
% use_pix = 20:50;
% [Xinds,Finds] = meshgrid(1:nPix*2,1:flen); %position and time-lag indices associated with each filter element
% raw_stim = Xcell{1}(:,Finds == 1 & ismember(mod(Xinds,nPix),use_pix));
% clear Xcell
% 
% nPix = length(use_pix);
% 
% TB_dsfac = 1;
% flen = round(10/TB_dsfac);
% mod_stim_params = NIM.create_stim_params([flen nPix*2],'tent_spacing',TB_dsfac,...
%     'boundary_conds',[0 0 0],'split_pts',[2 nPix 0]);
% TB_stim = NIM.create_time_embedding(raw_stim,mod_stim_params);
% Xcell{1} = TB_stim;
% 
% anti_set = find(corrs == -1);
% corr_set = find(corrs == 1);
% uncorr_set = find(abs(corrs) ~= 1);
% 
% neural_delay = 3;
% % fit_set = sort([anti_set corr_set] + neural_delay);
% % fit_set = sort([corr_set] + neural_delay);
% % fit_set(fit_set > length(corrs)) = [];
% fit_set = 1:length(corrs);
% corfit_set = sort([corr_set] + neural_delay);
% corfit_set(corfit_set > length(corrs)) = [];
%%
flen = params_stim.dims(1);
nPix = params_stim.dims(2)/2;

use_pix = 20:50;
[Xinds,Finds] = meshgrid(1:nPix*2,1:flen); %position and time-lag indices associated with each filter element
use_lag = 4;
Xcell{1} = Xcell{1}(:,Finds == 4 & ismember(mod(Xinds,nPix),use_pix));
nPix = length(use_pix);

flen = 1;
mod_stim_params = NIM.create_stim_params([nPix*2 flen],...
    'boundary_conds',[0 0 0],'split_pts',[1 nPix 0]);

anti_set = find(corrs == -1);
corr_set = find(corrs == 1);
uncorr_set = find(abs(corrs) ~= 1);

% fit_set = 1:length(corrs);
neural_delay = use_lag;
% fit_set = sort([anti_set corr_set] + neural_delay);
fit_set = sort([corr_set] + neural_delay);
fit_set(fit_set > length(corrs)) = [];

%%
%baseline regularization strengths
lambda_d2XT = 0;
lambda_d2X = 0;
lambda_d2T = 200;
lambda_L1 = 0;

%stimulus subunits
n_Efilts = 4; %number of excitatory squared filters
n_Ifilts = 4; %number of inhibitory squared filters
mod_signs = [ones(1,n_Efilts) -1*ones(1,n_Ifilts)];
NL_types = repmat({'rectpow'},1,n_Efilts + n_Ifilts); %NL types for full model
NL_params = 1.5; %exponent

optim_params.optTol = 1e-4;
optim_params.progTol = 1e-10;
optim_params.maxIter = 2e3;
%% estimate initial models and adjust base reg strength
% init_filts = cell(length(mod_signs),1);
% for ii = 1:length(mod_signs)
%     init_filts{ii} = randn(flen*nPix*2,1)/1e5;
% end
cur_mod = NIM(mod_stim_params,NL_types,mod_signs,'nlparams',NL_params,'spknl','softplus',...
    'd2xt',lambda_d2XT,'d2x',lambda_d2X','d2t',lambda_d2T,'l1',lambda_L1);
% for ii = 1:length(mod_signs)
%    if mod_signs(ii) == 1
% %    if ii <= 3
%        cur_mod.subunits(ii).NLparams = 2;
%    else
%        cur_mod.subunits(ii).NLparams = 1;
%    end
% end
cur_mod = cur_mod.fit_filters(Robs,Xcell,fit_set,'fit_offsets',1,'optim_params',optim_params);
% [~,~,mod_internals] = cur_mod.eval_model(Robs,Xcell,fit_set);
% new_lambda_d2xt = lambda_d2XT./var(mod_internals.gint);
% new_lambda_l1 = lambda_L1./std(mod_internals.gint);
% cur_mod = cur_mod.set_reg_params('d2xt',new_lambda_d2xt,'l1',new_lambda_l1);
% cur_mod = cur_mod.fit_filters(Robs,Xcell,fit_set,'fit_offsets',1,'optim_params',optim_params);
% 
%%
cur_mod2 = cur_mod;
for ii = 1:length(mod_signs)
    if mod_signs(ii) == -1
% if ii > 3
    cur_mod2.subunits(ii).NLparams = 0.75;
end
end
cur_mod2 = cur_mod2.fit_filters(Robs,Xcell,fit_set,'fit_offsets',1,'optim_params',optim_params);

%%
cur_mod2 = cur_mod.fit_NLparams(Robs,Xcell,fit_set);

for ii = 1:3
cur_mod2 = cur_mod2.fit_filters(Robs,Xcell,fit_set,'fit_offsets',1,'optim_params',optim_params);
cur_mod2 = cur_mod2.fit_NLparams(Robs,Xcell,fit_set);
cur_mod2 = cur_mod2.fit_spkNL(Robs,Xcell,fit_set);
end
% cur_mod2 = cur_mod2.fit_spkNL(Robs,Xcell);

%%
% cur_mod3 = cur_mod.fit_spkNL(Robs,Xcell,fit_set);

[~,pred_rate,mod_internals] = cur_mod2.eval_model(Robs,Xcell);
% [~,pred_rate] = cur_mod3.eval_model(Robs,Xcell);

g_exc = sum(mod_internals.fgint(:,mod_signs == 1),2);
g_inh = sum(mod_internals.fgint(:,mod_signs == -1),2);
g_exc = g_exc - mean(g_exc);
g_inh = g_inh - mean(g_inh);

%%
nbins = 100;
G = mod_internals.G;
G_prc = prctile(G,linspace(0,100,nbins));
[n,bin_ids] = histc(G,G_prc);
np_rate = nan(nbins,1);
for ii = 1:nbins
    np_rate(ii) = mean(Robs(bin_ids == ii));
%     np_rate(ii) = mean(Robs(bin_ids == ii & ismember((1:length(bin_ids))',corfit_set )));
end
pred_rate = nan(length(corrs),1);
pred_rate = np_rate(bin_ids);
%%
% cur_mod_rp = cur_mod;
% cur_mod_rp.spkNL.type = 'rectpow';
% cur_mod_rp.spkNL.params = [1 2];
% cur_mod_rp = cur_mod_rp.fit_spkNL(Robs,Xcell,fit_set);
% 
% cur_mod_rp = cur_mod_rp.fit_NLparams(Robs,Xcell,fit_set);
% cur_mod_rp = cur_mod_rp.fit_filters(Robs,Xcell,fit_set);
% 
% [~,pred_rate] = cur_mod_rp.eval_model(Robs,Xcell);

%%
% % cur_mod = cur_mod.fit_NLparams(Robs,Xcell);
% cur_mod2 = cur_mod.init_nonpar_NLs(Xcell,'lambda_nld2',10);
% cur_mod2 = cur_mod2.fit_upstreamNLs(Robs,Xcell);
% cur_mod2 = cur_mod2.fit_spkNL(Robs,Xcell);
% 
% [~,pred_rate] = cur_mod2.eval_model(Robs,Xcell);

%%

maxlag = 10;
pad_pred_rate = cat(1,pred_rate,nan(maxlag+1,1));
% pad_pred_rate2 = cat(1,pred_rate2,nan(maxlag+1,1));
pad_Robs = cat(1,Robs',nan(maxlag+1,1));
poss_disp = unique(dxs(~isnan(dxs)));
pad_gexc = cat(1,g_exc,nan(maxlag+1,1));
pad_gsup = cat(1,g_inh,nan(maxlag+1,1));

[pred_anti_dtune,pred_corr_dtune,emp_corr_dtune,emp_anti_dtune,pred_anti2,pred_corr2] = deal(nan(length(poss_disp),maxlag));
[anti_gsup,anti_gexc,corr_gsup,corr_gexc] = deal(nan(length(poss_disp),maxlag));
for ii = 1:length(poss_disp)
   cur_set = corr_set(dxs(corr_set) == poss_disp(ii));
   for ff = 1:maxlag
      pred_corr_dtune(ii,ff) = nanmean(pad_pred_rate(cur_set + ff - 1));
%       pred_corr2(ii,ff) = nanmean(pad_pred_rate2(cur_set + ff - 1));
      emp_corr_dtune(ii,ff) = nanmean(pad_Robs(cur_set + ff - 1));
      corr_gexc(ii,ff) = nanmean(pad_gexc(cur_set + ff - 1));
      corr_gsup(ii,ff) = nanmean(pad_gsup(cur_set + ff - 1));
   end
   cur_set = anti_set(dxs(anti_set) == poss_disp(ii));
   for ff = 1:maxlag
      pred_anti_dtune(ii,ff) = nanmean(pad_pred_rate(cur_set + ff - 1));
%       pred_anti2(ii,ff) = nanmean(pad_pred_rate2(cur_set + ff - 1));
      emp_anti_dtune(ii,ff) = nanmean(pad_Robs(cur_set + ff - 1));
      anti_gexc(ii,ff) = nanmean(pad_gexc(cur_set + ff - 1));
      anti_gsup(ii,ff) = nanmean(pad_gsup(cur_set + ff - 1));
   end
end

%%
figure
poss_lags = [2 3 4 5];
for ii = 1:4
    subplot(2,2,ii)
    plot(emp_anti_dtune(:,poss_lags(ii)))
    hold on
    plot(emp_corr_dtune(:,poss_lags(ii)),'r')
    plot(pred_anti_dtune(:,poss_lags(ii)),'--')
    plot(pred_corr_dtune(:,poss_lags(ii)),'r--')
%     plot(pred_anti2(:,poss_lags(ii)),'o-')
%     plot(pred_corr2(:,poss_lags(ii)),'ro-')
    title(sprintf('lag %d',poss_lags(ii)));
end

% figure
% poss_lags = [2 3 4 5];
% for ii = 1:4
%     subplot(2,2,ii)
%     plot(corr_gexc(:,poss_lags(ii)))
%     hold on
%     plot(corr_gsup(:,poss_lags(ii)),'k')
%     plot(anti_gexc(:,poss_lags(ii)),'g')
%     plot(anti_gsup(:,poss_lags(ii)),'r')
%     title(sprintf('lag %d',poss_lags(ii)));
% end