clear all;
close all;

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
defmod = pars.defmod;
defmod.h(1:end-1) = [];%restrict PSC term to delta function
flen = pars.flen;
hilen = length(defmod.h);
foff = flen + pars.hilen;

datdir = '~/Data/rust/stcbar/DATA/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);

cd ~/James_scripts/stc_sparse_test_figs/stac_allcells
load ./used_stcims.mat


ncells = length(uids);
nn = 5; flen = 14;


fprintf('ANALYZING CELL %d OF %d\n\n',nn,ncells);

%% load data
eval(['load ',['~/Data/rust/stcbar/DATA/',fnames{nn}]]);
Robs = spikes_per_frm(:);

%%
flen = 12;
[NT,nPix] = size(stim);
nmm_stim_params = NMMcreate_stim_params([flen nPix]);
Xmat = create_time_embedding(stim,nmm_stim_params);

tr_set = 1:round(NT/2);
%% INITIAL MODEL FIT
base_lambda_d2XT = 20;
base_lambda_L1 = 20;
fit0 = NMMinitialize_model(nmm_stim_params,[1 1],{'lin','quad'});
fit0 = NMMfit_filters(fit0,Robs(tr_set),Xmat(tr_set,:),[],[],0);

%SET REG PARAMS
[LL, penLL, pred_rate, G, gint] = NMMmodel_eval(fit0,Robs(tr_set),Xmat(tr_set,:));
Xtargs = [fit0.mods(:).Xtarget];
fit0 = NMMadjust_regularization(fit0,find(Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
fit0 = NMMadjust_regularization(fit0,find(Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
fit0 = NMMfit_filters(fit0,Robs(tr_set),Xmat(tr_set,:),[],[],0);
[fit0_LL, fit0_penLL, pred_rate0, G, fgint, nullLL] = NMMmodel_eval( fit0, Robs(tr_set), Xmat(tr_set,:));
%% FIT CONV KERNELS
addpath('~/James_scripts/nmmc/')
fit1 = NMMCinitialize_Xconv( fit0, [2], 8, 25);
fit1 = NMMCfit_Xconv(fit1, Robs(tr_set), Xmat(tr_set,:),[],[],0);
[fit1_LL, fit1_penLL, pred_rate1, G, fgint, nullLL] = NMMCmodel_eval( fit1, Robs(tr_set), Xmat(tr_set,:));
%% REFIT FILTERS
fit2 = NMMCfit_filters( fit1, Robs(tr_set), Xmat(tr_set,:), [], [], 0);

[fit2_LL, fit2_penLL, pred_rate, G, fgint, nullLL] = NMMCmodel_eval( fit2, Robs(tr_set), Xmat(tr_set,:));