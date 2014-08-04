%%
clear all
load monoc_datadump

flen = 13; %number of lags
dt = 0.01; %in sec
dx = 0.0283; %in deg
spatial_usfac = 2; %spatial up-sampling
use_nPix_us = size(cor_stim,2); %number of (upsampled) pixels
stim_params = NMMcreate_stim_params([flen use_nPix_us]);
Xmat = create_time_embedding(cor_stim,stim_params);

%%
cc = 103; %SUs are cc > 96 (103 is best)

%get relevant spike and stim variables
cur_Robs = Robs_mat(:,cc);
cc_uinds = find(~isnan(cur_Robs));
cur_Robs = cur_Robs(cc_uinds);
cur_Xmat = Xmat(used_inds(cc_uinds),:);

nim_stim_params = NMMcreate_stim_params([flen use_nPix_us]);
init_lambda_d2XT = 20;
nim_reg_params = NMMcreate_reg_params('lambda_d2XT',init_lambda_d2XT);

silent = 0;
mod_signs = [1 1 1];
NL_types = {'lin','quad','quad'};
GQM = NMMinitialize_model(nim_stim_params,mod_signs,NL_types,nim_reg_params);
GQM = NMMfit_filters(GQM,cur_Robs,cur_Xmat,[],[],silent);

optim_params.optTol = 1e-4;
optim_params.progTol = 1e-8;

%adjust regularization and refit
base_lambda_d2XT = 200; %300
base_lambda_L1 = 20; %15
[LL, penLL, pred_rate, G, gint] = NMMmodel_eval(GQM,cur_Robs,cur_Xmat);
GQM = NMMadjust_regularization(GQM,[],'lambda_d2XT',base_lambda_d2XT./var(gint)');
GQM = NMMadjust_regularization(GQM,[],'lambda_L1',base_lambda_L1./std(gint)');
GQM = NMMfit_filters(GQM,cur_Robs,cur_Xmat,[],[],0,optim_params);

NMMdisplay_model(GQM);


