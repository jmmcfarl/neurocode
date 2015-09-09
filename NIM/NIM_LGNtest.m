%% NIM -- LGN temporal demo
% Make sure the NIMtoolbox is in your Matlab path
clear all

%% Load and format data
% MAT-file contains: stimulus (FFstim), spike times (FFspks), and time res of stim (DTstim)
% Also includes cross-validation stim (FFstimR) and 64 repeats (FFspksR)
load LGN_FFdata.mat

% Format data for fitting
up_samp_fac = 16; % temporal up-sampling factor applied to stimulus 
tent_basis_spacing = 4; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
nLags = 50; % number of time lags for estimating stimulus filters
params_stim = NIMcreate_stim_params( nLags, DTstim, up_samp_fac, tent_basis_spacing );

% Create T x nLags 'Xmatrix' representing the relevant stimulus history at each time point
Xstim = create_time_embedding( FFstim, params_stim ); 
Xcell{1} = Xstim;

% Bin spikes at analysis resolution
Robs = histc( FFspks, (0:(size(Xstim,1)-1))*params_stim.dt ); % Prob w conventions

%% Fit a single-filter LN model
NL_types = {'lin'}; % define subunit as linear 
mod_signs = [1]; % determines whether input is exc or sup (doesn't matter in the linear case)

stim_params.dims = [nLags]; 
stim_params.dt = DTstim/up_samp_fac*tent_basis_spacing; %time res of filter coeffs

% Initialize NIM (use 'help NIMinitialize_model' for more details about the 
% model struct components).
% params_reg = NIMcreate_reg_params( 'lambda_d2T', 1 ); % this is purposely set too low here 
% fit0 = NIMinitialize_model( params_stim, mod_signs, NL_types, params_reg ); 

new_fit = NIM(stim_params,NL_types,mod_signs,'spkNL','softplus');
new_fit = new_fit.set_reg_params('d2t',1);
new_optim_params.optTol = 1e-4; 
new_optim_params.progTol = 1e-8;

% Fit stimulus filter
tic
silent = 0;
% fit0 = NIMfit_filters( fit0, Robs, Xstim, [], [], silent ); 
toc
tic
% new_fit = new_fit.fit_filters(Robs,Xcell,'optim_params',new_optim_params);
toc
% Plot filter -- should be less smooth than what we would like (our 'prior')
% filtfig = figure; hold on
% plot(fit0.mods(1).filtK,'r')
% plot(new_fit.subunits(1).filtK,'k--')
% 
% Increase smoothness regularization and refit
%   you dould do this: fit0.mods(1).reg_params.lambda_d2T = 80;
%   but we made a function:
% fit0 = NIMadjust_regularization( fit0, 1, 'lambda_d2T', 80);
% fit0 = NIMfit_filters( fit0, Robs, Xstim, [], [], silent ); 

new_fit = new_fit.set_reg_params('d2t',80);
new_fit = new_fit.fit_filters(Robs,Xcell,'optim_params',new_optim_params);

% figure(filtfig); hold on
% plot(fit0.mods(1).filtK,'r')
% plot(new_fit.subunits(1).filtK,'k--')

%%
gain_funs = ones(size(Robs'))*2;
new_fit2 = new_fit.fit_filters(Robs,Xcell,'optim_params',new_optim_params,'gain_funs',gain_funs);
%%
% fitS = NIMfit_logexp_spkNL( fit0, Robs, Xstim );
% new_fit.spkNL.theta = fitS.spk_NL_params(1);
% new_fit.spkNL.params = [fitS.spk_NL_params(2:3)];
new_fitS = new_fit.fit_spkNL(Robs,Xcell);

%% GLM: add spike-history term and refit
n_bins = 20;
fitH = NIMinit_spkhist( fit0, n_bins, 1, 5 );

new_fitH = new_fit.init_spkhist(n_bins,'init_spacing',1,'doubling_time',5);
new_fitH.subunits(1).filtK = fitH.mods(1).filtK;

fitH = NIMfit_filters( fitH, Robs, Xstim, [], [], silent ); 

% new_fitH.spk_hist.coefs = fitH.spk_hist.coefs;
% new_fitH.spkNL.theta = fitH.spk_NL_params(1);
% new_fitH.subunits(1).filtK = fitH.mods(1).filtK;
new_fitH = new_fitH.fit_filters(Robs,Xcell,'optim_params',new_optim_params);

%% NIM: linear plus suppressive; like Butts et al (2011)
% Add an inhibitory input and fit 
%fit1 = NIMadd_module( fitH, {'threshlin'}, -1, 1, 4 );
% Make delayed suppression an initial guess (note that not delaying will
% result in a different local minima

ksup = fit0.mods(1).filtK;
ksup(4:end) = ksup(1:(end-3));
fit1 = NIMadd_NLinput( fitH, 'threshlin', -1, ksup );

new_fit1 = new_fitH.add_subunits('rectlin',-1,'init_filts',ksup);

% Fit 
fit1 = NIMfit_filters( fit1, Robs, Xstim, [], [], silent );

new_fit1 = new_fit1.fit_filters(Robs,Xcell,'optim_params',new_optim_params);

% Compare two kernels
figure(filtfig); hold on; 
plot(fit1.mods(1).filtK,'b')
plot(fit1.mods(2).filtK,'r')
plot(new_fit1.subunits(1).filtK,'k--')
plot(new_fit1.subunits(2).filtK,'g--')

% Add 'Xstim' to see how nonlinearities compare with stim distributions
NIMdisplay_model( fit1, Xstim )

new_fit1.display_model(Xcell);

%% Convert upstream NLs to nonparametric (monotonic) functions and fit
to_nonpar = [2]; % which subunits to convert to tent-basis (nonparametric) upstream NLs (just nonlinear term)
lambda_NLd2 = 20; % smoothness regularization on the tent-basis coefs
fit2 = NIMinitialize_upstreamNLs( fit1, Xstim, to_nonpar, lambda_NLd2 ); % initializes the tent-basis representation
new_fit2 = new_fit1.init_nonpar_NLs(Xcell,'sub_inds',to_nonpar,'lambda_nld2',lambda_NLd2);

% Now fit the NLs (+ spk history term) -- where is spike history term?
fit2 = NIMfit_upstreamNLs( fit2, Robs, Xstim, [],[],[], silent ); % fit the upstream NL
% new_fit2 = new_fit2.fi
% new_fit2.subunits(2).TBx = fit2.mods(2).NLx;
% new_fit2.subunits(2).TBy = fit2.mods(2).NLy;
% new_fit2.subunits(2).filtK = fit2.mods(2).filtK;
% new_fit2.spk_hist = fit2.spk_hist;
% new_fit2.spkNL.theta = fit2.spk_NL_params(1);
% new_fit2.subunits(1).filtK = fit1.mods(1).filtK;
new_fit2 = new_fit2.fit_upstreamNLs(Robs,Xcell,'optim_params',new_optim_params);

%% Do another iteration of fitting filters and upstream NLs
silent = 1; % switching to stealth mode
fit2 = NIMfit_filters( fit2, Robs, Xstim, [],[], silent );
fit2 = NIMfit_upstreamNLs( fit2, Robs, Xstim, [],[],[], silent );

new_fit2 = new_fit2.fit_filters(Robs,Xcell,'silent');
new_fit2 = new_fit2.fit_upstreamNLs(Robs,Xcell,'silent');

% Note that LLs aren't improving much with these additional optimizations
fprintf('Log-likelihoods so far: \n');
disp(fit2.LL_seq);

%% Fit spiking NL params
% Note that in this case re-estimating the spiking NL shape doesnt improve things 
% although this is not always the case.
spk_NL_display = 1; % turn on display to get a before-and-after picture of the spk NL
fit3 = NIMfit_logexp_spkNL( fit2, Robs, Xstim, [], spk_NL_display );

%% Fit quadratic model (with linear and squared terms)
fitQ = NIMadd_NLinput( fitH, 'quad', -1, ksup );
fitQ = NIMfit_filters( fitQ, Robs, Xstim, [],[], silent );
fitQ = NIMfit_logexp_spkNL( fitQ, Robs, Xstim, [], spk_NL_display );

%% View different models with NIMdisplay_model
NIMdisplay_model(fitQ,Xstim,[],Robs)
NIMdisplay_model(fit3,Xstim,[],Robs)


%%%%%%%%%% CROSS-VALIDATION / MODEL PERFORMANCE %%%%%%%%%%%%%
% Fit spiking nonlinearities on LN and GLM too
fit0 = NIMfit_logexp_spkNL( fit0, Robs, Xstim, [], spk_NL_display );
fitH = NIMfit_logexp_spkNL( fitH, Robs, Xstim, [], spk_NL_display );

% First, check likelihood comparisons on same data relative to null model
avfiringrate = length(FFspks)/(length(FFstim)*DTstim);
LLnull = log(avfiringrate*params_stim.dt)-1;  % formula for LL/spk of null model (in units of nats)

LLsame = [fit0.LL_seq(end) fitH.LL_seq(end) fit3.LL_seq(end) fitQ.LL_seq(end)] - LLnull

%% Next, calculate on cross-validated data (repeats)
XstimR = create_time_embedding( FFstimR, params_stim );

% Histogram repeat data into RobsR matrix (NT x Nreps)
%   FFspksR is a list of spike times with '-1' in between each repeat
Nreps = length(find(FFspksR < 0));
RobsR = zeros(size(XstimR,1),Nreps);

rep_locs = [0 find(FFspksR < 0)];
for rep = 1:Nreps
	rspks = FFspksR( (rep_locs(rep)+1):(rep_locs(rep+1)-1) );
	onetrialhist = histc(rspks, (0:(size(XstimR,1)-1))*params_stim.dt );
	RobsR(:,rep) = onetrialhist';
end

% Calculate likelihoods: one for each repeat
LLs0 = NIMeval_LLreps( fit0, XstimR, RobsR );
LLsH = NIMeval_LLreps( fitH, XstimR, RobsR );
LLs3 = NIMeval_LLreps( fit3, XstimR, RobsR );
LLsQ = NIMeval_LLreps( fitQ, XstimR, RobsR );

[median(LLs0) median(LLsH) median(LLs3) median(LLsQ)]/log(2) % converts to bits/spk

%% Finally, look at PSTH comparisons to calculate R-squared
% Simulate and compare predicted firing rates
dtbin = DTstim/8;  % roughly 1 ms resolution (double analysis res)
T = length(FFstimR)*DTstim;

% Experiment PSTH
Nreps = sum(FFspksR < 0)
psth = histc( FFspksR, 0:dtbin:T )' / Nreps / dtbin;

% Simulate different models (1000 reps to average over spike history)
Nrepsim = 1000;

[r0 spks0] = NIMsimulate( fit0, XstimR, [], Nrepsim );
pred0 = histc( spks0, 0:dtbin:T ) / Nrepsim / dtbin;

% Note: these take a while due to spike history term
[rH spksH] = NIMsimulate( fitH, XstimR, [], Nrepsim );
[r3 spks3] = NIMsimulate( fit3, XstimR, [], Nrepsim );
[rQ spksQ] = NIMsimulate( fitQ, XstimR, [], Nrepsim );

predH = histc( spksH, 0:dtbin:T ) / Nrepsim / dtbin;
pred3 = histc( spks3, 0:dtbin:T ) / Nrepsim / dtbin;
predQ = histc( spksQ, 0:dtbin:T ) / Nrepsim / dtbin;

% Figure
twindow = [1.1 1.3]; % sec
tbins = floor(twindow(1)/dtbin):floor(twindow(2)/dtbin);
figure; hold on
plot(tbins*dtbin,psth(tbins),'k')
plot(tbins*dtbin,pred0(tbins),'r')
plot(tbins*dtbin,predH(tbins),'g')
plot(tbins*dtbin,pred3(tbins),'b')
plot(tbins*dtbin,predQ(tbins),'m')
axis tight

% Check traditional measures of R-squared
R2s = [
1- var(psth-pred0)/var(psth)
1- var(psth-predH)/var(psth)
1- var(psth-pred3)/var(psth)
1- var(psth-predQ)/var(psth)]
