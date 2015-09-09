%% NIM -- LGN temporal demo... Comparing NIM class and NIMtoolbox

% Make sure the NIMtoolbox and NIM class are in your Matlab path
clear all

%I will denote old models in lower case and new models (objects) as upper case

%% Load and format data
% MAT-file contains: stimulus (FFstim), spike times (FFspks), and time res of stim (DTstim)
% Also includes cross-validation stim (FFstimR) and 64 repeats (FFspksR)
load LGN_FFdata.mat

% Format data for fitting
up_samp_fac = 16; % temporal up-sampling factor applied to stimulus 
tent_basis_spacing = 4; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
nLags = 50; % number of time lags for estimating stimulus filters
params_stim = NIMcreate_stim_params( nLags, DTstim, up_samp_fac, tent_basis_spacing ); %for nim

%make STIM_PARAMS struct for NIM
STIM_PARAMS = NIM.create_stim_params(nLags,'stim_dt',DTstim,'up_samp_fac',up_samp_fac,'tent_spacing',tent_basis_spacing);

NL_types = {'lin'}; % define subunit as linear 
mod_signs = [1]; % determines whether input is exc or sup (doesn't matter in the linear case)

% Create T x nLags 'Xmatrix' representing the relevant stimulus history at each time point
Xstim = create_time_embedding( FFstim, params_stim ); 
Xcell{1} = NIM.create_time_embedding(FFstim,STIM_PARAMS); %for the NIM class, need the stimulus formatted as a cell array

% Bin spikes at analysis resolution
Robs = histc( FFspks, (0:(size(Xstim,1)-1))*params_stim.dt ); % Prob w conventions

%% Fit a single-filter LN model

% Initialize nim 
params_reg = NIMcreate_reg_params( 'lambda_d2T', 80); % this is purposely set too low here 
fit0 = NIMinitialize_model( params_stim, mod_signs, NL_types, params_reg ); 

FIT0 = NIM(STIM_PARAMS,NL_types,mod_signs); %first make model
FIT0 = FIT0.set_reg_params('d2t',80); %set regularization parameters
%the default optimization parameters are a bit different in the NIM class, so setting them to old
%values for more direct comparison.
new_optim_params.optTol = 1e-4; 
new_optim_params.progTol = 1e-8;

% Fit stimulus filter
silent = 0;
fit0 = NIMfit_filters(fit0, Robs, Xstim, [], [], silent ); 
FIT0 = FIT0.fit_filters(Robs,Xcell,'optim_params',new_optim_params);

filtfig = figure(); hold on
plot(fit0.mods(1).filtK,'r') %nim filter
plot(FIT0.subunits(1).filtK,'k--') %NIM filter


%% Add spike-history term and refit
n_bins = 20;
fitH = NIMinit_spkhist( fit0, n_bins, 1, 5 );
fitH = NIMfit_filters( fitH, Robs, Xstim, [], [], silent ); 

FITH = FIT0.init_spkhist(n_bins,'init_spacing',1,'doubling_time',5);
FITH = FITH.fit_filters(Robs,Xcell,'optim_params',new_optim_params);

%% Make 2-filter model, linear plus suppressive; like Butts et al (2011)
% Add an inhibitory input and fit 
% Make delayed suppression an initial guess (note that not delaying will
% result in a different local minima
ksup = fit0.mods(1).filtK;
ksup(4:end) = ksup(1:(end-3));
fit1 = NIMadd_NLinput( fitH, 'threshlin', -1, ksup );
fit1 = NIMfit_filters( fit1, Robs, Xstim, [], [], silent );

ksup = FIT0.subunits(1).filtK;
ksup(4:end) = ksup(1:(end-3));
FIT1 = FITH.add_subunits('rectlin',-1,'init_filts',ksup);
FIT1 = FIT1.fit_filters(Robs,Xcell,'optim_params',new_optim_params);

NIMdisplay_model( fit1, Xstim )
FIT1.display_model(Robs,Xcell);

%% Convert upstream NLs to nonparametric (monotonic) functions and fit
to_nonpar = [2]; % which subunits to convert to tent-basis (nonparametric) upstream NLs (just nonlinear term)
lambda_NLd2 = 20; % smoothness regularization on the tent-basis coefs
fit2 = NIMinitialize_upstreamNLs( fit1, Xstim, to_nonpar, lambda_NLd2 ); % initializes the tent-basis representation
fit2 = NIMfit_upstreamNLs( fit2, Robs, Xstim, [],[],[], silent ); % fit the upstream NL

NL_optim_params.TolX = 1e-6;%these are the default tolerances used in NIMfit_upstreamNLs
NL_optim_params.TolFun = 1e-6;
FIT2 = FIT1.init_nonpar_NLs(Xcell,'sub_inds',to_nonpar,'lambda_nld2',lambda_NLd2);
FIT2 = FIT2.fit_upstreamNLs(Robs,Xcell,'optim_params',NL_optim_params);
%% Do another iteration of fitting filters and upstream NLs
silent = 1; % switching to stealth mode
fit2 = NIMfit_filters( fit2, Robs, Xstim, [],[], silent );
fit2 = NIMfit_upstreamNLs( fit2, Robs, Xstim, [],[],[], silent );

FIT2 = FIT2.fit_filters(Robs,Xcell,'silent','optim_params',new_optim_params);
FIT2 = FIT2.fit_upstreamNLs(Robs,Xcell,'silent','optim_params',NL_optim_params);

% Note that LLs aren't improving much with these additional optimizations
fprintf('nim Log-likelihoods so far: \n');
disp(fit2.LL_seq);
%the exact LLs are a bit different for the NIM after introducing the nonpar upstream NLs. I believe
%this is because the old algorithm wasnt 'rescaling' the NL coefs by default (though it probably
%should have been).
fprintf('NIM Log-likelihoods so far: \n');
disp([FIT2.fit_hist(:).LL]');

%plot LL sequences
f2 = figure(); hold on
plot(fit2.LL_seq,'o-');
plot([FIT2.fit_hist(:).LL],'ro-');

%% Fit spiking NL params
% Note that in this case re-estimating the spiking NL shape doesnt improve things 
% although this is not always the case.
fit3 = NIMfit_logexp_spkNL( fit2, Robs, Xstim);
FIT3 = FIT2.fit_spkNL(Robs, Xcell);

