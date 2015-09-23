%% NIM -- LGN temporal demo... Comparing NIM class and NIMtoolbox
% Make sure the NIMtoolbox and NIM class are in your Matlab path
% I will denote old models in lower case and new models (objects) as upper case

%% Listing available methods
%can list all the methods associated with the NIM class
methods(NIM)

%similarly for the SUBUNIT class
methods(SUBUNIT)

%% GETTING HELP WITH SPECIFIC METHODS (fit_filters)
help NIM.fit_filters

%% GETTING HELP WITH SPECIFIC METHODS (NIM constructor)
help NIM.NIM

%% Load and format data
% MAT-file contains: stimulus (FFstim), spike times (FFspks), and time res of stim (DTstim)
% Also includes cross-validation stim (FFstimR) and 64 repeats (FFspksR)
load LGN_FFdata.mat

% Format data for fitting
up_samp_fac = 16; % temporal up-sampling factor applied to stimulus 
tent_basis_spacing = 4; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
nLags = 50; % number of time lags for estimating stimulus filters

%make stim params struct
stim_params = NIMcreate_stim_params(nLags, DTstim, up_samp_fac, tent_basis_spacing); %for nim
STIM_PARAMS = NIM.create_stim_params(nLags,'stim_dt',DTstim,'up_samp_fac',up_samp_fac,'tent_spacing',tent_basis_spacing); %make STIM_PARAMS struct for NIM

% Create T x nLags 'Xmatrix' representing the relevant stimulus history at each time point
Xstim = create_time_embedding( FFstim, stim_params ); 
Xcell{1} = NIM.create_time_embedding(FFstim,STIM_PARAMS); %for the NIM class, need the stimulus formatted as a cell array

% Bin spikes at analysis resolution
Robs = histc( FFspks, (0:(size(Xstim,1)-1))*stim_params.dt ); % Prob w conventions

%% Fit a single-filter LN model
%specify model structure
NL_types = {'lin'}; %define subunit as linear 
mod_signs = [1]; %determines whether input is exc or sup (doesn't matter in the linear case)

% Initialize nim 
params_reg = NIMcreate_reg_params('lambda_d2T',80); % this is purposely set too low here 
fit0 = NIMinitialize_model( stim_params, mod_signs, NL_types, params_reg ); 
FIT0 = NIM(STIM_PARAMS,NL_types,mod_signs,'d2t',80); %first make model

%the default optimization parameters are a bit different in the NIM class, so setting them to old
%values for more direct comparison.
new_optim_params.optTol = 1e-4; 
new_optim_params.progTol = 1e-8;

% Fit stimulus filter
silent = 1; %don't print iterative optimization output
fit0 = NIMfit_filters(fit0, Robs, Xstim, [], [], silent); 
FIT0 = FIT0.fit_filters(Robs,Xcell,'optim_params',new_optim_params,'silent',silent);

filtfig = figure(); hold on
plot(fit0.mods(1).filtK,'r') %nim filter
plot(FIT0.subunits(1).filtK,'k--') %NIM filter

%% Add spike-history term and refit
n_bins = 20; %number of bins to use for representing spike-hist filter
fitH = NIMinit_spkhist( fit0, n_bins, 1, 5 );
fitH = NIMfit_filters( fitH, Robs, Xstim, [], [], silent ); 

FITH = FIT0.init_spkhist(n_bins,'init_spacing',1,'doubling_time',5);
FITH = FITH.fit_filters(Robs,Xcell,'optim_params',new_optim_params,'silent',silent);

%% Make 2-filter model, linear plus suppressive; like Butts et al (2011)
% Add an inhibitory input and fit 

% Make delayed suppression an initial guess (note that not delaying will
% result in a different local minima
ksup = fit0.mods(1).filtK;
ksup(4:end) = ksup(1:(end-3));
fit1 = NIMadd_NLinput( fitH, 'threshlin', -1, ksup );
FIT1 = FITH.add_subunits('rectlin',-1,'init_filts',ksup);

%fit filters
fit1 = NIMfit_filters( fit1, Robs, Xstim, [], [], silent );
FIT1 = FIT1.fit_filters(Robs,Xcell,'optim_params',new_optim_params,'silent',silent);

%display results
NIMdisplay_model( fit1, Xstim )
FIT1.display_model(Robs,Xcell);

%% Convert upstream NLs to nonparametric (monotonic) functions and fit
to_nonpar = [2]; % which subunits to convert to tent-basis (nonparametric) upstream NLs (just nonlinear term)
lambda_NLd2 = 20; % smoothness regularization on the tent-basis coefs
fit2 = NIMinitialize_upstreamNLs( fit1, Xstim, to_nonpar, lambda_NLd2 ); % initializes the tent-basis representation
FIT2 = FIT1.init_nonpar_NLs(Xcell,'sub_inds',to_nonpar,'lambda_nld2',lambda_NLd2);

%these are the default tolerances used in NIMfit_upstreamNLs
NL_optim_params.TolX = 1e-6;
NL_optim_params.TolFun = 1e-6;

%now fit upsteram NLs
fit2 = NIMfit_upstreamNLs( fit2, Robs, Xstim, [],[],[], silent ); % fit the upstream NL
FIT2 = FIT2.fit_upstreamNLs(Robs,Xcell,'optim_params',NL_optim_params,'silent',silent);

%% Do another iteration of fitting filters and upstream NLs
fit2 = NIMfit_filters( fit2, Robs, Xstim, [],[], silent );
FIT2 = FIT2.fit_filters(Robs,Xcell,'silent',silent,'optim_params',new_optim_params);

fit2 = NIMfit_upstreamNLs( fit2, Robs, Xstim, [],[],[], silent );
FIT2 = FIT2.fit_upstreamNLs(Robs,Xcell,'silent',silent,'optim_params',NL_optim_params);

%plot LL sequences
f2 = figure(); hold on
plot(fit2.LL_seq,'o-');
plot([FIT2.fit_hist(:).LL],'ro-');

%% Fit spiking NL params
% Note that in this case re-estimating the spiking NL shape doesnt improve things 
% although this is not always the case.
fit3 = NIMfit_logexp_spkNL( fit2, Robs, Xstim,[],0);
FIT3 = FIT2.fit_spkNL(Robs, Xcell,'silent',silent);

