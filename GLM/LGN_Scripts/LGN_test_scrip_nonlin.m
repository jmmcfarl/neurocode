%%
clear all
load /Users/James/Data/LGN/LGN_FFdata.mat

% loaded variables:
% DTstim: stimulus temporal dt (s)
% FFspks: spike times (in seconds)
% FFspksR: spike time data with repeats
% FFstim: stimulus time series
% FFstimR: stim time series for repeat

%%
%%% SCRIPT FOR LGN temporal fits

Tp0.stimres = DTstim*1000;
Tp0.frac = 8;
Tp0.NK = 25;
Tp0.dK = 4;
Tp0.latK = 0;
Tp0.NW = 20;
Tp0.dW = (0:2:40);
Tp0.NRP = 0;
Tp0.dRP = [];

% also
%Tp0 = paramsFF( 8, 4, 1 );

fit0 = RGfit_h( FFstim, FFspks, [], [], [], [10 0], Tp0 );
fit0.lambdas(1) = 100;
fit0 = RGfit_h( FFstim, FFspks, fit0 );

TpH = Tp0;
TpH.dRP = [0 1 2 3 4 5 6 7 8 10 12 16 20 24 30 40 50 60];
TpH.NRP = 17;
fitH = RGfit_h( FFstim, FFspks, fit0, [], [], [], TpH );
fp = RGfitprocess( fitH);

%% test reverse correlation technique for effective filters
noise_stim = randn(1e6,1);
NT = length(noise_stim)*1;

lin_out = conv(noise_stim,fitH.K,'same');
temp = xcov(lin_out,noise_stim,TpH.NK);
temp(1:ceil(TpH.NK/2)) = []; temp(end-ceil(TpH.NK/2)+1:end) = [];
figure
plot(temp/norm(temp))
hold on
plot(fitH.K/norm(fitH.K),'r')

%% Add nonlinear "suppressive" "module"
WmI = NLmod_createFF( fitH, -1, FFstim );    % Note: this is old function with lots of crap

fitH.lambdas = [100 100 1 1];

% Add module and fit PSC term 
fit1 = RGfit_h( FFstim, FFspks, fitH, 1, WmI );
% Fit nonlinearity
fit1 = RGfit_nl( FFstim, FFspks, fit1, 1 );
% Fit both (alternating)
fit1 = RGfit_alt( FFstim, FFspks, fit1, 1 );

% % "Old method: refine RF"
% fit1r = RGrefine( FFstim, FFspks, fit1, 1 );

% New method to fit RF -- note this has not been modified to use log(1+e) nonlinearity 
% -- so LLs will be slightly inconsistent (this function uses exp spiking nonlin)
fit1 = RGfit_RFalt( FFstim, FFspks, fit1, 1 );

%% test reverse correlation technique for effective filters (NL module)
noise_stim = randn(1e6,1);
NT = length(noise_stim)*1;

Nkbc = length(fit1.Wmod.KBc);
kb = keat_basis(fit1.Wmod.tf,TpH.frac,Nkbc,fit1.Wmod.stimres);
kC = filter_recompose( fit1.Wmod.KBc, kb );

kern_out = conv(noise_stim,kC,'same');
NL_kern_out = nlin_proc_stim(kern_out,fit1.Wmod.NL,fit1.Wmod.NLx);
WNL_kern_out = conv(NL_kern_out,fit1.Wmod.kw,'same');

temp = xcov(WNL_kern_out,noise_stim,length(kC));
temp(1:ceil(length(kC)/2)) = []; temp(end-ceil(length(kC)/2)+1:end) = [];

% temp2 = xcov(kern_out,noise_stim,length(kC));
% temp2(1:ceil(length(kC)/2)) = []; temp2(end-ceil(length(kC)/2)+1:end) = [];

figure
plot(-temp/norm(temp))
hold on
plot(temp2/norm(temp2),'r')

k = Wmod_filter(fit1.Wmod); %recomposed filter kernel
w = wpiece( fit1.Wmod.kw, fit1.Wmod.dW ); %PSC kernel recomposed
k = g_convolve(k,w,1); %convolve filter kernel and PSC kernel


%%
%%% USING KEAT BASIS -- most modern
TpX = TpH;
TpX.tf = 16*DTstim*1000;
TpX.NK = 12;
TpX.KBshift = 0;

Xstim = KBXstim( FFstim, TpX ); %convolve stim with each KB function
fitHkb = KBfit_h( Xstim, FFspks, [], [], [], [100 100 0 10], TpX );

WmI = RNLmod_create( fitHkb.KBc, TpX, -1, Xstim );    % Note: this is old function with lots of crap
fit1kb = KBfit_h( Xstim, FFspks, fitHkb, 1, WmI );
fit1kb = KBfit_nl( Xstim, FFspks, fit1kb, 1 );
fit1kb.lambdas(2) = 100;
fit1kb = KBfit_alt( Xstim, FFspks, fit1kb, 1 );

% Delay inhibitory filter
fit1kb.Wmod(1) = KBshift_filter( fit1kb.Wmod(1), -6 );
fit1kb = KBfit_h( Xstim, FFspks, fit1kb, 1 );
fit1kb = KBfit_alt( Xstim, FFspks, fit1kb, 1 );

fit1kb = KBfit_RF( Xstim, FFspks, fit1kb, 1 );
% or alternating...
fit1kb = KBfit_RFalt( Xstim, FFspks, fit1kb, 1 );



plot((1:length(fit0.K))*4, fit0.K/max(fit0.K),'r')
plot(fit0kb.K/max(fit0kb.K),'b')