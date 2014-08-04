%%
clear all
load /Users/James/Data/LGN/LGN_FFdata.mat

% loaded variables:
% DTstim: stimulus temporal dt (s)
% FFspks: spike times (in seconds)
% FFspksR: spike time data with repeats
% FFstim: stimulus time series
% FFstimR: stim time series for repeats

%% SET PARAMETERS
Tp0.stimres = DTstim*1000; %stimulus dt in ms
Tp0.frac = 8; %upsample factor
Tp0.NK = 25; %number of temporal lags for kernel
Tp0.dK = 4; %up-sample factor for temporal kernel ??
Tp0.latK = 0; %??
Tp0.NW = 20; %??
Tp0.dW = (0:20)*Tp0.frac; %??
% Tp0.NRP = 0; %number of history dep terms
% Tp0.dRP = []; %temporal binning for hist dep terms
Tp0.dRP = [0 1 2 3 4 5 6 7 8 10 12 16 20 24 30 40 50 60];
Tp0.NRP = length(Tp0.dRP)-1; %? is this always true

%reg params
lambdas = [100 0 0]; %[kern ? hist]

fitH = RGfit_h( FFstim, FFspks, [], [], [], lambdas, Tp0);
fp = RGfitprocess( fitH);

%% PLOT spike times with overlaid firing rate prediction
% NT = length(FFstim)*Tp0.frac;
% X = zeros(NT,Tp0.NK+Tp0.NRP); %initialize X
% 
% X(:,1:Tp0.NK) = matrix_FFkernel(FFstim, Tp0.frac, Tp0.NK, Tp0.dK, Tp0.latK ); %add the kernel portion of X (series of lagged stim vectors)
% %  Nonmodern version: X(:,1:Tp.NK) = zmatrix_FFkernel( FFstim, Tp );
% 
% %%% Spike history (and spikes-process)
% [Xt spksN] = matrix_spike_history( FFspks, Tp0.dRP, Tp0.stimres/1000/Tp0.frac, NT );
% 
% K(1:Tp0.NK) = fitH.K;
% K(Tp0.NK+(1:length(fitH.wRP))) = fitH.wRP;
% K(end) = fitH.offset;
% 
% k = K(1:end-1);
% b = K(end);
% 
% r = exp(k*X(:,1:end-1)'+b);
% 
% T = Tp0.stimres/1000/Tp0.frac*(1:NT);
% 
% figure
% plot(T,r)
% hold on
% plot(FFspks,ones(size(FFspks)),'k.'),shg

%% validate with repeat data
NT = length(FFstimR)*Tp0.frac;
X = zeros(NT,Tp0.NK+Tp0.NRP); %initialize X

X(:,1:Tp0.NK) = matrix_FFkernel(FFstimR, Tp0.frac, Tp0.NK, Tp0.dK, Tp0.latK ); %add the kernel portion of X (series of lagged stim vectors)
%  Nonmodern version: X(:,1:Tp.NK) = zmatrix_FFkernel( FFstim, Tp );

%%% Spike history (and spikes-process)
[Xt spksN] = matrix_spike_history( FFspksR, Tp0.dRP, Tp0.stimres/1000/Tp0.frac, NT );

K(1:Tp0.NK) = fitH.K;
K(Tp0.NK+(1:length(fitH.wRP))) = fitH.wRP;
K(end) = fitH.offset;

k = K(1:end-1);
b = K(end);

r = log(1+exp(k*X(:,1:end-1)'+b));

T = Tp0.stimres/1000/Tp0.frac*(1:NT);
tstim = (1:length(FFstimR))*Tp0.stimres/1000;
stim_up = interp1(tstim,FFstimR,T,'nearest');

%% PLOT COMPARISON TO REPEAT TRIALS
raster(FFspksR,[1 3])
hold on
plot(T,-15*r,'r')
plot(T,stim_up*3-30,'k')
xlim([1.2 1.4])

tres = Tp0.stimres/Tp0.frac;
tkern = tres*(0:(length(fp.kern)-1));
figure
plot(tkern,fp.kern)

%%
% % %% also
% % Tp0 = paramsFF( 8, 4, 1 ); %(frac, KorDK, dW, dRP, MAX_RP, period, stimres)
% 
% % %elimate a random subset of spks
% % to_elim = ceil(rand(ceil(length(FFspks)*2/3),1)*ceil(length(FFspks)));
% % FFspks(to_elim) = [];
% 
% % fits model (FFstim, spks, fit0, target, Wmods, lambdas, Tp, silent, uncon_wRP, Xmods )
% fit0 = RGfit_h( FFstim, FFspks, [], [], [], [10 0], Tp0 );
% 
% %%
% fit0.lambdas(1) = 100;
% fit0 = RGfit_h( FFstim, FFspks, fit0 );
% 
% %%
% fit0.lambdas(3) = 10;
% TpH = Tp0;
% TpH.dRP = [0 1 2 3 4 5 6 7 8 10 12 16 20 24 30 40 50 60];
% TpH.NRP = 17;
% fitH = RGfit_h( FFstim, FFspks, fit0, [], [], [], TpH );
% fp = RGfitprocess( fitH);
% 
% %% to be continued