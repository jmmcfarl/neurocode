
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% David P. Nguyen <dpnguyen@neurostat.mit.edu>
% September 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST SCRIPT
% Provides examples of how to use the code in this directory
%
%
% PART A) SIMULATE DATA
%
%
% PART B) FREQ ESTIMATION
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART A)
%
% SAMPLING FREQ 200 Hz

simopt.T = 20; % sec
simopt.fs   = 200;  % Hz

%% compute trajectory for frequency
simopt.ft.f = [30 70 30];
simopt.ft.t = [0 5 20];

%% compute trajectory for amplitude
simopt.at.a = [1];
simopt.at.t = [0 20];

%% compute trajectory for noise
simopt.nt.n = [0.1];
simopt.nt.t = [0 20];

y = sine_amfm(simopt);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART B)
%
% FREQ ESTIMATION

[An,e] = dnf_arkal(y, simopt.fs, 'ARorder', 2);
plnf = dnf_arkal_coef2poleinfo(An);

figure(1);clf;
plot(plnf.polephase(2,:).*(1/pi).*(100));
title('Frequency tracking');

specgram = dnf_arkal_specgram(An, 'fs', 200, 'psdres', 100);
figure(2);clf;
imagesc(specgram.timestamp, specgram.w, specgram.specgram);
caxis([0 20]);
axis xy;
ylabel('Freq (Hz)');
title('Spectrogram: Adaptive filtering');


