clear 
close all
cd /Users/James/Data/bruce/2_27_12/ExptA
imID = 54;
if imID <10
    filename = ['ExptA0000' num2str(imID) '.png'];
elseif imID >=10 && imID <100
    filename = ['ExptA000' num2str(imID) '.png'];
elseif imID >= 100
    filename = ['ExptA00' num2str(imID) '.png'];
end

IMAGEorg = imread(filename);
IMAGEorg = double(IMAGEorg); % convert to double format
IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]

[Ny Nx] = size(IMAGEorg);

IMAGE = IMAGEorg;
% dsfrac = 4;
dsfrac = 1;
IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
[Nyd Nxd] = size(IMAGE); %down-sampled image size
% extra = Nxd - 256;
% extra = Nxd - 512;
extra = Nxd - 1024;
IMAGE(:,[1:extra/2 end-extra/2+1:end]) = []; %crop to 256x256
figure
showIm(IMAGE, 'auto', 'full', 'Original texture');

%%
ImSize = size(IMAGE);

RandomPhase = angle(fft2(rand(ImSize(1), ImSize(2))));
%generate random phase structure

ImFourier = fft2(IMAGE);
%Fast-Fourier transform
Amp = abs(ImFourier);
%amplitude spectrum
Phase = angle(ImFourier);
%phase spectrum
Phase = Phase + RandomPhase;
%add random phase to original phase
ImScrambled = ifft2(Amp.*exp(sqrt(-1)*(Phase)));
%combine Amp and Phase then perform inverse Fourier

ImScrambled = real(ImScrambled); %get rid of imaginery
figure
showIm(ImScrambled, 'auto', 'full', 'Original texture');

figure
% temp = spatialPattern([512 512],-2);
temp = spatialPattern([1024 1024],-2);
showIm(temp, 'auto', 'full', 'Original texture');
%%
cd /Users/James/Data/bruce/2_27_12/

addpath(genpath('~/James_scripts/textureSynth/'));
addpath(genpath('~/James_scripts/matlabPyrTools/'));
Nsc = 6; % Number of scales
Nor = 4; % Number of orientations
Na = 5;  % Spatial neighborhood is Na x Na coefficients
	 % It must be an odd number!

params = textureAnalysis(IMAGE, Nsc, Nor, Na);

Niter = 10;	% Number of iterations of synthesis loop
% Nsx = 256;	% Size of synthetic image is Nsy x Nsx
% Nsy = 256;	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)
% Nsx = 512;	% Size of synthetic image is Nsy x Nsx
% Nsy = 512;	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)
Nsx = 1024;	% Size of synthetic image is Nsy x Nsx
Nsy = 1024;	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)

res = textureSynthesis(params, [Nsy Nsx], Niter);

% close all
figure
showIm(res, 'auto', 'full', 'Synthesized texture');
