Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
siz = [1280 1280];
Fs = 1/Pix2Deg;
niqf = Fs/2;
xax = linspace(-Nxp/2,Nxp/2,Nxp)/Fs; yax = linspace(-Nyp/2,Nyp/2,Nyp)/Fs;
[XAX,YAX] = meshgrid(xax,yax);

%%
freq_cent = 4.5; %cyc/deg
freq_a = 9;
or_width = 35;

or_cent = 45;
[noise,w,Fx,Fy] = create_bandfilt_noise(siz,freq_cent/niqf,freq_a/niqf,[-or_width or_width]/2+or_cent);
noise(Nyp+1:end,:) = [];

or_cent = -45;
[noise2,w,Fx,Fy] = create_bandfilt_noise(siz,freq_cent/niqf,freq_a/niqf,[-or_width or_width]/2+or_cent);
Fx = Fx*niqf;
Fy = Fy*niqf;
noise2(Nyp+1:end,:) = [];

% dec_const = 2/niqf;
% [noise2,w2,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,orange,dec_const);
% Fx = Fx*niqf;
% Fy = Fy*niqf;
% noise2(Nyp+1:end,:) = [];

%%
% freq_cent = 4.5; %cyc/deg
% gauss_a = 1.5;
% or0 = 45;
% [noise,w,Fx,Fy] = create_gaussfilt_noise(siz,freq_cent/niqf,or0,gauss_a/niqf);
% noise(Nyp+1:end,:) = [];
% 
% freq_cent = 4.5; %cyc/deg
% gauss_a = 1.5;
% or0 = -45;
% [noise2,w,Fx,Fy] = create_gaussfilt_noise(siz,freq_cent/niqf,or0,gauss_a/niqf);
% Fx = Fx*niqf;
% Fy = Fy*niqf;
% noise2(Nyp+1:end,:) = [];

freq_cent = 4.5; %cyc/deg
gauss_a = 2;
or0 = 45;
dec_const = 2/niqf;
[noise,w,Fx,Fy] = create_gaussexpfilt_noise(siz,freq_cent/niqf,or0,gauss_a/niqf,dec_const);
noise(Nyp+1:end,:) = [];

freq_cent = 4.5; %cyc/deg
gauss_a = 2;
or0 = -45;
dec_const = 2/niqf;
[noise2,w2,Fx,Fy] = create_gaussexpfilt_noise(siz,freq_cent/niqf,or0,gauss_a/niqf,dec_const);
Fx = Fx*niqf;
Fy = Fy*niqf;
noise2(Nyp+1:end,:) = [];

%%
% fore_cent = [0.35 -0.4];
blackColor = 0;
whiteColor = 255;
cbounds = [5 95];
fore_cent = [-2 0];
fore_width = [3 3];
foreground = ones(size(noise));
d = XAX > fore_cent(1) + fore_width(1)/2 | XAX < fore_cent(1)-fore_width(1)/2 | ...
    YAX > fore_cent(2) + fore_width(2)/2 | YAX < fore_cent(2) - fore_width(2)/2;
foreground(d) = 0;

text_image = noise;
text_image(foreground==1) = noise2(foreground==1);

bounds = prctile(text_image(:),cbounds);
text_image = (text_image-bounds(1))/range(bounds);
text_image = text_image*whiteColor;
text_image(text_image < blackColor) = blackColor;
text_image(text_image > whiteColor) = whiteColor;
