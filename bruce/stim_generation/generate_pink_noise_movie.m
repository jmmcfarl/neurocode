close all
clear all

save_dir = '~/Data/bruce/isotropic_noise_test';
blackColor = 0;
whiteColor = 255;
% cbounds = [25 75];
cbounds = [5 95];
% cbounds = [10 90];
grayColor = (blackColor+whiteColor)/2;
Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
siz = [1280 1280];
Fs = 1/Pix2Deg;
niqf = Fs/2;
xax = linspace(-Nxp/2,Nxp/2,Nxp)/Fs; yax = linspace(-Nyp/2,Nyp/2,Nyp)/Fs;
[XAX,YAX] = meshgrid(xax,yax);


%%
space_std = 10;
temp_std = 10;
h = fspecial3('gaussian',[space_std space_std temp_std]);

%%
Npix = 200;
Nframes = 50;
noise = randn(Npix,Npix,Nframes);

filt_noise = imfilter(noise,h);

%%
imagesc(squeeze(filt_noise(:,:,1)));colormap(gray);
caxis([-0.04 0.04])
pause
for ff = 1:Nframes
    imagesc(squeeze(filt_noise(:,:,ff)));colormap(gray);
caxis([-0.04 0.04])
    pause(0.025)
    
end
