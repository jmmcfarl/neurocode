clear all
close all

Pix2Deg = 0.018837;
Xsize = 6.01; %in deg
Ysize = 6.01; %in deg
xax = -Xsize/2:Pix2Deg:Xsize/2;
yax = -Ysize/2:Pix2Deg:Ysize/2;
[XX,YY] = meshgrid(xax,yax);

sigma_mean = 0.065;
sigma_std = 0.02;
sigma_min = 0.035;
sigma_max = 0.125;

n_frames = 500;
avg_n_gabors = 5*Xsize*Ysize;

blackColor = 0;
whiteColor = 255;
cbounds = [2 98];
% cbounds = [0 100];
thresh = 0.75;

print_dir = '/home/james/James_scripts/bruce/stim_generation/gabor_noise';
%%
cd(print_dir);
% gabor_sequence = zeros(n_frames,length(yax),length(xax));
for jj = 1:n_frames
% jj=1
fprintf('Frame %d of %d\n',jj,n_frames);
    n_gabors = poissrnd(avg_n_gabors);
    
    x0 = rand(n_gabors,1)*Xsize-Xsize/2;
    y0 = rand(n_gabors,1)*Ysize-Ysize/2;
    
    sigma = randn(n_gabors,1)*sigma_std+sigma_mean;
    sigma(sigma < sigma_min) = sigma_min; sigma(sigma > sigma_max) = sigma_max;

%     sigma = 0.14*ones(n_gabors,1);
    bandwidth = 0.3;
    lambda = sigma/bandwidth;
    theta = rand(n_gabors,1)*2*pi;
    phase = rand(n_gabors,1)*2*pi;
    
    gabor_noise = zeros(size(XX));
    for ii = 1:n_gabors
        Gfun = generate_gabor(XX,YY,x0(ii),y0(ii),theta(ii),sigma(ii),phase(ii),lambda(ii));
        Gfun = Gfun/max(abs(Gfun(:)));
        gabor_noise = gabor_noise + Gfun;
    end
    
    %     bounds = prctile(gabor_noise(:),cbounds);
    bounds = [-thresh thresh];
    gabor_noise = (gabor_noise-bounds(1))/range(bounds);
    gabor_noise = gabor_noise*whiteColor;
    gabor_noise(gabor_noise < blackColor) = blackColor;
    gabor_noise(gabor_noise > whiteColor) = whiteColor;
        
    gabor_png = gabor_noise/whiteColor;
    fname = sprintf('Gabor_noise_%.3d.pgm',jj);
    imwrite(gabor_png,fname,'pgm');
    close
    
    gabor_sequence(jj,:,:) = gabor_noise;
end

%%
tgabor_sequence = gabor_sequence;
thresh = 0.75;
tgabor_sequence(tgabor_sequence > thresh) = thresh; tgabor_sequence(tgabor_sequence < -thresh) = -thresh;

%%
n_frames = size(tgabor_sequence,1);
close all
imagesc(xax,yax,squeeze(tgabor_sequence(1,:,:))); colormap(gray)
pause
cloc = 0;tic
for jj = 1:n_frames
     clf
    fprintf('Frame %d of %d at %.4f sec\n',jj,n_frames,toc);
    imagesc(xax,yax,squeeze(tgabor_sequence(jj,:,:))); colormap(gray)
    pause(0.02)
end

%%
clear all

load examp_gabor_noise

n_frames = size(tgabor_sequence,1);
close all
imagesc(xax,yax,squeeze(tgabor_sequence(1,:,:))); colormap(gray)
pause
cloc = 0;tic
for jj = 1:n_frames
     clf
    fprintf('Frame %d of %d at %.4f sec\n',jj,n_frames,toc);
    imagesc(xax,yax,squeeze(tgabor_sequence(jj,:,:))); colormap(gray)
    pause(0.02)
end

%%
 F = linspace(-0.5,0.5,length(xax));
filt = fspecial('gaussian',10,2);
pspec = abs(fftshift(fft2(gabor_noise)));
pspec = imfilter(pspec,filt);