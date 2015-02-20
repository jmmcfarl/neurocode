close all
clear all

save_dir = '~/Data/bruce/Expt_1_9_13_imfolder';
blackColor = 0;
whiteColor = 255;
cbounds = [5 95];
grayColor = (blackColor+whiteColor)/2;
Pix2Deg = 0.018837;
dt = 1/60;

niqf_x = 1/Pix2Deg/2;
niqf_t = 1/dt/2;

Nx = 400;
use_radius = Nx/2;
Nt = 250;

%%
close all
[fx,ft] = freqspace([Nt Nx]);
[Fx,Ft] = meshgrid(fx,ft);
fx = fx*niqf_x; Fx = Fx*niqf_x;
ft = ft*niqf_t; Ft = Ft*niqf_t;

win_smooth = 4;
lower_cutoff = [1 1];
upper_cutoff = [15 4];
win = fspecial('gaussian',[Nx Nt],win_smooth);

A = atan2(Ft,Fx);
X = [Fx(:) Ft(:)];
Hd = ones(Nt,Nx);
Hd(X(:,1).^2/lower_cutoff(1)^2 + X(:,2).^2/lower_cutoff(2)^2 < 1) = 0;
Hd(X(:,1).^2/upper_cutoff(1)^2 + X(:,2).^2/upper_cutoff(2)^2 > 1) = 0;

% ori_range = ([-45 45]+45)*pi/180;
% Hdmult = (A > ori_range(1) & A < ori_range(2));
% Hd = Hd.*Hdmult;

win = win./max(win(:));
Hd = conv2(Hd,win,'same');
noise = randn(Nt,Nx);
FT=Hd.*fftshift(fft2(noise));
noise=real(ifft2(ifftshift(FT)));
noise = noise';

% %original method of filtering
% win = fspecial('gaussian',[Nx Nt],10);
% win = win./max(win(:));
% h = fwind2(Hd',win);
% noise = randn(Nx,Nt);
% noise = filter2(h,noise,'same');

%normalize 
noise = noise/std(noise(:));

%set up constraints for circular aperture
xax = (1:Nx); xax = xax - mean(xax);
[XX,YY] = meshgrid(xax);
RR = sqrt(XX.^2 + YY.^2);
use_pix = RR < use_radius;

%do a rotation of the images of desired
or = 0;
noise_ims = repmat(noise,[1 1 Nx]);
noise_ims = permute(noise_ims,[1 3 2]);
if or ~= 0
    for ii = 1:Nt
        noise_ims(:,:,ii) = imrotate(noise_ims(:,:,ii),or,'bilinear','crop');
    end
end

%enforce circular aperture
noise_ims = reshape(noise_ims,[],Nt);
use_pix = use_pix(:);
noise_ims(~use_pix,:) = 0;
noise_ims = reshape(noise_ims,Nx,Nx,Nt);


% sat_val = prctile(abs(noise(:)),90);
% figure;colormap(gray)
% for ii = 1:Nt
%    imagesc(squeeze(noise_ims(:,:,ii))');caxis([-sat_val sat_val])
%    pause(0.01);
% end
% 


%%
cd ~/Desktop/
writerObj = VideoWriter('newfile.avi');
writerObj.FrameRate = 60;
open(writerObj);

figure;colormap(gray)
sat_val = prctile(abs(noise(:)),90);
for ii = 1:Nt
   imagesc(squeeze(noise_ims(:,:,ii))');caxis([-sat_val sat_val])
frame = getframe;
writeVideo(writerObj,frame);
end
close(writerObj);