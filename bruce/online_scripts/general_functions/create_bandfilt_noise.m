function [noise,w,Fx,Fy] = create_bandfilt_noise(siz,f0,a,orange)

nx = siz(2); ny = siz(1);

Fy = ((0:ny-1)-floor(ny/2))*(2/(ny));
Fx = ((0:nx-1)-floor(nx/2))*(2/(nx));
[fx,fy] = meshgrid(Fx,Fy);

% a = 0.2;
% f0 = 0.3;

f_p=atan2(fy,fx)*180/pi; % direction, -180 to 180 deg.
f_m=sqrt(fx.*fx+fy.*fy); % radial frequency squared

w=1/f0/a^3*exp(-2*pi/a^2*(f_m-f0).^2); % start with all-pass tefilter
d=(f_p < orange(1) | f_p > orange(2)); % find out-of-band frequencies
w(d) = 0;

zf = find(fx==0 & fy==0);
w(zf) = 0;

noise = randn(ny,nx);
% noise = trnd(5,ny,nx);
ft=w.*fftshift(fft2(noise));
noise=real(ifft2(ifftshift(ft)));
noise = noise/std(noise(:));