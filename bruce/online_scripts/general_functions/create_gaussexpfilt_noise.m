function [noise,w,Fx,Fy] = create_gaussexpfilt_noise(siz,f0,or0,a,dec_const)

nx = siz(2); ny = siz(1);

Fy = ((0:ny-1)-floor(ny/2))*(2/(ny));
Fx = ((0:nx-1)-floor(nx/2))*(2/(nx));
[fx,fy] = meshgrid(Fx,Fy);

fx0 = f0*cos(or0);
fy0 = f0*sin(or0);
% f_p=atan2(fy,fx)*180/pi; % direction, -180 to 180 deg.
f_m=sqrt(fx.*fx+fy.*fy); % radial frequency squared

w = exp(-((fx-fx0).^2+(fy-fy0).^2)/(2*a^2));
w = w.*(1-exp(-f_m/dec_const));

noise = randn(ny,nx);
ft=w.*fftshift(fft2(noise));
noise=real(ifft2(ifftshift(ft)));
noise = noise/std(noise(:));