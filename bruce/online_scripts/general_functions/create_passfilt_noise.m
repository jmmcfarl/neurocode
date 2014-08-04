function noise = create_passfilt_noise(siz,frange,orange)

nx = siz(2); ny = siz(1);

fy = ((0:ny-1)-floor(ny/2))*(2/(ny));
fx = ((0:nx-1)-floor(nx/2))*(2/(nx));
[fx,fy] = meshgrid(fx,fy);

f_p=atan2(fy,fx)*180/pi; % direction, -180 to 180 deg.
f_m=sqrt(fx.*fx+fy.*fy); % radial frequency squared

w=ones(size(fx)); % start with all-pass filter
d=(f_m<frange(1) | f_m>frange(2) | f_p < orange(1) | f_p > orange(2)); % find out-of-band frequencies
w(d) = 0;

noise = randn(ny,nx);
ft=w.*fftshift(fft2(noise));
noise=real(ifft2(ifftshift(ft)));
noise = noise/std(noise(:));