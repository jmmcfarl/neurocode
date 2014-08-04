close all
clear all
%% SET PARAMETERS
L = 10; %length in deg
T = 10; %duration of sim in sec
Fs = 50; %spatial sampling freq (1/deg)
Ft = 100; %temporal sampling freq (Hz)

sf1 = 1; %spatial freq 1 (cyc/deg)
sf2 = 2; %spatial freq 2 (cyc/deg)
sf3 = 3; %spatial freq 2 (cyc/deg)

%phase offsets
phase0_1 = 0;
phase0_2 = pi/3;
phase0_3 = pi/5;
amp_1 = 1; %amplitude of first grating
amp_2 = 1; %amplitude of second grating
amp_3 = 1; %amplitude of second grating

tf = 4; %temporal freq (Hz)
%% CREATE GRATINGS
xx = 0:1/Fs:L;
tt = 0:1/Ft:T;

phases_1 = bsxfun(@plus,xx*sf1*2*pi + phase0_1,tt'*tf*2*pi);
grate_1 = amp_1*sin(phases_1);

phases_2 = bsxfun(@plus,xx*sf2*2*pi + phase0_2,tt'*tf*2*pi);
grate_2 = amp_2*sin(phases_2);

phases_3 = bsxfun(@plus,xx*sf3*2*pi + phase0_3,tt'*tf*2*pi);
grate_3 = amp_3*sin(phases_3);

grate_sum = grate_1 + grate_2 + grate_3;

% figure
% subplot(2,2,1)
% imagesc(xx,tt,grate_1); colormap(gray); set(gca,'ydir','normal');
% xlabel('Position');ylabel('Time');title('Grating 1');
% subplot(2,2,2)
% imagesc(xx,tt,grate_2); colormap(gray); set(gca,'ydir','normal');
% xlabel('Position');ylabel('Time');title('Grating 2');
% subplot(2,2,3)
% imagesc(xx,tt,grate_3); colormap(gray); set(gca,'ydir','normal');
% xlabel('Position');ylabel('Time');title('Grating 3');
% subplot(2,2,4)
% imagesc(xx,tt,grate_sum); colormap(gray); set(gca,'ydir','normal');
% xlabel('Position');ylabel('Time');title('Grating Sum');

%% VIEW ANIMATION OF GRATING COMPONENTS 
use_inds = find(xx < 4);
figure;
for ii = 1:length(tt)
% subplot(4,1,1)
% plot(xx(use_inds),grate_1(ii,use_inds));
% subplot(4,1,2)
% plot(xx(use_inds),grate_2(ii,use_inds));
% subplot(4,1,3)
% plot(xx(use_inds),grate_3(ii,use_inds));
% subplot(4,1,4)
clf;plot(xx(use_inds),grate_sum(ii,use_inds));
hold on
plot(1,0,'ro');
ylim([-3 3]);
pause(1/Ft);
end

%% ANIMATE 2d grating in time
close all
use_inds = find(xx < 2);
figure;colormap(gray);
grate_portion = grate_sum(:,use_inds);
for ii = 1:length(tt)
    imagesc(repmat(grate_portion(ii,:),length(use_inds),1));caxis([-3 3]);
    pause(1/Ft)
end

%% compute spatiotemporal power spectra
window = bsxfun(@times,hamming(length(xx))',hamming(length(tt)));
pow_spec_1 = fftshift(abs(fft2(grate_1.*window)));
pow_spec_2 = fftshift(abs(fft2(grate_2.*window)));
pow_spec_3 = fftshift(abs(fft2(grate_3.*window)));
pow_spec_sum = fftshift(abs(fft2(grate_sum.*window)));
ff_t = linspace(-Ft/2,Ft/2,length(tt));
ff_x = linspace(-Fs/2,Fs/2,length(xx));

figure
subplot(2,2,1)
imagesc(ff_x,ff_t,log10(pow_spec_1)); colormap(gray); xlim([0 Fs/2]); ylim([-Ft/2 Ft/2]); caxis([1 3]); set(gca,'ydir','normal');
xlabel('Spatial freq');ylabel('Temporal freq');title('Grating 1');
subplot(2,2,2)
imagesc(ff_x,ff_t,log10(pow_spec_2)); colormap(gray); xlim([0 Fs/2]); ylim([-Ft/2 Ft/2]);caxis([1 3]);set(gca,'ydir','normal');
xlabel('Spatial freq');ylabel('Temporal freq');title('Grating 2');
subplot(2,2,3)
imagesc(ff_x,ff_t,log10(pow_spec_3)); colormap(gray); xlim([0 Fs/2]); ylim([-Ft/2 Ft/2]);caxis([1 3]);set(gca,'ydir','normal');
xlabel('Spatial freq');ylabel('Temporal freq');title('Grating 3');
subplot(2,2,4)
imagesc(ff_x,ff_t,log10(pow_spec_sum)); colormap(gray); xlim([0 Fs/2]); ylim([-Ft/2 Ft/2]);caxis([1 3]);set(gca,'ydir','normal');
xlabel('Spatial freq');ylabel('Temporal freq');title('Grating Sum');

%% compute spatial power spectra in each frame
spatial_powspectra_1 = fftshift(abs(fft(bsxfun(@times,grate_1,hamming(length(xx))'),[],2)));
spatial_powspectra_2 = fftshift(abs(fft(bsxfun(@times,grate_2,hamming(length(xx))'),[],2)));
spatial_powspectra_3 = fftshift(abs(fft(bsxfun(@times,grate_3,hamming(length(xx))'),[],2)));
spatial_powspectra_sum = fftshift(abs(fft(bsxfun(@times,grate_sum,hamming(length(xx))'),[],2)));

[~,sf1_ind] = min(abs(ff_x-sf1));
[~,sf2_ind] = min(abs(ff_x-sf2));
[~,sf3_ind] = min(abs(ff_x-sf3));

figure
subplot(2,2,1)
imagesc(ff_x,tt,log10(spatial_powspectra_1)); colormap(gray); xlim([0 Fs/2]); caxis([0 3]); set(gca,'ydir','normal');
xlabel('Spatial freq');ylabel('Time');title('Grating 1');
subplot(2,2,2)
imagesc(ff_x,tt,log10(spatial_powspectra_2)); colormap(gray); xlim([0 Fs/2]); caxis([0 3]); set(gca,'ydir','normal');
xlabel('Spatial freq');ylabel('Time');title('Grating 2');
subplot(2,2,3)
imagesc(ff_x,tt,log10(spatial_powspectra_3)); colormap(gray); xlim([0 Fs/2]); caxis([0 3]); set(gca,'ydir','normal');
xlabel('Spatial freq');ylabel('Time');title('Grating 3');
subplot(2,2,4)
imagesc(ff_x,tt,log10(spatial_powspectra_sum)); colormap(gray); xlim([0 Fs/2]); caxis([0 3]); set(gca,'ydir','normal');
xlabel('Spatial freq');ylabel('Time');title('Grating Sum');

%plot the time course of power in each fundamental spatial freq
figure
plot(tt,log10(spatial_powspectra_sum(:,sf1_ind)));
hold on
plot(tt,log10(spatial_powspectra_sum(:,sf2_ind)),'r');
plot(tt,log10(spatial_powspectra_sum(:,sf3_ind)),'k');
