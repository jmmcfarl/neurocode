% close all
clear all

fig_dir = '/home/james/Analysis/bruce/variability/stoch_proc_sim/';

stim_type = 'grating';
% stim_type = 'rls';
%%
if strcmp(stim_type,'grating')
gr_sf = 4;
gr_tf = 4;

pix_width = .01;
usfac = 1;
pix_dx = pix_width*usfac;
xax = -2.5:pix_dx:2.5;
nsamps = 5e3; %number of random pattern samples ('trials')
dt = .03;
tax = (1:nsamps)*dt;

[XX,TT] = meshgrid(xax,tax);
Xstim_up = sin(2*pi*(gr_sf*XX + gr_tf*TT)); %base (static EP) drift grate
Xstim_up = Xstim_up';
end
%% make spatial stim
if strcmp(stim_type,'rls')
pix_width = 0.0565; %bar width in deg
npix = round(5/pix_width); %number of bars 
nsamps = 5e3; %number of random pattern samples ('trials')
dds = 67; %dot density
usfac = 5; %spatial up-sampling factor for calculations
npix_us = npix*usfac; %number of pixels in up-sampled stim

%create RLS stim
Xstim = randi(2,nsamps,npix);
Xstim(Xstim == 2) = -1;
Xiszero = rand(nsamps,npix);
Xstim(Xiszero > dds/100) = 0;

%spatial up-sampling
if usfac > 1
    Xstim_up = zeros(nsamps,npix_us);
    for ii = 1:size(Xstim,2)
        for jj = 1:usfac
            Xstim_up(:,usfac*(ii-1)+jj) = Xstim(:,ii);
        end
    end
elseif usfac == 1
    Xstim_up = Xstim;
end

%convert to dims [space X trials]
Xstim_up = Xstim_up';

%position axis for stimulus
xax = (1:npix_us)*pix_width/usfac;
xax = xax - mean(xax);
end

%% create gaussian filter
sf = 4; %spatial freq
bandwidth = 1/2; 
env_sigma = 1/sf*bandwidth; %spatial envelope SD
spatial_phase = 0;

gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf + spatial_phase);

%% compute amplitude spectrum of stimulus
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));
params.Fs = 1/pix_width*usfac;
params.tapers = [1 1];
params.trialave = 1;

[P_FT,F] = mtspectrumc(Xstim_up,params);
A_FT = sqrt(P_FT); %convert from pow-spec to amp-spec 

%% get amplitude spectrum of the gabor filter (analytically)
fft_bandwidth = 1/(2*pi*env_sigma); %frequency-domain width (SD) of gabor
gabor_fft = exp(-(F-sf).^2/(2*fft_bandwidth^2)); %gaussian centered on SF

%% compute eye position distribution
ep_ax = -2:(pix_width/usfac):2;
eye_pos_sigma = 0.1;
eye_pos_fft_sigma = 1/(2*pi*eye_pos_sigma);
eye_pos_fft = exp(-F.^2/(2*eye_pos_fft_sigma^2));
ep_dist = exp(-ep_ax.^2/(2*eye_pos_sigma^2));

%% get fft of ep distribution
ep_Fs = 1/median(diff(ep_ax));
ep_fft = 2*abs(fftshift(fft(ep_dist))); %one-sided amplitude spec
ep_fax = linspace(-ep_Fs/2,ep_Fs/2,length(ep_dist)); %eye position dist frequency axis

ep_dist_convfilt = ep_dist/sum(ep_dist); %normalize for convolution

%% compute firing rate outputs of these models
%spk NL parameters
spkNL_alpha = 1;
spkNL_beta = 2;
spkNL_theta = 0;

%output of gabor filter (complex-valued)
Xconv = convn(Xstim_up,gabor_filt(:),'same');
Xconv = Xconv/std(real(Xconv(:))); %normalize to unit SD

%firing rate output of simple-cell model
simple_rate = spkNL_alpha*log(1+exp((real(Xconv)+spkNL_theta)*spkNL_beta));
simple_srate = std(simple_rate(:));
simple_rate = simple_rate/simple_srate;

%energy model output
Xout = (real(Xconv).^2 + imag(Xconv).^2);
complex_rate = spkNL_alpha*log(1+exp((Xout+spkNL_theta)*spkNL_beta));
complex_srate = std(complex_rate(:));
complex_rate = complex_rate/complex_srate;

%%
simple_psth = convn(simple_rate,ep_dist_convfilt(:),'same');
complex_psth = convn(complex_rate,ep_dist_convfilt(:),'same');
%% compute power spectra of different model outputs
% %for gabor filter output
% [Pxc_FT,F] = mtspectrumc(real(Xconv),params);
% Axc_FT = sqrt(Pxc_FT); %convert from pow-spec to amp-spec 

%simple cell firing rate output
[Prate,F] = mtspectrumc(simple_rate - mean(simple_rate(:)),params);
Arate = sqrt(Prate); %convert from pow-spec to amp-spec 

%energy model firing rate output
[Prate_EN,F] = mtspectrumc(complex_rate - mean(complex_rate(:)),params);
Arate_EN = sqrt(Prate_EN); %convert from pow-spec to amp-spec 

%simple cell PSTH
[Ppsth,F] = mtspectrumc(simple_psth - mean(simple_psth(:)),params);
Apsth = sqrt(Ppsth); %convert from pow-spec to amp-spec 

%energy model PSTH
[Ppsth_EN,F] = mtspectrumc(complex_psth - mean(complex_psth(:)),params);
Apsth_EN = sqrt(Ppsth_EN); %convert from pow-spec to amp-spec 

%% calculate alphas
simple_tot_var = var(simple_rate(:));
simple_psth_var = var(simple_psth(:));
energy_tot_var = var(complex_rate(:));
energy_psth_var = var(complex_psth(:));

fprintf('Simple alpha: %.4f\n',1-simple_psth_var/simple_tot_var);
fprintf('Energy alpha: %.4f\n',1-energy_psth_var/energy_tot_var);

%% plot amplitude spectra
close all
freq_range = [0 10];

f1 = figure();
subplot(2,1,1);
plot(F,A_FT/max(A_FT));
hold on
plot(ep_fax,ep_fft/max(ep_fft),'k');
plot(F,gabor_fft/max(gabor_fft),'r')
plot(F,Arate/max(Arate),'m')
plot(F,Apsth/max(Arate),'m--');
xlim(freq_range);
xlabel('Frequency (cyc/deg)');
ylabel('Relative amplitude');

subplot(2,1,2);
plot(F,A_FT/max(A_FT));
hold on
plot(ep_fax,ep_fft/max(ep_fft),'k');
plot(F,gabor_fft/max(gabor_fft),'r')
plot(F,Arate_EN/max(Arate_EN),'m')
plot(F,Apsth_EN/max(Arate_EN),'m--');
xlim(freq_range);
xlabel('Frequency (cyc/deg)');
ylabel('Relative amplitude');


% fig_width = 4; rel_height = 1.6;
% figufy(f1);
% fname = [fig_dir sprintf('sim_pspecs_%s.pdf',stim_type)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% plot example rate(E) functions
f2 = figure();
subplot(2,1,1);
plot(xax,simple_rate(:,1));
hold on
plot(xax,simple_psth(:,1),'k');
xlabel('Eye position');
ylabel('Rate');
xlim([-1 1]);

subplot(2,1,2);
plot(xax,complex_rate(:,1));
hold on
plot(xax,complex_psth(:,1),'k');
xlabel('Eye position');
ylabel('Rate');
xlim([-1 1]);

fig_width = 4; rel_height = 1.6;
figufy(f2);
fname = [fig_dir sprintf('sim_ratefuns_%s.pdf',stim_type)];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);


%% eye position distribution
f3 = figure();
plot(ep_ax,ep_dist,'k');
xlim([-0.5 0.5]);
xlabel('Eye position');
ylabel('Probability');

fig_width = 4; rel_height = 0.8;
figufy(f3);
fname = [fig_dir 'EP_dist.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

%% plot gabor filter
f4 = figure();
plot(xax,real(gabor_filt),'k');
hold on
plot(xax,imag(gabor_filt),'r');
xlim([-0.5 0.5]);
xlabel('Eye position');
ylabel('Filter amplitude');
xl = xlim(); line(xl,[0 0],'color','k','linestyle','--');

fig_width = 4; rel_height = 0.8;
figufy(f4);
fname = [fig_dir 'Gabor_filt.pdf'];
exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%%
if strcmp(stim_type,'grating');
examp_grating = Xstim_up(:,1);
examp_grating = repmat(examp_grating,1,length(examp_grating));

f5 = figure();
imagesc(xax,xax,examp_grating');
colormap(gray);
xlim([-1 1]); 
ylim([-1 1]);
set(gca,'xtick',[],'ytick',[]);

fig_width = 4; rel_height = 0.8;
figufy(f5);
fname = [fig_dir 'Grating_stim.pdf'];
exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f5);

end

%% estimate dependence of alpha on EP SD and grating frequency
grating_SFs = linspace(0,6,100);
sigma_vals = 0.05:.025:0.2;

f6 = figure(); hold on

cmap = jet(length(sigma_vals));
for ii = 1:length(sigma_vals)
    sigma = sigma_vals(ii);
    alphas = 1 - (normpdf(grating_SFs,0,1/(2*pi*sigma))/sqrt(2*pi)/sigma).^2;
    
    plot(grating_SFs,alphas,'color',cmap(ii,:),'linewidth',2);
end
    
ylim([0 1]);
legendCell = cellstr(num2str(sigma_vals', 'N=%.3f'));
legend(legendCell,'location','southeast')
xlabel('Spatial frequency (cyc/deg)');
ylabel('Alpha');
