% close all
clear all

fig_dir = '/home/james/Analysis/bruce/variability/stoch_proc_sim/';
% fig_dir = '/Users/james/Analysis/bruce/variability/stoch_proc_sim/';

eye_pos_sigma = 0.11; %SD (deg) of eye position distribution (assume gaussian)

%spk NL parameters of model neuron
spkNL_alpha = 1;
spkNL_beta = 2;
spkNL_theta = 0.5;

%gabor params
bandwidth = 1/2;
spatial_phase = 0;

%% make spatial stim
rel_pwidth = 1; %pixels per bar
pix_width = 0.01; %bar width in deg
npix = round(10/pix_width); %number of bars
nsamps = 5e3; %number of random pattern samples ('trials')
dds = 67; %dot density
usfac = 1; %spatial up-sampling factor for calculations
npix_us = npix*usfac; %number of pixels in up-sampled stim
Xstim_up = generate_RLS_stim(nsamps,npix,dds,usfac);

%convert to dims [space X trials]
Xstim_up = Xstim_up';

%position axis for stimulus
xax = (1:npix_us)*pix_width/usfac;
xax = xax - mean(xax);
Fs = 1/(pix_width/usfac); %spatial sample freq
N = length(xax); %number of spatial samples
fax = 0:Fs/N:Fs/2; %frequency axis (single sided)

% compute eye position distribution
ep_dist = exp(-xax.^2/(2*eye_pos_sigma^2)); %EP distribution
ep_dist = ep_dist/sum(ep_dist);
ep_dist_convfilt = ep_dist/sum(ep_dist); %normalize for convolution
ep_dist_fft = abs(fft(ep_dist));
ep_dist_fft = ep_dist_fft(1:N/2+1);

%% run analysis using a range of different filter SFs
% poss_SFs = [0.25 0.5 1 2 3 4 5 6 7 8];
poss_SFs = [0.25 0.5:0.5:8];
[simple_alpha,complex_alpha,filt_alpha] = deal(nan(length(poss_SFs),1));
for ii = 1:length(poss_SFs)
    
    % make gabor filter
    fprintf('SF %d/%d\n',ii,length(poss_SFs));
    sf = poss_SFs(ii); %spatial freq
    env_sigma = 1/sf*bandwidth; %spatial envelope SD
    
    %make gabor filter
    gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf + spatial_phase);
    
    %output of gabor filter (complex-valued)
    Xconv = convn(Xstim_up,gabor_filt(:),'same');
    Xconv = Xconv/std(real(Xconv(:))); %normalize to unit SD
    
    %firing rate output of simple-cell model
    simple_rate = spkNL_alpha*log(1+exp((real(Xconv)+spkNL_theta)*spkNL_beta));
    
    %energy model output
    Xout = (real(Xconv).^2 + imag(Xconv).^2);
    complex_rate = spkNL_alpha*log(1+exp((Xout+spkNL_theta)*spkNL_beta));
    
    %subtract out avg rates
    simple_rate = simple_rate - mean(simple_rate(:));
    complex_rate = complex_rate - mean(complex_rate(:));
    
    
    filt_out = real(Xconv); 
    filt_fft = mean(abs(fft(filt_out)),2);
    simple_fft = mean(abs(fft(simple_rate)),2);
    complex_fft = mean(abs(fft(complex_rate )),2);
    filt_fft = filt_fft(1:N/2 + 1);
    simple_fft = simple_fft(1:N/2 + 1);
    complex_fft = complex_fft(1:N/2 + 1);
    
    % compute psths by convolving the model-rates with the EP distribution
    filt_psth = convn(filt_out,ep_dist_convfilt(:),'same');
    simple_psth = convn(simple_rate,ep_dist_convfilt(:),'same');
    complex_psth = convn(complex_rate,ep_dist_convfilt(:),'same');
    
    % calculate alphas
    simple_tot_var = var(simple_rate(:));
    simple_psth_var = var(simple_psth(:));
    energy_tot_var = var(complex_rate(:));
    energy_psth_var = var(complex_psth(:));
    
%     simple_alpha(ii) = 1 - simple_psth_var/simple_tot_var;
%     complex_alpha(ii) = 1 - energy_psth_var/energy_tot_var; 
    filt_alpha(ii) = 1 - trapz(fax,(filt_fft.*ep_dist_fft').^2)/trapz(fax,filt_fft.^2);
    simple_alpha(ii) = 1 - trapz(fax,(simple_fft.*ep_dist_fft').^2)/trapz(fax,simple_fft.^2);
    complex_alpha(ii) = 1 - trapz(fax,(complex_fft.*ep_dist_fft').^2)/trapz(fax,complex_fft.^2);
end

% plot filter SF vs alpha for simple and complex models
% close all
f1 = figure(); hold on
plot(poss_SFs,simple_alpha,'linewidth',2);
plot(poss_SFs,complex_alpha,'r','linewidth',2);
plot(poss_SFs,filt_alpha,'k','linewidth',2);
xlabel('Preferred spatial frequency (cyc/deg)');
ylabel('Alpha');

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'sim_alpha_vs_SF.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% now run analysis using an example gabor filter
examp_sf = 2; %spatial freq

% make spatial stim
rel_pwidth = 1;
pix_width = 0.01; %bar width in deg
npix = round(10/pix_width); %number of bars
nsamps = 1e4; %number of random pattern samples ('trials')
dds = 67; %dot density
usfac = rel_pwidth; %spatial up-sampling factor for calculations
npix_us = npix*usfac; %number of pixels in up-sampled stim

Xstim_up = generate_RLS_stim(nsamps,npix,dds,usfac);

%convert to dims [space X trials]
Xstim_up = Xstim_up';

%position axis for stimulus
xax = (1:npix_us)*pix_width/usfac;
xax = xax - mean(xax);
Fs = 1/(pix_width/usfac); %spatial sample freq
N = length(xax); %number of spatial samples
fax = 0:Fs/N:Fs/2; %frequency axis (single sided)

env_sigma = 1/examp_sf*bandwidth; %spatial envelope SD
gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*examp_sf + spatial_phase);

% get amplitude spectrum of the gabor filter (analytically)
fft_bandwidth = 1/(2*pi*env_sigma); %frequency-domain width (SD) of gabor
gabor_fft = exp(-(fax-examp_sf).^2/(2*fft_bandwidth^2)); %gaussian centered on SF

ep_dist = exp(-xax.^2/(2*eye_pos_sigma^2)); %EP distribution
ep_dist = ep_dist/sum(ep_dist);
ep_dist_convfilt = ep_dist/sum(ep_dist); %normalize for convolution
ep_dist_fft = abs(fft(ep_dist));
ep_dist_fft = ep_dist_fft(1:N/2+1);

% stim_fft = mean(abs(fftshift(fft(Xstim_up))),2);
stim_fft = mean(abs(fft(Xstim_up)),2);
stim_fft = stim_fft(1:N/2+1);

% compute firing rate outputs of these models
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

% compute psths through convolution with EP dist
simple_psth = convn(simple_rate,ep_dist_convfilt(:),'same');
complex_psth = convn(complex_rate,ep_dist_convfilt(:),'same');

%overall avg rates
simp_mrate = mean(simple_rate(:));
comp_mrate = mean(complex_rate(:));

% compute power spectra of different model outputs
%for gabor filter output
gout_fft = mean(abs(fft(real(Xconv))),2);
gout_fft = gout_fft(1:N/2+1);

%simple cell firing rate output
simprate_fft = mean(abs(fft(simple_rate - simp_mrate)),2);
simprate_fft = simprate_fft(1:N/2+1);

%energy model firing rate output
comprate_fft = mean(abs(fft(complex_rate - comp_mrate)),2);
comprate_fft = comprate_fft(1:N/2+1);

%simple cell PSTH
simppsth_fft = mean(abs(fft(simple_psth - simp_mrate)),2);
simppsth_fft = simppsth_fft(1:N/2+1);

%energy model PSTH
comppsth_fft = mean(abs(fft(complex_psth - comp_mrate)),2);
comppsth_fft = comppsth_fft(1:N/2+1);
  
ex_simple_alpha = 1 - trapz(fax,(simprate_fft.*ep_dist_fft').^2)/trapz(fax,simprate_fft.^2);
ex_complex_alpha = 1 - trapz(fax,(comprate_fft.*ep_dist_fft').^2)/trapz(fax,comprate_fft.^2);

%% plot amplitude spectra
freq_range = [0 6];

f1 = figure();
subplot(2,1,1);
hold on
% plot(fax,stim_fft/max(stim_fft),'b');
plot(fax,ep_dist_fft/max(ep_dist_fft),'k');
plot(fax,gabor_fft/max(gabor_fft),'r')
plot(fax,simprate_fft/max(simprate_fft),'m')
plot(fax,simppsth_fft/max(simprate_fft),'m--');
xlim(freq_range);
xlabel('Frequency (cyc/deg)');
ylabel('Relative amplitude');

subplot(2,1,2);
% plot(fax,stim_fft/max(stim_fft),'b');
hold on
plot(fax,ep_dist_fft/max(ep_dist_fft),'k');
plot(fax,gabor_fft/max(gabor_fft),'r')
plot(fax,comprate_fft/max(comprate_fft),'m')
plot(fax,comppsth_fft/max(comprate_fft),'m--');
xlim(freq_range);
xlabel('Frequency (cyc/deg)');
ylabel('Relative amplitude');


% fig_width = 4; rel_height = 1.6;
% figufy(f1);
% fname = [fig_dir 'sim_pspecs_grate.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% plot example rates and EP dist
xr = [-2 2];
ex_trial = 2;
close all

f1 = figure();
subplot(2,1,1);
hold on
plot(xax,simple_rate(:,ex_trial));
plot(xax,simple_psth(:,ex_trial),'r')
xlim(xr);
xlabel('Eye position (deg)');
ylabel('Firing rate');

subplot(2,1,2);
hold on
plot(xax,complex_rate(:,ex_trial));
plot(xax,complex_psth(:,ex_trial),'r')
xlim(xr);
xlabel('Eye position (deg)');
ylabel('Firing rate');

f2 = figure();
plot(xax,ep_dist,'k');
xlim(xr);
xlabel('Eye position');
ylabel('Probability');

fig_width = 4; rel_height = 1.6;
figufy(f1);
fname = [fig_dir 'sim_rates.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

fig_width = 4; rel_height = 0.8;
figufy(f2);
fname = [fig_dir 'sim_EP_dist.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% plot gabor filter
f4 = figure();
plot(xax,real(gabor_filt),'k');
hold on
plot(xax,imag(gabor_filt),'r');
xlim([-1 1]);
xlabel('Eye position');
ylabel('Filter amplitude');
xl = xlim(); line(xl,[0 0],'color','k','linestyle','--');

% fig_width = 4; rel_height = 0.8;
% figufy(f4);
% fname = [fig_dir 'Gabor_filt.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f4);
% 
%% analysis for grating stim
gr_sf = 2;
gr_tf = 4;

pix_dx = .01;
xax = -5:pix_dx:5;
if mod(length(xax),2) ~= 0
    xax(end) = [];
end
nsamps = 5e3; %number of random pattern samples ('trials')
dt = .01;
tax = (1:nsamps)*dt;
N = length(xax);
Fs = 1/pix_dx;
fax = 0:Fs/N:Fs/2; %frequency axis (single sided)

[XX,TT] = meshgrid(xax,tax);
Xstim_up = sin(2*pi*(gr_sf*XX + gr_tf*TT)); %base (static EP) drift grate
Xstim_up = Xstim_up';

env_sigma = 1/examp_sf*bandwidth; %spatial envelope SD
gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*examp_sf + spatial_phase);

% get amplitude spectrum of the gabor filter (analytically)
fft_bandwidth = 1/(2*pi*env_sigma); %frequency-domain width (SD) of gabor
gabor_fft = exp(-(fax-examp_sf).^2/(2*fft_bandwidth^2)); %gaussian centered on SF

ep_dist = exp(-xax.^2/(2*eye_pos_sigma^2)); %EP distribution
ep_dist = ep_dist/sum(ep_dist);
ep_dist_convfilt = ep_dist/sum(ep_dist); %normalize for convolution
ep_dist_fft = abs(fft(ep_dist));
ep_dist_fft = ep_dist_fft(1:N/2+1);

stim_fft = mean(abs(fft(Xstim_up)),2);
stim_fft = stim_fft(1:N/2+1);

% compute firing rate outputs of these models
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

% compute psths through convolution with EP dist
simple_psth = convn(simple_rate,ep_dist_convfilt(:),'same');
complex_psth = convn(complex_rate,ep_dist_convfilt(:),'same');

%overall avg rates
simp_mrate = mean(simple_rate(:));
comp_mrate = mean(complex_rate(:));

% compute power spectra of different model outputs
%for gabor filter output
gout_fft = mean(abs(fft(real(Xconv))),2);
gout_fft = gout_fft(1:N/2+1);

%simple cell firing rate output
simprate_fft = mean(abs(fft(simple_rate - simp_mrate)),2);
simprate_fft = simprate_fft(1:N/2+1);

%energy model firing rate output
comprate_fft = mean(abs(fft(complex_rate - comp_mrate)),2);
comprate_fft = comprate_fft(1:N/2+1);

%simple cell PSTH
simppsth_fft = mean(abs(fft(simple_psth - simp_mrate)),2);
simppsth_fft = simppsth_fft(1:N/2+1);

%energy model PSTH
comppsth_fft = mean(abs(fft(complex_psth - comp_mrate)),2);
comppsth_fft = comppsth_fft(1:N/2+1);
  
ex_simple_alpha = 1 - trapz(fax,(simprate_fft.*ep_dist_fft').^2)/trapz(fax,simprate_fft.^2);
ex_complex_alpha = 1 - trapz(fax,(comprate_fft.*ep_dist_fft').^2)/trapz(fax,comprate_fft.^2);



%     examp_grating = Xstim_up(:,1);
%     examp_grating = repmat(examp_grating,1,length(examp_grating));
%
%     f5 = figure();
%     imagesc(xax,xax,examp_grating');
%     colormap(gray);
%     xlim([-1 1]);
%     ylim([-1 1]);
%     set(gca,'xtick',[],'ytick',[]);
%
%     %     fig_width = 4; rel_height = 0.8;
%     %     figufy(f5);
%     %     fname = [fig_dir 'Grating_stim.pdf'];
%     %     exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
%     %     % close(f5);
%

%% estimate dependence of alpha on EP SD and grating frequency
grating_SFs = linspace(0,6,100);
% sigma_vals = 0.05:.025:0.2;
sigma = eye_pos_sigma;
ex_SF = 2;
sigma_vals = 0.0:.01:0.2;


alphas_SF = 1 - (normpdf(grating_SFs,0,1/(2*pi*sigma))/sqrt(2*pi)/sigma).^2;

grating_SF_set = [0.5 1 2 4];
alpha_SD = nan(length(sigma_vals),length(grating_SF_set));
for ii = 1:length(sigma_vals)
    alpha_SD(ii,:) = 1 - (normpdf(grating_SF_set,0,1/(2*pi*sigma_vals(ii)))/sqrt(2*pi)/sigma_vals(ii)).^2;
end
alpha_SD(sigma_vals == 0,:) = 0;

f1 = figure(); hold on
plot(grating_SFs,alphas_SF,'k','linewidth',2);
ylim([0 1]);
xlabel('Spatial frequency (cyc/deg)');
ylabel('Alpha');


f2 = figure(); hold on
line_widths = [0.5 1 2 3];
for jj = 1:length(grating_SF_set)
plot(sigma_vals,alpha_SD(:,jj),'k','linewidth',line_widths(jj));
end
ylim([0 1]);
xlabel('Eye pos SD (deg)');
ylabel('Alpha');

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Grating_alphas.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);
% 
% fig_width = 4; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Grating_SD.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);
% 
%%
% stim_type = 'white';
stim_type = 'NS';

sf = 2; %spatial freq
bandwidth = 0.5;
env_sigma = 1/sf*bandwidth; %spatial envelope SD
spatial_phase = 0;
spatial_AR = 0.75; %spatial aspect ratio
orientation = 45; %grating orientation
  
pix_dx = .01;
xax = -10:pix_dx:10;
[X,Y] = meshgrid(xax);
Xp = X*cosd(orientation) + Y*sind(orientation);
Yp = -X*sind(orientation) + Y*cosd(orientation);

gabor_filt = exp(-(Xp.^2 + spatial_AR*Yp.^2)/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*Xp*sf + spatial_phase);
gabor_filt = gabor_filt - mean(gabor_filt(:));

gabor_fft = abs(fftshift(fft2(real(gabor_filt))));
niqf = 1/pix_dx/2;
fax = linspace(-niqf,niqf,length(xax));

[Fx,Fy] = meshgrid(fax);
FF = sqrt(Fx.^2 + Fy.^2);
if strcmp(stim_type,'NS')
    stim_aspec = 1./FF;
    stim_aspec(FF==0) = nan; stim_aspec(FF==0) = nanmax(stim_aspec(:));
    stim_aspec = stim_aspec/max(stim_aspec(:));
elseif strcmp(stim_type,'white');
    stim_aspec = ones(size(FF));
end  

eye_pos_sigma_2d = sqrt(2*eye_pos_sigma.^2);
eye_pos_fft_sigma = 1/(2*pi*eye_pos_sigma_2d);
eye_pos_fft = exp(-(Fx.^2 + Fy.^2)/(2*eye_pos_fft_sigma^2));

gabor_output_spec = gabor_fft.*stim_aspec;
gabor_psth_spec = gabor_output_spec.*eye_pos_fft;

freq_range = [-4 4];
alpha_2d = 1 - trapz(trapz(gabor_psth_spec.^2))/trapz(trapz(gabor_output_spec.^2));

close all

%plot spectrum of gabor filter
f1 = figure();
imagesc(fax,fax,gabor_fft);
set(gca,'ydir','normal');
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');
xlabel('Frequency (cyc/deg)');
ylabel('Frequency (cyc/deg)');

%plot stimulus power spec
f2 = figure();
imagesc(fax,fax,log10(stim_aspec.^2));
set(gca,'ydir','normal');
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');
colorbar
caxis([-5 0]);
xlabel('Frequency (cyc/deg)');
ylabel('Frequency (cyc/deg)');

f3 = figure();
imagesc(fax,fax,eye_pos_fft);
set(gca,'ydir','normal');
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');
xlabel('Frequency (cyc/deg)');
ylabel('Frequency (cyc/deg)');

f4 = figure();
imagesc(fax,fax,gabor_output_spec);
set(gca,'ydir','normal');
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');
ca = caxis();
xlabel('Frequency (cyc/deg)');
ylabel('Frequency (cyc/deg)');

f5 = figure();
imagesc(fax,fax,gabor_psth_spec);
set(gca,'ydir','normal');
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');
caxis(ca);
xlabel('Frequency (cyc/deg)');
ylabel('Frequency (cyc/deg)');

f6 = figure();
imagesc(xax,xax,real(gabor_filt));
set(gca,'ydir','normal');
colormap(jet(1e3));
xlim([-1 1]); ylim([-1 1]);
caxis([-1 1]);
xlabel('Position (deg)');
ylabel('Position (deg)');

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'NS2d_gaborfft.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'NS2d_stimfft.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'NS2d_EPfft.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f3);
% 
% figufy(f4);
% fname = [fig_dir 'NS2d_goutfft.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f4);
% 
% figufy(f5);
% fname = [fig_dir 'NS2d_psthfft.pdf'];
% exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f5);

figufy(f6);
fname = [fig_dir 'NS2d_gaborfilt.pdf'];
exportfig(f6,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f6);

%% check dependence of alpha on spectral properties
poss_exps = [-4:0.25:4];
[gabor_alpha_anal] = deal(nan(length(poss_exps),1));

xax = -6:.01:6; %pixel axis (deg)
if mod(length(xax),2) ~= 0
    xax(end) = [];
end
Fs = 1/(pix_width/usfac); %spatial sample freq
N = length(xax); %number of spatial samples
fax = 0:Fs/N:Fs/2;

% make gabor filter
sf = 2; %spatial freq
env_sigma = 1/sf*bandwidth; %spatial envelope SD

%make gabor filter and its fft
gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf + spatial_phase);
gabor_fft = abs(fft(gabor_filt));
gabor_fft = gabor_fft(1:N/2+1);

% compute eye position distribution
ep_dist = exp(-xax.^2/(2*eye_pos_sigma^2)); %EP distribution
ep_dist = ep_dist/sum(ep_dist);
ep_dist_convfilt = ep_dist/sum(ep_dist); %normalize for convolution
ep_dist_fft = abs(fft(ep_dist));
ep_dist_fft = ep_dist_fft(1:N/2+1);

for ii = 1:length(poss_exps)    
    stim_fft = sqrt(1./fax.^(poss_exps(ii)));
    stim_fft(1) = stim_fft(2);
    
    %spectrum of the output of linear gabor filter
    gabor_out_fft = stim_fft.*gabor_fft;
    
    %alpha
    gabor_alpha_anal(ii) = 1 - trapz(fax,(gabor_out_fft.*ep_dist_fft).^2)/trapz(fax,gabor_out_fft.^2);
end

f1 = figure();
plot(-poss_exps,gabor_alpha_anal,'linewidth',2);
xlabel('power law exponent');
ylabel('Alpha');
ylim([0 1]);

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'alpha_vs_exp.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% check dependence of alpha on EP SD
poss_SDs = 0:0.01:0.2;
poss_SFs = [0.5 1 2 4];

rel_pwidth = 1; %pixels per bar
pix_width = 0.01; %bar width in deg
npix = round(10/pix_width); %number of bars
nsamps = 5e3; %number of random pattern samples ('trials')
dds = 67; %dot density
usfac = 1; %spatial up-sampling factor for calculations
npix_us = npix*usfac; %number of pixels in up-sampled stim
Xstim_up = generate_RLS_stim(nsamps,npix,dds,usfac);

%convert to dims [space X trials]
Xstim_up = Xstim_up';

%position axis for stimulus
xax = (1:npix_us)*pix_width/usfac;
xax = xax - mean(xax);
Fs = 1/(pix_width/usfac); %spatial sample freq
N = length(xax); %number of spatial samples
fax = 0:Fs/N:Fs/2; %frequency axis (single sided)

% make gabor filter
[gabor_alpha,simple_alpha,complex_alpha] = deal(nan(length(poss_SDs),length(poss_SDs)));
for jj = 1:length(poss_SFs)
    jj
    sf = poss_SFs(jj); %spatial freq
    env_sigma = 1/sf*bandwidth; %spatial envelope SD
    
    %make gabor filter and its fft
    gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf + spatial_phase);
    gabor_fft = abs(fft(gabor_filt));
    gabor_fft = gabor_fft(1:N/2+1);
    
    %output of gabor filter (complex-valued)
    Xconv = convn(Xstim_up,gabor_filt(:),'same');
    Xconv = Xconv/std(real(Xconv(:))); %normalize to unit SD
    
    %firing rate output of simple-cell model
    simple_rate = spkNL_alpha*log(1+exp((real(Xconv)+spkNL_theta)*spkNL_beta));
    
    %energy model output
    Xout = (real(Xconv).^2 + imag(Xconv).^2);
    complex_rate = spkNL_alpha*log(1+exp((Xout+spkNL_theta)*spkNL_beta));
    
    %subtract out avg rates
    simple_rate = simple_rate - mean(simple_rate(:));
    complex_rate = complex_rate - mean(complex_rate(:));
    
    simple_fft = mean(abs(fft(simple_rate)),2);
    complex_fft = mean(abs(fft(complex_rate )),2);
    simple_fft = simple_fft(1:N/2 + 1);
    complex_fft = complex_fft(1:N/2 + 1);
    simple_tot_var = var(simple_rate(:));
    energy_tot_var = var(complex_rate(:));
    
    for ii = 1:length(poss_SDs)
        % compute eye position distribution
        ep_dist = exp(-xax.^2/(2*poss_SDs(ii)^2)); %EP distribution
        ep_dist = ep_dist/sum(ep_dist);
        ep_dist_convfilt = ep_dist/sum(ep_dist); %normalize for convolution
        ep_dist_fft = abs(fft(ep_dist));
        ep_dist_fft = ep_dist_fft(1:N/2+1);
        
        %     %spectrum of the output of linear gabor filter (assume white noise
        %     gabor_out_fft = gabor_fft;
        
        % compute psths by convolving the model-rates with the EP distribution
        simple_psth = convn(simple_rate,ep_dist_convfilt(:),'same');
        complex_psth = convn(complex_rate,ep_dist_convfilt(:),'same');
        
        % calculate alphas
        simple_psth_var = var(simple_psth(:));
        energy_psth_var = var(complex_psth(:));
        
        %alpha
        %     gabor_alpha(ii) = 1 - trapz(fax,(gabor_out_fft.*ep_dist_fft).^2)/trapz(fax,gabor_out_fft.^2);
        simple_alpha(ii,jj) = 1 - trapz(fax,(simple_fft.*ep_dist_fft').^2)/trapz(fax,simple_fft.^2);
        complex_alpha(ii,jj) = 1 - trapz(fax,(complex_fft.*ep_dist_fft').^2)/trapz(fax,complex_fft.^2);
    end
end
% gabor_alpha(poss_SDs==0) = 0;
simple_alpha(poss_SDs==0,:) = 0;
complex_alpha(poss_SDs==0,:) = 0;

line_widths = [0.5 1 2 3];
f1 = figure(); hold on
for jj = 1:length(poss_SFs)
    % plot(poss_SDs,gabor_alpha,'linewidth',2);
    plot(poss_SDs,simple_alpha(:,jj),'k','linewidth',line_widths(jj));
    plot(poss_SDs,complex_alpha(:,jj),'r','linewidth',line_widths(jj));
end
xlabel('Eye pos SD');
ylabel('Alpha');
ylim([0 1]);

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'alpha_vs_EPSD.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% check dependence of alpha on bar width with white noise expts
% poss_bar_widths = [0.01 0.025 0.05 0.075 0.1 0.2];
% [simple_alpha_anal,complex_alpha_anal,simple_alpha,complex_alpha] = deal(nan(length(poss_bar_widths),1));
% for ii = 1:length(poss_bar_widths)
%     ii
%     %% make spatial stim
%     
%     pix_width = poss_bar_widths(ii); %bar width in deg
%     rel_pwidth = round(pix_width/.01);
%     npix = round(5/pix_width); %number of bars
%     nsamps = 5e3; %number of random pattern samples ('trials')
%     dds = 67; %dot density
%     usfac = rel_pwidth; %spatial up-sampling factor for calculations
%     npix_us = npix*usfac; %number of pixels in up-sampled stim
%     
%     Xstim_up = generate_RLS_stim(nsamps,npix,dds,usfac);
%     
%     %convert to dims [space X trials]
%     Xstim_up = Xstim_up';
%     
%     %position axis for stimulus
%     xax = (1:npix_us)*pix_width/usfac;
%     xax = xax - mean(xax);
%     
%     %% compute eye position distribution
%     ep_dist = exp(-xax.^2/(2*eye_pos_sigma^2)); %EP distribution
%     ep_dist = ep_dist/sum(ep_dist);
%     ep_dist_convfilt = ep_dist/sum(ep_dist); %normalize for convolution
%     ep_dist_fft = abs(fftshift(fft(ep_dist)));
%     niqf = 1/(pix_width/usfac)/2;
%     fax = linspace(-niqf,niqf,length(xax));
%     
%     %% make gabor filter
%     sf = 2; %spatial freq
%     bandwidth = 1/2;
%     env_sigma = 1/sf*bandwidth; %spatial envelope SD
%     spatial_phase = 0;
%     
%     gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf + spatial_phase);
%     
%     %% compute firing rate outputs of these models
%     %output of gabor filter (complex-valued)
%     Xconv = convn(Xstim_up,gabor_filt(:),'same');
%     Xconv = Xconv/std(real(Xconv(:))); %normalize to unit SD
%     
%     %firing rate output of simple-cell model
%     simple_rate = spkNL_alpha*log(1+exp((real(Xconv)+spkNL_theta)*spkNL_beta));
%     
%     %energy model output
%     Xout = (real(Xconv).^2 + imag(Xconv).^2);
%     complex_rate = spkNL_alpha*log(1+exp((Xout+spkNL_theta)*spkNL_beta));
%     
%     %%
%     simple_fft = mean(abs(fftshift(fft(simple_rate - mean(simple_rate(:))))),2);
%     complex_fft = mean(abs(fftshift(fft(complex_rate - mean(complex_rate(:))))),2);
%     
%     %% compute psths by convolving the model-rates with the EP distribution
%     simple_psth = convn(simple_rate,ep_dist_convfilt(:),'same');
%     complex_psth = convn(complex_rate,ep_dist_convfilt(:),'same');
%     
%     %% calculate alphas
%     simple_tot_var = var(simple_rate(:));
%     simple_psth_var = var(simple_psth(:));
%     energy_tot_var = var(complex_rate(:));
%     energy_psth_var = var(complex_psth(:));
%     
%     simple_alpha(ii) = 1 - simple_psth_var/simple_tot_var;
%     complex_alpha(ii) = 1 - energy_psth_var/energy_tot_var;
%     
%     simple_alpha_anal(ii) = 1 - trapz(fax,(simple_fft.*ep_dist_fft').^2)/trapz(fax,simple_fft.^2);
%     %     complex_alpha_anal(ii) = 1 - sum((complex_fft.*ep_dist_fft').^2)/sum(complex_fft.^2);
%     complex_alpha_anal(ii) = 1 - trapz(fax,(complex_fft.*ep_dist_fft').^2)/trapz(fax,complex_fft.^2);
% end

