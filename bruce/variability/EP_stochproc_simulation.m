% close all
clear all

% fig_dir = '/home/james/Analysis/bruce/variability/stoch_proc_sim/';
fig_dir = '/Users/james/Analysis/bruce/variability/stoch_proc_sim/';

% stim_type = 'grating';
stim_type = 'rls';

eye_pos_sigma = 0.11; %SD (deg) of eye position distribution (assume gaussian)

%spk NL parameters of model neuron
spkNL_alpha = 1;
spkNL_beta = 2;
spkNL_theta = 0;

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
    rel_pwidth = 1;
    pix_width = 0.01; %bar width in deg
    npix = round(5/pix_width); %number of bars
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
end

%% compute eye position distribution
ep_dist = exp(-xax.^2/(2*eye_pos_sigma^2)); %EP distribution
ep_dist = ep_dist/sum(ep_dist);
ep_dist_convfilt = ep_dist/sum(ep_dist); %normalize for convolution
ep_dist_fft = abs(fftshift(fft(ep_dist)));
niqf = 1/(pix_width/usfac)/2;
fax = linspace(-niqf,niqf,length(xax));

%% run analysis using a range of different filter SFs
% poss_SFs = [0.25 0.5 1 2 3 4 5 6 7 8];
poss_SFs = 0.25:0.25:8;
[simple_alpha,complex_alpha] = deal(nan(length(poss_SFs),1));
for ii = 1:length(poss_SFs)
    
    %% make gabor filter
    fprintf('SF %d/%d\n',ii,length(poss_SFs));
    sf = poss_SFs(ii); %spatial freq
    bandwidth = 1/2;
    env_sigma = 1/sf*bandwidth; %spatial envelope SD
    spatial_phase = 0;
    
    gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf + spatial_phase);
    
    %% compute firing rate outputs of these models
    %output of gabor filter (complex-valued)
    Xconv = convn(Xstim_up,gabor_filt(:),'same');
    Xconv = Xconv/std(real(Xconv(:))); %normalize to unit SD
    
    %firing rate output of simple-cell model
    simple_rate = spkNL_alpha*log(1+exp((real(Xconv)+spkNL_theta)*spkNL_beta));
    
    %energy model output
    Xout = (real(Xconv).^2 + imag(Xconv).^2);
    complex_rate = spkNL_alpha*log(1+exp((Xout+spkNL_theta)*spkNL_beta));
    
    %%
    simple_fft = mean(abs(fftshift(fft(simple_rate - mean(simple_rate(:))))),2);
    complex_fft = mean(abs(fftshift(fft(complex_rate - mean(complex_rate(:))))),2);
    
    %% compute psths by convolving the model-rates with the EP distribution
    simple_psth = convn(simple_rate,ep_dist_convfilt(:),'same');
    complex_psth = convn(complex_rate,ep_dist_convfilt(:),'same');
    
    %% calculate alphas
    simple_tot_var = var(simple_rate(:));
    simple_psth_var = var(simple_psth(:));
    energy_tot_var = var(complex_rate(:));
    energy_psth_var = var(complex_psth(:));
    
    simple_alpha(ii) = 1 - simple_psth_var/simple_tot_var;
    complex_alpha(ii) = 1 - energy_psth_var/energy_tot_var; 
    
    simple_alpha_anal(ii) = 1 - trapz(fax,(simple_fft.*ep_dist_fft').^2)/trapz(fax,simple_fft.^2);
%     complex_alpha_anal(ii) = 1 - sum((complex_fft.*ep_dist_fft').^2)/sum(complex_fft.^2);
    complex_alpha_anal(ii) = 1 - trapz(fax,(complex_fft.*ep_dist_fft').^2)/trapz(fax,complex_fft.^2);
end

% plot filter SF vs alpha for simple and complex models
close all
f1 = figure(); hold on
% plot(poss_SFs,simple_alpha,'o-');
% plot(poss_SFs,complex_alpha,'ro-');
plot(poss_SFs,simple_alpha_anal,'-');
plot(poss_SFs,complex_alpha_anal,'r-');
xlabel('Preferred spatial frequency (cyc/deg)');
ylabel('Alpha');



%% now run analysis using an example gabor filter
sf = 4; %spatial freq
bandwidth = 1/2;
env_sigma = 1/sf*bandwidth; %spatial envelope SD
spatial_phase = 0;
gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf + spatial_phase);

% get amplitude spectrum of the gabor filter (analytically)
fft_bandwidth = 1/(2*pi*env_sigma); %frequency-domain width (SD) of gabor
gabor_fft = exp(-(fax-sf).^2/(2*fft_bandwidth^2)); %gaussian centered on SF

%% make spatial stim
if strcmp(stim_type,'rls')
    rel_pwidth = 1;
    pix_width = 0.01; %bar width in deg
    npix = round(5/pix_width); %number of bars
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
end
stim_fft = mean(abs(fftshift(fft(Xstim_up))),2);

%% compute firing rate outputs of these models
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
gout_fft = mean(abs(fftshift(fft(real(Xconv)))),2);

%simple cell firing rate output
simprate_fft = mean(abs(fftshift(fft(simple_rate - mean(simple_rate(:))))),2);

%energy model firing rate output
comprate_fft = mean(abs(fftshift(fft(complex_rate - mean(complex_rate(:))))),2);

%simple cell PSTH
simppsth_fft = mean(abs(fftshift(fft(simple_psth - mean(simple_psth(:))))),2);

%energy model PSTH
comppsth_fft = mean(abs(fftshift(fft(complex_psth - mean(complex_psth(:))))),2);
   
%% plot amplitude spectra
freq_range = [0 10];

f1 = figure();
subplot(2,1,1);
hold on
plot(fax,stim_fft/max(stim_fft),'b');
plot(fax,ep_dist_fft/max(ep_dist_fft),'k');
plot(fax,gabor_fft/max(gabor_fft),'r')
plot(fax,simprate_fft/max(simprate_fft),'m')
plot(fax,simppsth_fft/max(simprate_fft),'m--');
xlim(freq_range);
xlabel('Frequency (cyc/deg)');
ylabel('Relative amplitude');

subplot(2,1,2);
plot(fax,stim_fft/max(stim_fft),'b');
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
% fname = [fig_dir sprintf('sim_pspecs_%s.pdf',stim_type)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


%% check dependence of alpha on bar width with white noise expts
poss_bar_widths = [0.01 0.025 0.05 0.075 0.1 0.2];
[simple_alpha_anal,complex_alpha_anal,simple_alpha,complex_alpha] = deal(nan(length(poss_bar_widths),1));
for ii = 1:length(poss_bar_widths)
    ii
    %% make spatial stim
    
    pix_width = poss_bar_widths(ii); %bar width in deg
    rel_pwidth = round(pix_width/.01);
    npix = round(5/pix_width); %number of bars
    nsamps = 5e3; %number of random pattern samples ('trials')
    dds = 67; %dot density
    usfac = rel_pwidth; %spatial up-sampling factor for calculations
    npix_us = npix*usfac; %number of pixels in up-sampled stim
    
    Xstim_up = generate_RLS_stim(nsamps,npix,dds,usfac);
    
    %convert to dims [space X trials]
    Xstim_up = Xstim_up';
    
    %position axis for stimulus
    xax = (1:npix_us)*pix_width/usfac;
    xax = xax - mean(xax);
    
    %% compute eye position distribution
    ep_dist = exp(-xax.^2/(2*eye_pos_sigma^2)); %EP distribution
    ep_dist = ep_dist/sum(ep_dist);
    ep_dist_convfilt = ep_dist/sum(ep_dist); %normalize for convolution
    ep_dist_fft = abs(fftshift(fft(ep_dist)));
    niqf = 1/(pix_width/usfac)/2;
    fax = linspace(-niqf,niqf,length(xax));
    
    %% make gabor filter
    sf = 2; %spatial freq
    bandwidth = 1/2;
    env_sigma = 1/sf*bandwidth; %spatial envelope SD
    spatial_phase = 0;
    
    gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf + spatial_phase);
    
    %% compute firing rate outputs of these models
    %output of gabor filter (complex-valued)
    Xconv = convn(Xstim_up,gabor_filt(:),'same');
    Xconv = Xconv/std(real(Xconv(:))); %normalize to unit SD
    
    %firing rate output of simple-cell model
    simple_rate = spkNL_alpha*log(1+exp((real(Xconv)+spkNL_theta)*spkNL_beta));
    
    %energy model output
    Xout = (real(Xconv).^2 + imag(Xconv).^2);
    complex_rate = spkNL_alpha*log(1+exp((Xout+spkNL_theta)*spkNL_beta));
    
    %%
    simple_fft = mean(abs(fftshift(fft(simple_rate - mean(simple_rate(:))))),2);
    complex_fft = mean(abs(fftshift(fft(complex_rate - mean(complex_rate(:))))),2);
    
    %% compute psths by convolving the model-rates with the EP distribution
    simple_psth = convn(simple_rate,ep_dist_convfilt(:),'same');
    complex_psth = convn(complex_rate,ep_dist_convfilt(:),'same');
    
    %% calculate alphas
    simple_tot_var = var(simple_rate(:));
    simple_psth_var = var(simple_psth(:));
    energy_tot_var = var(complex_rate(:));
    energy_psth_var = var(complex_psth(:));
    
    simple_alpha(ii) = 1 - simple_psth_var/simple_tot_var;
    complex_alpha(ii) = 1 - energy_psth_var/energy_tot_var;
    
    simple_alpha_anal(ii) = 1 - trapz(fax,(simple_fft.*ep_dist_fft').^2)/trapz(fax,simple_fft.^2);
    %     complex_alpha_anal(ii) = 1 - sum((complex_fft.*ep_dist_fft').^2)/sum(complex_fft.^2);
    complex_alpha_anal(ii) = 1 - trapz(fax,(complex_fft.*ep_dist_fft').^2)/trapz(fax,complex_fft.^2);
end

f1 = figure(); hold on
plot(poss_bar_widths,simple_alpha_anal,'o-');
plot(poss_bar_widths,complex_alpha_anal,'ro-');
% plot(poss_bar_widths,simple_alpha,'o-');
% plot(poss_bar_widths,complex_alpha,'go-');



%% eye position distribution
f3 = figure();
plot(ep_ax,ep_dist,'k');
xlim([-0.5 0.5]);
xlabel('Eye position');
ylabel('Probability');

% fig_width = 4; rel_height = 0.8;
% figufy(f3);
% fname = [fig_dir 'EP_dist.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f3);

%% plot gabor filter
f4 = figure();
plot(xax,real(gabor_filt),'k');
hold on
plot(xax,imag(gabor_filt),'r');
xlim([-0.5 0.5]);
xlabel('Eye position');
ylabel('Filter amplitude');
xl = xlim(); line(xl,[0 0],'color','k','linestyle','--');

% fig_width = 4; rel_height = 0.8;
% figufy(f4);
% fname = [fig_dir 'Gabor_filt.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f4);

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
    
    %     fig_width = 4; rel_height = 0.8;
    %     figufy(f5);
    %     fname = [fig_dir 'Grating_stim.pdf'];
    %     exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    %     % close(f5);
    
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

% fig_width = 4; rel_height = 0.8;
% figufy(f6);
% fname = [fig_dir 'Grating_alphas.pdf'];
% exportfig(f6,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f6);


%%
stim_type = 'white';
% stim_type = 'NS';

sf = 2; %spatial freq
bandwidth = 0.5;
env_sigma = 1/sf*bandwidth; %spatial envelope SD
spatial_phase = 0;
spatial_AR = 0.75;
orientation = 45;
  
pix_dx = .01;
xax = -4:pix_dx:4;
[X,Y] = meshgrid(xax);
Xp = X*cosd(orientation) + Y*sind(orientation);
Yp = -X*sind(orientation) + Y*cosd(orientation);

gabor_filt = exp(-(Xp.^2 + spatial_AR*Yp.^2)/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*Xp*sf + spatial_phase);
gabor_filt = gabor_filt - mean(gabor_filt(:));

gabor_fft = abs(fftshift(fft2(real(gabor_filt))));
niqf = 1/pix_dx/2;
fax = linspace(-niqf,niqf,length(xax));

freq_range = [-5 5];

f1 = figure();
imagesc(fax,fax,gabor_fft);
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');

[Fx,Fy] = meshgrid(fax);
FF = sqrt(Fx.^2 + Fy.^2);
if strcmp(stim_type,'NS')
stim_aspec = 1./(FF.^2); 
stim_aspec(FF==0) = nan; stim_aspec(FF==0) = nanmax(stim_aspec(:));
stim_aspec = stim_aspec/max(stim_aspec(:));
f2 = figure();
imagesc(fax,fax,log10(stim_aspec));
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');
elseif strcmp(stim_type,'white');
    stim_aspec = ones(size(FF));
% f2 = figure();
% imagesc(fax,fax,log10(stim_aspec));
% xlim(freq_range); ylim(freq_range);
% line(freq_range,[0 0],'color','w')
% line([0 0],freq_range,'color','w');
end  

f3 = figure();
eye_pos_sigma_2d = sqrt(2*eye_pos_sigma.^2);
eye_pos_fft_sigma = 1/(2*pi*eye_pos_sigma_2d);
eye_pos_fft = exp(-(Fx.^2 + Fy.^2)/(2*eye_pos_fft_sigma^2));
imagescnan(fax,fax,eye_pos_fft);
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');

gabor_output_spec = gabor_fft.*stim_aspec;
f4 = figure();
imagesc(fax,fax,gabor_output_spec);
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');
ca = caxis();

gabor_psth_spec = gabor_output_spec.*eye_pos_fft;
f5 = figure();
imagesc(fax,fax,gabor_psth_spec);
xlim(freq_range); ylim(freq_range);
line(freq_range,[0 0],'color','w')
line([0 0],freq_range,'color','w');
caxis(ca);

alpha_2d = 1 - trapz(trapz(gabor_psth_spec.^2))/trapz(trapz(gabor_output_spec.^2))