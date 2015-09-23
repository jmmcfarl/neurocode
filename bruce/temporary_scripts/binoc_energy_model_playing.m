clear all

%spk NL parameters of model neuron
spkNL_alpha = 1;
% spkNL_beta = 1;
% spkNL_theta = 0;

%gabor params
bandwidth = 1/2; %spatial freq bandwidth
spatial_phase = 0; %spatial phase
sf = 2;
env_sigma = 1/sf*bandwidth; %spatial envelope SD

%% make RLS stim
pix_width = 0.1; %bar width in deg (using closer approx of white noise)
npix = round(10/pix_width); %number of bars
nsamps = 1e4; %number of random pattern samples ('trials')
dds = 100; %dot density

usfac = 6; %spatial up-sampling factor for calculations
npix_us = npix*usfac; %number of pixels in up-sampled stim
pix_dx = pix_width/usfac;

%position axis for stimulus
xax = (1:npix_us)*pix_width/usfac;
xax = xax - mean(xax); %center at 0
Fs = 1/(pix_width/usfac); %spatial sample freq
N = length(xax); %number of spatial samples
fax = 0:Fs/N:Fs/2; %frequency axis (single sided)

%% make gabor filters
preferred_disp = 0.1;
gabors_left = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf + spatial_phase);
gabors_right = exp(-(xax - preferred_disp).^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*(xax - preferred_disp)*sf + spatial_phase);

%%
Xstim_left = generate_RLS_stim(nsamps,npix,dds,usfac);
real_outs = Xstim_left*real(gabors_left') + Xstim_left*real(gabors_right)';
% real_outs = abs(real_outs);
real_outs = (real_outs).^2;


real_sup = Xstim_left*real(gabors_left') - Xstim_left*real(gabors_right)';
real_sup = (real_sup).^2;

% imag_outs = Xstim_left*imag(gabors_left') + Xstim_left*imag(gabors_right)';
% imag_outs = abs(imag_outs);
% imag_outs = zeros(size(real_outs));
% imag_outs = (imag_outs).^2;
% G = (imag_outs + real_outs );
G = (real_outs - real_sup);
spkNL_beta = 1/std(G);
spkNL_theta = -mean(G);
%% run analysis using a range of different filter SFs
disp_axis = [-0.8:pix_dx*2:0.8];
cor_disp_tune = nan(size(disp_axis));
acor_disp_tune = nan(size(disp_axis));
for pp = 1:length(disp_axis)
    pp
Xstim_left = generate_RLS_stim(nsamps,npix,dds,usfac);
    cur_disp = disp_axis(pp);
    Xstim_right = shift_matrix_Nd(Xstim_left,round(cur_disp/pix_dx),2);
    
    real_outs = Xstim_left*real(gabors_left') + Xstim_right*real(gabors_right)';
    %     real_outs = abs(real_outs);
    real_outs = (real_outs).^2;
    
    real_sup = Xstim_left*real(gabors_left') - Xstim_right*real(gabors_right)';
    real_sup(abs(real_sup) < 5) = 0;
    real_sup = (real_sup).^2;

    
% imag_outs = zeros(size(real_outs));
%     imag_outs = Xstim_left*imag(gabors_left') + Xstim_right*imag(gabors_right)';
%     imag_outs = abs(imag_outs);
%      imag_outs = (imag_outs).^2;
   
%     G = (imag_outs + real_outs + spkNL_theta)*spkNL_beta;
    G = (real_outs - 5*real_sup + spkNL_theta)*spkNL_beta;
    r = spkNL_alpha*log(1+exp(G));
%     r(G > 50) = spkNL_alpha*G(G > 50);
    cor_disp_tune(pp) = mean(r);


    Xstim_right = -Xstim_right;
    
    real_outs = Xstim_left*real(gabors_left') + Xstim_right*real(gabors_right)';
    real_outs = (real_outs).^2;
    
    real_sup = Xstim_left*real(gabors_left') - Xstim_right*real(gabors_right)';
    real_sup(abs(real_sup) < 5) = 0;
    real_sup = (real_sup).^2;

    G = (real_outs - 5*real_sup + spkNL_theta)*spkNL_beta;
    r = spkNL_alpha*log(1+exp(G));
%     r(G > 50) = spkNL_alpha*G(G > 50);
    acor_disp_tune(pp) = mean(r);

end

figure
plot(disp_axis,cor_disp_tune,'.-');
hold on
plot(disp_axis,acor_disp_tune,'r.-');
