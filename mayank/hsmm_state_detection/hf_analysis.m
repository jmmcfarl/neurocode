clear all
close all

load G:\WC_Germany\persistent_2010\pers_revised_dir_2010
addpath('G:\Code\informationTheoryToolboxv1.0')
Fs=  2016;
dsf = 8;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[0.1/niqf 100/niqf]);

min_freq = 3; max_freq = 100; delta_j = 0.075;
k0 = 6; %wavenumber for morlet wavelet
fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));
min_scale = 1/max_freq/fourier_factor;
max_scale = 1/min_freq/fourier_factor;
n_scales = round(log2(max_scale/min_scale)/delta_j);

%%
for d = 1:29
    
    cd(dir_array{d})
    disp(dir_array{d})
    
    load used_data wcv_minus_spike lf8
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = zscore(downsample(wcv_f,dsf));
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = zscore(downsample(lf8_f,dsf));
    
    [wav_trans,periods,scales,COI] = wavelet(wcv_f,1/Fsd,0,delta_j,min_scale,n_scales);
    wfreqs = 1./periods;
    inv_scales = 1./scales';
    raw_scalogram = abs(wav_trans).^2;
    n = size(raw_scalogram,2);
%     sm_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*raw_scalogram,1,delta_j,scales);
        
    [wav_trans8,periods,scales,COI] = wavelet(lf8_f,1/Fsd,0,delta_j,min_scale,n_scales);
    raw_scalogram8 = abs(wav_trans8).^2;
%     sm_scalogram8 = smoothwavelet(inv_scales(:,ones(1,n)).*raw_scalogram8,1,delta_j,scales);

raw_scalogram = zscore(raw_scalogram')';
raw_scalogram8 = zscore(raw_scalogram8')';
split = 3:60;
for i = 1:length(split)
    vec1 = mean(raw_scalogram(1:split(i),:));
    vec2 = mean(raw_scalogram(split(i)+1:end,:));
    temp = corrcoef(vec1,vec2);
    c(i) = temp(2,1);
    i
end

    %transform to Nxp data matrix
%     raw_scalogram = log(raw_scalogram');
%     raw_scalogram8 = log(raw_scalogram8');
    raw_scalogram = raw_scalogram';
    raw_scalogram8 = raw_scalogram8';
%     raw_scalogram = zscore(raw_scalogram);
%     raw_scalogram = log(raw_scalogram - repmat(min(raw_scalogram),length(wcv_f),1) + 1);
%     raw_scalogram8 = zscore(raw_scalogram8);
%     raw_scalogram8 = log(raw_scalogram8 - repmat(min(raw_scalogram8),length(lf8_f),1) + 1);
    
%     sm_scalogram = sm_scalogram';
%      sm_scalogram8 = sm_scalogram8';
   
    c = corrcoef([raw_scalogram raw_scalogram8]);
%     c = corrcoef(zscore([raw_scalogram raw_scalogram8]));
    c1 = c(1:length(wfreqs),1:length(wfreqs));
    c2 = c(1:length(wfreqs),length(wfreqs)+1:end);
    c3 = c(length(wfreqs)+1:end,length(wfreqs)+1:end);
    figure
    subplot(1,3,1)
    pcolor(wfreqs,wfreqs,c1);shading flat, caxis([0.1 0.7])
    set(gca,'xscale','log','yscale','log')
    subplot(1,3,2)
    pcolor(wfreqs,wfreqs,c2);shading flat, caxis([0.1 0.4])
    set(gca,'xscale','log','yscale','log')
    subplot(1,3,3)
    pcolor(wfreqs,wfreqs,c3);shading flat, caxis([0.1 0.6])
    set(gca,'xscale','log','yscale','log')
    
    [COEFF, SCORE, LATENT] = princomp(raw_scalogram);
    [COEFF8, SCORE8, LATENT8] = princomp(raw_scalogram8);
%      [zCOEFF, zSCORE, zLATENT] = princomp(zscore(sm_scalogram));
%     [zCOEFF8, zSCORE8, zLATENT8] = princomp(zscore(sm_scalogram8));
%    zSCORE = zscore(zSCORE);
%    zSCORE8 = zscore(zSCORE8);
   
  figure
  subplot(2,1,1)
    plot(wfreqs,COEFF(:,1))
hold on
plot(wfreqs,COEFF(:,2),'r')
plot(wfreqs,COEFF(:,3),'k')
set(gca,'xscale','log'), axis tight
    subplot(2,1,2)
    plot(wfreqs,COEFF8(:,1))
hold on
plot(wfreqs,COEFF8(:,2),'r')
plot(wfreqs,COEFF8(:,3),'k')
set(gca,'xscale','log'), axis tight

end



