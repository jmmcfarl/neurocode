clear all
close all

load ./random_bar_phase_models_lfp_ds
% load ./random_bar_phase_models_lfp_nsm3.mat
% load ./random_bar_phase_models_lfp_nsm4.mat
% load ./random_bar_phase_models_csd_nsm3.mat
%%

% %for full spatial sampling
% [WF,CH] = ndgrid(wfreqs,(1:24)*0.05);
% linch = (1:.4:24)*0.05; %new depth axis for interpolation
% n_chs = 24;
% 
%for 100um spatial sampling
[WF,CH] = ndgrid(wfreqs,(1:12)*0.1);
linch = (1:.4:12)*0.1;
n_chs = 12;

WF = WF'; CH = CH';

%create new frequency and space axes for interpolation
linw = linspace(wfreqs(end),wfreqs(1),250);
% linw = logspace(log10(wfreqs(end)),log10(wfreqs(1)),250);
[lWF,lCH] = ndgrid(linw,linch);

use_freqs = 2:(length(linw)-1); %don't plot lowest and highest freqs after interpolation to help minimize edge artifacts of interpolation
sat_factor = 0.9; %saturate amplitude kernel at this fraction of the maximum

for cc = 1:24;
    fprintf('Unit %d of %d\n',cc,24);
    
    cur_ampkern = post_ampkern(cc,:); %with stimulus
%     cur_ampkern = [reshape(cur_ampkern',length(wfreqs),n_chs)'; zeros(1,length(wfreqs))]; %pad with a zero row for pcolor plotting
    cur_ampkern = reshape(cur_ampkern',length(wfreqs),n_chs)'; 
    cur_ampkern = interp2(WF,CH,cur_ampkern,lWF,lCH); %interpolate amplitude kernel 
    cur_ampkern = cur_ampkern(use_freqs,:); %throw out lowest and highest frequencies
    cur_ampkern = cur_ampkern/max(cur_ampkern(:))/sat_factor;
    cur_ampkern(cur_ampkern > 1) = 1; %set max value to 1
    
    cur_phasekern = post_phasekern(cc,:); %with stimulus
%     cur_phasekern = [reshape(cur_phasekern',length(wfreqs),n_chs)'; zeros(1,length(wfreqs))]; %pad with a zero row for pcolor plotting
    cur_phasekern = reshape(cur_phasekern',length(wfreqs),n_chs)'; 
    cur_phasekern = cur_phasekern/360;
    cur_phasekern = interp2(WF,CH,cur_phasekern,lWF,lCH);
    cur_phasekern = cur_phasekern(use_freqs,:);
    
    C = ones([size(cur_ampkern') 3]); %initialize color matrix in HSV
    C(:,:,1) = cur_phasekern'; %set hue values to phase kernel
    C(:,:,3) = cur_ampkern'; %set value to amplitude kernel    
    C = hsv2rgb(C); %convert to RGB
    
    subplot(2,1,1)
    imagesc(linw(2:end-2),linch,C)
    % set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Depth (um)','fontsize',16)
    colorbar
    subplot(2,1,2)
    imagesc(linw(2:end-2),linch,(cur_ampkern'));shading interp
    % set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Depth (um)','fontsize',16)
    colorbar
    
    pause
    clf
end


%%
good_sus = [3 4 5 14 15 17];
stim_imp = (xv_null_LL - xv_so_LL)/log(2);
sac_imp = (xv_null_LL - xv_sac_LL)/log(2);
mulsac_imp = (xv_null_LL - xv_mulsac_LL)/log(2);
phase_imp = (xv_null_LL - xv_po_LL)/log(2);
stimphase_imp = (xv_null_LL - xv_post_LL)/log(2);

figure
boxplot([stim_imp(:) sac_imp(:) mulsac_imp(:) phase_imp(:) stimphase_imp(:)]);

sus = good_sus;
mus = setdiff(1:24,sus);
figure
subplot(1,2,1)
boxplot([stim_imp(sus)' sac_imp(sus)' mulsac_imp(sus)' phase_imp(sus)' stimphase_imp(sus)']);
ylim([0 1.3])
subplot(1,2,2)
boxplot([stim_imp(mus)' sac_imp(mus)' mulsac_imp(mus)' phase_imp(mus)' stimphase_imp(mus)']);
ylim([0 1.3])
