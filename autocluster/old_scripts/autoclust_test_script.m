clear all
close all
addpath('~/James_scripts/autocluster/');

%%
Expt_name = 'G086';
cd('~/James_scripts/autocluster/');

block_num = 12;
probe_num = 42;

%% LOAD VOLTAGE SIGNAL
sfile_name = sprintf('Expt%d.p%dFullV.mat',block_num,probe_num);

add_Vmean = 0; %this determines whether or not to include the across-electrode average voltage

filt_cutoff = [100 nan]; %[lcf hcf] in Hz. use nan when you dont want a cutoff.

%loads in (high-pass filtered) voltage signal
[V,Vtime,Fs] = Load_FullV(sfile_name, add_Vmean, filt_cutoff);

%%
thresh_sign = -1; %detect of negative or positive voltage peaks
target_rate = 100; %in Hz
target_Nspks = target_rate*range(Vtime); 

spk_pts = [-12:27]; %range of samples relative to each detected peak to use

[spk_id, sig_th] = triggerSpikes(V,thresh_sign,target_Nspks);
spk_id(spk_id <= abs(spk_pts(1)) | spk_id >= length(V)-spk_pts(end)) = []; %get rid of spikes at the edges

%extract spike snippets
Spikes = getSpikeSnippets(V,Vtime,spk_id,spk_pts);
[N_spks, N_samps, N_chs] = size(Spikes.V);

%% 2-COMPONENT GMM AUTOCLUSTERING USING MULTIPLE FEATURES AND MULTIPLE INITIALIZATIONS
params.outlier_thresh =7; %in Mah. distance units
params.verbose = 2;
params.use_best_only = 0; %if you want to use only the best spike waveform (and time-derivative) as a template
params.cluster_bias = 0.85; %if you want to control for type1/type2 errors

[su_inds, gmm_fit, cluster_details, spike_xy, spike_features] = autocluster_2comp(Spikes.V,params);
mu_inds = setdiff(1:N_spks,su_inds);
cluster_labels = cluster_details.cluster_labels;
cluster_idx = cluster_details.cluster_idx;

%get ISIs
su_spk_times = Spikes.times(su_inds);
isis = diff(su_spk_times)*1e3; %in ms

%two different measures of refractoriness
refractoriness(1) = sum(isis < 1)/length(isis)*100;
refractoriness(2) = sum(isis < 2)/length(isis)*100;
cluster_details.refract = refractoriness;

%compute two more measures of cluster quality
[Lratio,iso_dist] = compute_cluster_Lratio(spike_features,gmm_fit,cluster_idx,cluster_labels);

%% CREATE SUMMARY FIGURE
fprintf('Creating summary figure\n');
figure;

% DENSITY PLOT IN XY SPACE
subplot(2,3,1);
[handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
set(gca,'ytick',[]);
if cluster_details.used_it == -1
    title('Used PCs','fontsize',14);
elseif cluster_details.used_it == -2
    title('Used voltage proj','fontsize',14);
    
else
    title('Used template projection','fontsize',14);
end


% SCATTERPLOT IN XY SPACE
subplot(2,3,2); hold on
plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
plot(spike_xy(su_inds,1),spike_xy(su_inds,2),'r.');
for ii = 1:length(cluster_details.cluster_labels)
    if cluster_details.cluster_labels(ii) == 1
        h1 = plot_gaussian_2d(cluster_details.gmm_xyMeans(ii,:)',squeeze(cluster_details.gmm_xySigma(:,:,ii)),[1 2],'g',2);
    else
        h1 = plot_gaussian_2d(cluster_details.gmm_xyMeans(ii,:)',squeeze(cluster_details.gmm_xySigma(:,:,ii)),[1 2],'b',2);
    end
end
set(gca,'ytick',[]);
title(sprintf('D: %.3g  LL: %.3g  L: %.3g  I: %.3g',cluster_details.dprime,cluster_details.LL,Lratio,iso_dist));

% AVG WAVEFORMS
subplot(2,3,3); hold on
plot((1:N_samps)/Fs*1e3,Spikes.V(mu_inds,:),'k','linewidth',0.01);
plot((1:N_samps)/Fs*1e3,Spikes.V(su_inds,:),'r','linewidth',0.01);
errorbar((1:N_samps)/Fs*1e3,cluster_details.mean_spike(:,1),cluster_details.std_spike(:,1),'g','linewidth',1.5)
errorbar((1:N_samps)/Fs*1e3,cluster_details.mean_spike(:,2),cluster_details.std_spike(:,2),'b','linewidth',1.5)
xlabel('Time (ms)','fontsize',14);
ylabel('Amplitude (mV)','fontsize',14);
xlim([1 N_samps]/Fs*1e3);

% ISI DIST
subplot(2,3,4); hold on
isi_bins = [logspace(log10(5e-4),log10(1),100)]*1e3;
isi_hist = histc(isis,isi_bins);
isi_hist = isi_hist/sum(isi_hist);
stairs(isi_bins,isi_hist)
set(gca,'xscale','log');
yl = ylim();
line([2 2],yl,'color','r','linestyle','--');
line([1 1],yl,'color','r');
xlim([0.5 1000]);
xlabel('ISI (ms)','fontsize',14);
ylabel('Relative frequency','fontsize',14);
title(sprintf('%.2f - %.2f refractoriness',refractoriness(1),refractoriness(2)),'fontsize',20);

% Time v X plot
subplot(2,3,5); hold on
plot(Spikes.times(mu_inds),spike_xy(mu_inds,1),'k.')
plot(Spikes.times(su_inds),spike_xy(su_inds,1),'r.')
hold on
plot(Spikes.times(mu_inds),smooth(spike_xy(mu_inds,1),30,'rlowess'),'g','linewidth',2)
plot(Spikes.times(su_inds),smooth(spike_xy(su_inds,1),30,'rlowess'),'b','linewidth',2)
xlabel('Spike time (s)','fontsize',14)
ylabel('Feature projection','fontsize',14);
xlim(Vtime([1 end]));

% Trigger amp dist
trig_ax = linspace(min(Spikes.trig_vals)*1e3,max(Spikes.trig_vals)*1e3,200);
su_trig_hist = histc(Spikes.trig_vals(su_inds)*1e3,trig_ax)/length(su_inds);
mu_trig_hist = histc(Spikes.trig_vals(mu_inds)*1e3,trig_ax)/length(mu_inds);
subplot(2,3,6) ;hold on
stairs(trig_ax,su_trig_hist,'r','linewidth',1);
stairs(trig_ax,mu_trig_hist,'k','linewidth',1);
yl = ylim();
line([sig_th sig_th]*thresh_sign*1e3,yl,'color','b','linewidth',2);
xlabel('Trigger voltage (mV)','fontsize',14);
ylabel('Relative frequency','fontsize',14);



