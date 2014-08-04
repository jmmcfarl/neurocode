function sum_fig = create_summary_cluster_fig(clusterDetails,Spikes,spike_xy,params)


%% CREATE SUMMARY FIGURE
if params.summary_plot == 1
    sum_fig = figure('visible','off');
elseif params.summary_plot == 2
    sum_fig = figure;
end
fprintf('Creating summary figure\n');

[N_spks,N_samps,N_chs] = size(Spikes.V);
Fs = 3e4;

% spike_labels = zeros(size(clusterDetails.comp_idx));
% uids = find(clusterDetails.comp_idx > 0);
% spike_labels(uids) = clusterDetails.cluster_labels(clusterDetails.comp_idx(uids));
spike_labels = clusterDetails.spike_clusts;
uids = find(clusterDetails.comp_idx > 0);
mu_inds = find(spike_labels == 1);
out_inds = find(clusterDetails.comp_idx == -1);
N_sus = length(unique(spike_labels(uids))) - 1;
for ii = 1:N_sus
    su_inds{ii} = find(spike_labels == ii + 1);
end
cmap = cluster_cmap(N_sus);

% DENSITY PLOT IN XY SPACE
subplot(2,3,1);
[handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
set(gca,'ytick',[]);
if clusterDetails.template_it == -1
    title('Used PCs','fontsize',12);
elseif clusterDetails.template_it == -2
    title('Used voltage proj','fontsize',12);
elseif clusterDetails.template_it == -3
    title('Used energy proj','fontsize',12);
else
    title('Used template projection','fontsize',12);
end
xrange = minmax(details.x);
yrange = minmax(details.y);

% SCATTERPLOT IN XY SPACE
subplot(2,3,2); hold on
plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.','markersize',4);
for ii = 1:N_sus
    plot(spike_xy(spike_labels==ii+1,1),spike_xy(spike_labels==ii+1,2),'.','color',cmap(ii,:),'markersize',4);
end
plot(spike_xy(out_inds,1),spike_xy(out_inds,2),'r*','markersize',2);
for ii = 1:length(clusterDetails.cluster_labels)
    h1 = plot_gaussian_2d(clusterDetails.gmm_xyMeans(ii,:)',squeeze(clusterDetails.gmm_xySigma(:,:,ii)),[2],'r',2);
end
xlim(xrange); ylim(yrange);
set(gca,'ytick',[]);
%note this only displays SU stats for best SU cluster at this point
title(sprintf('D: %.3g  LL: %.3g  L: %.3g  I: %.3g',clusterDetails.dprime,clusterDetails.LL,clusterDetails.Lratios(1),clusterDetails.iso_dists(1)));

if N_chs > 1
    mean_spike = reshape(clusterDetails.mean_spike,[N_samps N_chs size(clusterDetails.mean_spike,2)]);
    std_spike = reshape(clusterDetails.std_spike,[N_samps N_chs size(clusterDetails.mean_spike,2)]);
    peak_amps = squeeze(max(abs(mean_spike)));
    [best_ch,best_clust] = ind2sub([N_chs size(clusterDetails.mean_spike,2)],find(peak_amps == max(peak_amps(:))));
    
    % AVG WAVEFORMS
    subplot(2,3,3); hold on
    %     n_show = 200;
    %     rand_mu_set = randperm(length(mu_inds));
    %     if length(mu_inds) > n_show; rand_mu_set = rand_mu_set(1:n_show); end;
    %     rand_mu_set = mu_inds(rand_mu_set);
    %     plot((1:N_samps)/Fs*1e3,squeeze(Spikes.V(rand_mu_set,:,best_ch)),'k','linewidth',0.01);
    %     for ii = 1:N_sus
    %         su_inds{ii} = find(spike_labels == ii + 1);
    %         rand_su_set = randperm(length(su_inds{ii}));
    %         if length(su_inds{ii}) > n_show; rand_su_set = rand_su_set(1:n_show); end;
    %         rand_su_set = su_inds{ii}(rand_su_set);
    %         plot((1:N_samps)/Fs*1e3,squeeze(Spikes.V(rand_su_set,:,best_ch)),'color',cmap(ii,:),'linewidth',0.01);
    %     end
    %     errorbar((1:N_samps)/Fs*1e3,squeeze(mean_spike(:,best_ch,1)),squeeze(std_spike(:,best_ch,1)),'g','linewidth',1.5)
    %     for ii = 1:N_sus
    %         errorbar((1:N_samps)/Fs*1e3,squeeze(mean_spike(:,best_ch,ii+1)),squeeze(std_spike(:,best_ch,ii+1)),'b','linewidth',1.5)
    %     end
    %     xlabel('Time (ms)','fontsize',14);
    %     ylabel('Amplitude (mV)','fontsize',14);
    %     xlim([1 N_samps]/Fs*1e3);
    
    cur_offset = 0;
    offset_spacing = -3*max(max(squeeze(std(mean_spike))));
    for cc = 1:N_chs
        errorbar((1:N_samps)/Fs*1e3,squeeze(mean_spike(:,cc,1))+cur_offset,squeeze(std_spike(:,cc,1)),'k','linewidth',1.5)
        hold on
        for ii = 1:N_sus
            errorbar((1:N_samps)/Fs*1e3,squeeze(mean_spike(:,cc,ii+1))+cur_offset,squeeze(std_spike(:,cc,ii+1)),'color',cmap(ii,:),'linewidth',1.5)
        end
        cur_offset = cur_offset + offset_spacing;
    end
    xlabel('Time (ms)','fontsize',14);
    ylabel('Amplitude (mV)','fontsize',14);
    xlim([1 N_samps]/Fs*1e3);
    
else
    % AVG WAVEFORMS
    subplot(2,3,3); hold on
    n_show = 250;
    rand_mu_set = randperm(length(mu_inds));
    if length(mu_inds) > n_show; rand_mu_set = rand_mu_set(1:n_show); end;
    plot((1:N_samps)/Fs*1e3,Spikes.V(mu_inds(rand_mu_set),:),'k','linewidth',0.1);
    for ii = 1:N_sus
        rand_su_set = randperm(length(su_inds{ii}));
        if length(rand_su_set) > n_show; rand_su_set = rand_su_set(1:n_show); end;
        plot((1:N_samps)/Fs*1e3,Spikes.V(su_inds{ii}(rand_su_set),:),'color',cmap(ii,:),'linewidth',0.1);
        %         errorbar((1:N_samps)/Fs*1e3,clusterDetails.mean_spike(:,ii+1),clusterDetails.std_spike(:,ii+1),'color',cmap(ii,:),'linewidth',2)
    end
    xlabel('Time (ms)','fontsize',14);
    ylabel('Amplitude (mV)','fontsize',14);
    xlim([1 N_samps]/Fs*1e3);
end

% ISI DIST
subplot(2,3,4); hold on
isi_bins = [logspace(log10(5e-4),log10(1),100)]*1e3;
%get ISIs
for ii = 1:N_sus
    su_spk_times = Spikes.times(su_inds{ii});
    isis = diff(su_spk_times)*1e3; %in ms
    isi_hist = histc(isis,isi_bins);
    isi_hist = isi_hist/sum(isi_hist);
    stairs(isi_bins,isi_hist,'color',cmap(ii,:))
end
set(gca,'xscale','log');
yl = ylim();
line([2 2],yl,'color','r','linestyle','--');
line([1 1],yl,'color','r');
xlim([0.5 1000]);
xlabel('ISI (ms)','fontsize',14);
ylabel('Relative frequency','fontsize',14);
%only displaying stats for best SU now
title(sprintf('%.2f - %.2f refractoriness',clusterDetails.refract(1,1),clusterDetails.refract(1,2)),'fontsize',14);

% Time v X plot
subplot(2,3,5); hold on
plot(Spikes.times(mu_inds),spike_xy(mu_inds,1),'k.')
for ii = 1:N_sus
    plot(Spikes.times(su_inds{ii}),spike_xy(su_inds{ii},1),'.','color',cmap(ii,:));
end
plot(Spikes.times(mu_inds),smooth(spike_xy(mu_inds,1),30,'rlowess'),'r','linewidth',2)
for ii = 1:N_sus
    plot(Spikes.times(su_inds{ii}),smooth(spike_xy(su_inds{ii},1),30,'rlowess'),'r','linewidth',2)
end
xlabel('Spike time (s)','fontsize',14)
ylabel('Feature projection','fontsize',14);
xlim(Spikes.times([1 end]));

% Trigger amp dist
trig_ax = linspace(min(Spikes.trig_vals)*1e3,max(Spikes.trig_vals)*1e3,200);
subplot(2,3,6)
mu_trig_hist = histc(Spikes.trig_vals(mu_inds)*1e3,trig_ax)/length(mu_inds);
stairs(trig_ax,mu_trig_hist,'k','linewidth',1);hold on
for ii = 1:N_sus
    su_trig_hist = histc(Spikes.trig_vals(su_inds{ii})*1e3,trig_ax)/length(su_inds{ii});
    stairs(trig_ax,su_trig_hist,'color',cmap(ii,:),'linewidth',1);
end
yl = ylim();
% line([sig_th sig_th]*params.thresh_sign*1e3,yl,'color','b','linewidth',2);
xlabel('Trigger voltage (mV)','fontsize',14);
ylabel('Relative frequency','fontsize',14);

