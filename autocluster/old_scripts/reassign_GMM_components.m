function [new_cluster,used_new] = reassign_GMM_components(sfile,cluster)


new_cluster = cluster;

[cluster,spike_features,spike_xy,Spikes] = apply_clustering(sfile,cluster,[],1);
N_spks = size(spike_xy,1);
spike_labels = zeros(size(new_cluster.comp_idx));
uids = find(new_cluster.comp_idx > 0);
spike_labels(uids) = new_cluster.cluster_labels(new_cluster.comp_idx(uids));

N_sus = length(unique(new_cluster.cluster_labels)) - 1;
mu_inds = find(spike_labels == 1);

clear h leg_labels
f1 = figure();
subplot(2,1,1)
hold on
plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
for ii = 1:N_sus
    plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'k.');
end
cmap = cluster_cmap(length(cluster.cluster_labels));
for ii = 1:length(cluster.cluster_labels)
    h(ii) = plot_gaussian_2d(cluster.gmm_xyMeans(ii,:)',squeeze(cluster.gmm_xySigma(:,:,ii)),[2],cmap(ii,:),2);
    leg_labels{ii} = sprintf('Comp %d',ii);
end
legend(h,leg_labels);
axis tight
xl = xlim(); yl = ylim();

subplot(2,1,2);hold on
[handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
set(gca,'ytick',[]);
if cluster.template_it == -1
    title('Used PCs','fontsize',12);
elseif cluster.template_it == -2
    title('Used voltage proj','fontsize',12);
elseif cluster.template_it == -3
    title('Used energy proj','fontsize',12);
else
    title('Used template projection','fontsize',12);
end
xlim(xl); ylim(yl);
% distFig()
fp = get(gcf,'Position'); fp(4) = fp(4) + 600; fp(3) = fp(3) + 100;
set(gcf,'Position',fp);




new_labels = input('Enter vector of cluster labels for each component?\n');
while length(new_labels) ~= length(cluster.cluster_labels)
    fprintf('Must specify vector of cluster labels for each component!\n');
    new_labels = input('Enter vector of cluster labels for each component?\n');
end

new_cluster.cluster_labels = new_labels;
N_sus = length(unique(new_cluster.cluster_labels)) - 1;
spike_labels = zeros(size(new_cluster.comp_idx));
uids = find(new_cluster.comp_idx > 0);
spike_labels(uids) = new_cluster.cluster_labels(new_cluster.comp_idx(uids));
mu_inds = find(spike_labels == 1);

f1 = figure();
subplot(2,1,1)
hold on
plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
cmap = cluster_cmap(N_sus);
for ii = 1:N_sus
    plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
end
for ii = 1:length(new_cluster.cluster_labels)
    h(ii) = plot_gaussian_2d(new_cluster.gmm_xyMeans(ii,:)',squeeze(new_cluster.gmm_xySigma(:,:,ii)),[2],'r',2);
    leg_labels{ii} = sprintf('Comp %d',ii);
end
legend(h,leg_labels);
axis tight
xl = xlim(); yl = ylim();

subplot(2,1,2);hold on
[handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
set(gca,'ytick',[]);
if cluster.template_it == -1
    title('Used PCs','fontsize',12);
elseif cluster.template_it == -2
    title('Used voltage proj','fontsize',12);
elseif cluster.template_it == -3
    title('Used energy proj','fontsize',12);
else
    title('Used template projection','fontsize',12);
end
xlim(xl); ylim(yl);
% distFig()
fp = get(gcf,'Position'); fp(4) = fp(4) + 600; fp(3) = fp(3) + 100;
set(gcf,'Position',fp);

keep = input('Use new cluster (y/n)?','s');
if ~strcmpi(keep,'y')
    new_cluster = cluster;
    used_new = 0;
else
    spike_stats = get_cluster_stats(Spikes.V,spike_labels);
    new_cluster.mean_spike = spike_stats.mean_spike;
    new_cluster.std_spike = spike_stats.std_spike;
    used_new = 1;
end

