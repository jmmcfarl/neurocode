function [new_cluster,used_new] = merge_GMM_cluster(sfile,cluster)


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
cmap = jet(N_sus);
for ii = 1:N_sus
    plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
end
cmap = jet(length(cluster.cluster_labels));
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

target_comps = input('Which components do you want to merge?');
while length(target_comps) ~= 2
    fprintf('Must specify vector of two components two merge!\n');
    target_comps = input('Which components do you want to merge?');
end

init_comp_idx = new_cluster.comp_idx;
init_comp_idx(init_comp_idx == target_comps(2)) = target_comps(1);
uids = init_comp_idx > 0;
[~,~,temp] = unique(init_comp_idx(uids));
init_comp_idx(uids) = temp;

init_cluster_labels = new_cluster.cluster_labels;
init_cluster_labels(target_comps(2)) = [];
[~,~,init_cluster_labels] = unique(init_cluster_labels);
n_clusters = length(unique(init_cluster_labels));

N_comps = new_cluster.Ncomps + 1;

[new_cluster.gmm_fit, new_cluster.gmm_distance, new_cluster.comp_idx, new_cluster.cluster_labels, new_cluster.cluster_stats, outliers] = ...
    GMM_fit(Spikes.V, spike_features, N_comps, cluster.params, N_comps,init_comp_idx,init_cluster_labels);
new_cluster.gmm_xyMeans = new_cluster.gmm_fit.mu*new_cluster.xy_projmat;
for ii = 1:size(new_cluster.gmm_fit.Sigma,3)
    new_cluster.gmm_xySigma(:,:,ii) = new_cluster.xy_projmat' * squeeze(new_cluster.gmm_fit.Sigma(:,:,ii)) * new_cluster.xy_projmat;
end
spike_labels = zeros(size(new_cluster.comp_idx));
uids = find(new_cluster.comp_idx > 0);
spike_labels(uids) = new_cluster.cluster_labels(new_cluster.comp_idx(uids));
mu_inds = find(spike_labels == 1);
N_sus = length(unique(new_cluster.cluster_labels)) - 1;
cmap = jet(N_sus);

subplot(2,1,1);hold off;
h(1) = plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
leg_labels{1} = sprintf('Cluster %d',1);
hold on
for ii = 1:N_sus
    h(ii+1) = plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
    leg_labels{ii+1} = sprintf('Cluster %d',ii+1);
end
cmap = jet(length(new_cluster.cluster_labels));
for ii = 1:length(new_cluster.cluster_labels)
    if ii == target_comps(1)
    h(ii+N_sus+1) = plot_gaussian_2d(new_cluster.gmm_xyMeans(ii,:)',squeeze(new_cluster.gmm_xySigma(:,:,ii)),[2],'r',2);
    else
    h(ii+N_sus+1) = plot_gaussian_2d(new_cluster.gmm_xyMeans(ii,:)',squeeze(new_cluster.gmm_xySigma(:,:,ii)),[2],'g',2);
    end
end
axis tight

new_label = input('What is the cluster labels for the new component (red)?\n');
new_cluster.cluster_labels(target_comps(1)) = new_label;

spike_labels(uids) = new_cluster.cluster_labels(new_cluster.comp_idx(uids));
mu_inds = find(spike_labels == 1);
N_sus = length(unique(new_cluster.cluster_labels)) - 1;
cmap = jet(N_sus);
clear h leg_labels;

hold off
h(1) = plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
hold on
leg_labels{1} = sprintf('Cluster %d',1);
hold on
for ii = 1:N_sus
    h(ii+1) = plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
    leg_labels{ii+1} = sprintf('Cluster %d',ii+1);
end
cmap = jet(length(new_cluster.cluster_labels));
for ii = 1:length(new_cluster.cluster_labels)
    plot_gaussian_2d(new_cluster.gmm_xyMeans(ii,:)',squeeze(new_cluster.gmm_xySigma(:,:,ii)),[2],cmap(ii,:),2);
%     leg_labels{ii+N_sus+1} = sprintf('Comp %d',ii);
end
legend(h,leg_labels);
axis tight

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

