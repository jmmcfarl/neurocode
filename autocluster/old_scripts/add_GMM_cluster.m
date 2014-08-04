function [new_cluster,used_new] = add_GMM_cluster(sfile,cluster)


new_cluster = cluster;

[cluster,spike_features,spike_xy,Spikes] = apply_clustering(sfile,cluster,[],1);
N_spks = size(spike_xy,1);
N_sus = length(unique(cluster.cluster_labels))-1;
mu_inds = find(ismember(cluster.comp_idx,find(cluster.cluster_labels == 1)));

f1 = figure();
subplot(2,1,1)
hold on
cmap = cluster_cmap(N_sus);
h(1) = plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
for ii = 1:N_sus
    su_inds = find(ismember(cluster.comp_idx,find(cluster.cluster_labels == ii+1)));
    h(ii+1) = plot(spike_xy(su_inds,1),spike_xy(su_inds,2),'.','color',cmap(ii,:));
end
clear leg_labels
for ii = 1:length(cluster.cluster_labels)
    hh(ii) = plot_gaussian_2d(cluster.gmm_xyMeans(ii,:)',squeeze(cluster.gmm_xySigma(:,:,ii)),[2],'r',2);
    leg_labels{ii} = sprintf('Cluster %d',ii);
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

target_clust = input('Which cluster do you want to split?');
target_comps = find(cluster.cluster_labels == target_clust);
target_inds = find(ismember(cluster.comp_idx,target_comps));
subplot(2,1,1)
fprintf('Create line with two points\n');
[x,y] = ginput(2);
sep_vec = [x(2)-x(1); y(2)-y(1)];
orth_vec = [0 1;-1 0]*sep_vec;
proj_data = spike_xy*orth_vec;
proj_thresh = [x(1) y(1)]*orth_vec;
above_thresh = target_inds(proj_data(target_inds) > proj_thresh);
below_thresh = target_inds(proj_data(target_inds) <= proj_thresh);

% N_comps = cluster.Ncomps + 1;
init_comp_idx = cluster.comp_idx;
    init_cluster_labels = cluster.cluster_labels;
if length(target_comps) == 1
    init_comp_idx(above_thresh) = max(init_comp_idx) + 1;
    init_cluster_labels(end+1) = max(init_cluster_labels) + 1;
else
    init_comp_idx(above_thresh) = target_comps(2);
    init_comp_idx(below_thresh) = target_comps(1);
    init_cluster_labels(target_comps(2)) = max(init_cluster_labels) + 1;
end

[new_cluster.gmm_fit, new_cluster.gmm_distance, new_cluster.comp_idx, new_cluster.cluster_labels, new_cluster.cluster_stats, outliers] = ...
    GMM_fit(Spikes.V, spike_features, [], cluster.params, [],init_comp_idx,init_cluster_labels);
new_cluster.gmm_xyMeans = new_cluster.gmm_fit.mu*new_cluster.xy_projmat;
for ii = 1:size(new_cluster.gmm_fit.Sigma,3)
    new_cluster.gmm_xySigma(:,:,ii) = new_cluster.xy_projmat' * squeeze(new_cluster.gmm_fit.Sigma(:,:,ii)) * new_cluster.xy_projmat;
end

spike_labels = zeros(size(new_cluster.comp_idx));
uids = find(new_cluster.comp_idx > 0);
spike_labels(uids) = new_cluster.cluster_labels(new_cluster.comp_idx(uids));
N_sus = length(unique(new_cluster.cluster_labels))-1;
cmap = cluster_cmap(N_sus);
subplot(2,1,1);hold off;
mu_inds = find(ismember(new_cluster.comp_idx,find(new_cluster.cluster_labels == 1)));
plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
hold on
for ii = 1:N_sus
    plot(spike_xy(spike_labels == ii+1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
end
clear h leg_labels
for ii = 1:length(new_cluster.cluster_labels)
    h(ii) = plot_gaussian_2d(new_cluster.gmm_xyMeans(ii,:)',squeeze(new_cluster.gmm_xySigma(:,:,ii)),[2],'r',2);
end

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

