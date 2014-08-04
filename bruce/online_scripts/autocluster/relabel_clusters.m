function [cluster_labels, cluster_stats] = relabel_clusters(SpikeV,comp_idx,cluster_labels)
% [cluster_labels, cluster_stats] = relabel_clusters(SpikeV,comp_idx,cluster_labels)
% determines which of two clusters is the background, and sets the cluster
% label to 1 for background and 2 for unit. 

if isnan(comp_idx)
    cluster_labels = nan;
    cluster_stats = nan;
    return;
end
cluster_assignments = zeros(size(comp_idx));
uspks = comp_idx > 0; %non-outlier spikes
unique_cids = unique(comp_idx(uspks)); %unique Gaussian components
n_clusters = nanmax(cluster_labels);
for ii = 1:n_clusters
    component_set = find(cluster_labels == ii);
    cluster_assignments(ismember(comp_idx,component_set)) = ii;
end

[cluster_stats] = get_cluster_stats(SpikeV,cluster_assignments);

max_amps = nanmax(abs(cluster_stats.mean_spike));
[~,label_order] = sort(max_amps); %sort cluster labels by increasing amplitude of average waveform

%if there's only one cluster make another as a placeholder (not sure if we
%need this)...
if size(cluster_stats.mean_spike,2) < 2
    label_order(2) = 2;
    cluster_stats.mean_spike(:,2) = nan;
    cluster_stats.std_spike(:,2) = nan;
end

cluster_labels = label_order(cluster_labels);

cluster_stats.mean_spike = cluster_stats.mean_spike(:,label_order);
cluster_stats.std_spike = cluster_stats.std_spike(:,label_order);