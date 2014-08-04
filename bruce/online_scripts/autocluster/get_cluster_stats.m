function [cluster_stats] = get_cluster_stats(X,cluster_assignments)
%[cluster_stats,new_cluster_ids] = get_cluster_stats(X,cluster_ids)

% uspks = cluster_assignments > 0;
n_clusts = nanmax(cluster_assignments);
unique_clusts = 1:n_clusts;
for cc = 1:n_clusts
    cids = find(cluster_assignments == unique_clusts(cc));
    if ~isempty(cids)
        cluster_stats.mean_spike(:,cc) = mean(X(cids,:),1);
        cluster_stats.std_spike(:,cc) = std(X(cids,:),[],1);
    else
        cluster_stats.mean_spike(:,cc) = nan(size(X,2)*size(X,3),1);
        cluster_stats.std_spike(:,cc) = nan(size(X,2)*size(X,3),1);
    end
end

