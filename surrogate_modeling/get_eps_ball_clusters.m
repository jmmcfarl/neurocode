function cluster_id = get_eps_ball_clusters(X,eps)

cluster_id = zeros(size(X));
cluster_id(1) = 1;
cur_n_clusts = 1;

clust_avgs = X;
for i = 2:length(X)
    clust_set = 1:(i-1);
    dist_vec = sqrt((X(i) - clust_avgs(clust_set)).^2);
    [a,b] = min(dist_vec);
%     b = clust_set(b);
    if a < eps & cluster_id(b) ~= 0
        cluster_id(i) = cluster_id(b);
        %move points in cluster to cluster avg
        cur_set = find(cluster_id == cluster_id(b)); 
        clust_avgs(cur_set) = mean(X(cur_set));
    else
        cur_n_clusts = cur_n_clusts + 1;
        cluster_id(i) = cur_n_clusts;
    end
end

