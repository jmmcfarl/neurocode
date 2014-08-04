function [new_gmm_fit, gmm_distance, new_comp_idx, new_clust_labels] = ...
    try_backgnd_splitting(gmm_fit,SpikeV,spike_features,spike_xy,xyproj_mat,comp_idx,clust_labels,params)


su_comp = find(clust_labels == 2);
su_set = find(comp_idx == su_comp);

med_pt = median(spike_xy(su_set,1));
% med_pt = prctile(spike_xy(su_set,1),90);

new_sucomp = su_set(spike_xy(su_set,1) > med_pt);
new_backcomp = su_set(spike_xy(su_set,1) < med_pt);

new_comp_idx = comp_idx;
new_comp_idx(new_sucomp) = length(clust_labels) + 1;
new_comp_idx(new_backcomp) = su_comp;
new_clust_labels = clust_labels;
new_clust_labels(su_comp) = 1;
new_clust_labels(end+1) = 2;

[new_gmm_fit, gmm_distance, new_comp_idx, new_clust_labels] = ...
    GMM_fit(SpikeV, spike_features, [], params, [],new_comp_idx,new_clust_labels);

if isobject(new_gmm_fit)
    %make sure smallest Pcomp is bigger than minPcomp
    if min(new_gmm_fit.PComponents) < params.min_Pcomp
        gmm_distance = nan;
    end
    
    [new_clust_labels,cluster_stats] = relabel_clusters(SpikeV,new_comp_idx,1:length(new_clust_labels));
    %     [~,new_su_comp] = max(new_clust_labels);
    new_su_comp = new_clust_labels(end);
    new_clust_labels = ones(size(new_clust_labels));
    new_clust_labels(new_su_comp) = 2;
%     gmm_distance = gmm_sample_dprime(spike_features,new_comp_idx,new_clust_labels);
    gmm_distance = gmm_dprime(new_gmm_fit,new_clust_labels);

%     [Lratio,iso_distance] = compute_cluster_Lratio(spike_features,new_gmm_fit,new_comp_idx,new_clust_labels);
%     gmm_distance = mean(iso_distance);

else
    gmm_distance = nan;
    new_clust_labels = nan;
    cluster_stats = nan;
end



   
    
    