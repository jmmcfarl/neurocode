function [xy,xy_projmat] = Project_GMMfeatures(features, gmm_obj, cluster_labels)
% [xy,xy_projmat] = Project_GMMfeatures(features, gmm_obj, cluster_labels)
% find best separating dimension, and another dimension for plotting

%%
[n_spks,D] = size(features);
n_clusters = nanmax(cluster_labels);
comp_mus = gmm_obj.mu;
comp_sigmas = gmm_obj.Sigma;
n_comps = size(comp_mus,1);

%if there are multiple components for any cluster, need to compute CLUSTER
%stats to find best separating dims. Do this as probability weighted
%averages across components of each cluster
if n_clusters < n_comps
    weights = gmm_obj.PComponents;
    cluster_mus = zeros(n_clusters,D);
    cluster_sigmas = zeros(D,D,n_clusters);
    for ii = 1:n_clusters
        comp_set = find(cluster_labels == ii);
        for jj = 1:length(comp_set)
            cluster_mus(ii,:) = cluster_mus(ii,:) + weights(comp_set(jj))*comp_mus(comp_set(jj),:);
            cluster_sigmas(:,:,ii) = cluster_sigmas(:,:,ii) + weights(comp_set(jj))*comp_sigmas(:,:,comp_set(jj));
        end
        cluster_mus(ii,:) = cluster_mus(ii,:)/sum(weights(comp_set));
        cluster_sigmas(:,:,ii) = cluster_sigmas(:,:,ii)/sum(weights(comp_set));
    end
else
    cluster_mus = comp_mus;
    cluster_sigmas = comp_sigmas;
end

%% gets 1xd best separating vector based on component means
if size(cluster_mus,1) > 2
    %use the range of cluster means to set first dim
    xy_projmat(:,1) = range(cluster_mus);
else
    %best separating dim between two gaussians
    md = diff(cluster_mus);
    S = squeeze(sum(cluster_sigmas,3));
    xy_projmat(:,1) = md/S;
end
xy_projmat(:,1) = xy_projmat(:,1)/norm(xy_projmat(:,1));

%for the second dimension pick from the first two PCs the one that is least
%correlated with the first dimension
[V,D] = eig(cov(features));
other_w = V(:,end-1:end);
cc = corr(other_w,xy_projmat);
[~,least] = min(abs(cc));
other_w = other_w(:,least);
xy_projmat(:,2) = other_w - dot(xy_projmat(:,1),other_w);
xy_projmat(:,2) = xy_projmat(:,2)/norm(xy_projmat(:,2));

%%
xy = features*xy_projmat;


