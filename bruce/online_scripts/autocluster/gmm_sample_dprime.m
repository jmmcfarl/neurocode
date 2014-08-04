function [d,Dmat]  = gmm_sample_dprime(X,comp_idx,cluster_labels)

uids = find(comp_idx > 0);
comp_ids = unique(comp_idx(uids));
n_comps = length(comp_ids);

if nargin < 3
    cluster_labels = 1:n_comps;
end
n_clusters = max(cluster_labels);

Dmat = zeros(n_clusters);
for ii = 1:(n_clusters-1)
    comp_set1 = find(cluster_labels == ii);
    for jj = (ii+1):n_clusters
        comp_set2 = find(cluster_labels == jj);
        spk_set1 = find(ismember(comp_idx,comp_set1));
        spk_set2 = find(ismember(comp_idx,comp_set2));
        Dmat(ii,jj) = mean(sqrt(mahal(X(spk_set1,:),X(spk_set2,:))));
        Dmat(jj,ii) = mean(sqrt(mahal(X(spk_set2,:),X(spk_set1,:))));
    end
end
d = mean(Dmat(Dmat > 0));