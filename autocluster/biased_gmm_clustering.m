function comp_idx = biased_gmm_clustering(features,gmm_fit,bias,cluster_labels,outlier_thresh)
% comp_idx = biased_gmm_clustering(features,gmm_fit,bias,cluster_labels,outlier_thresh)
% GMM clustering that sets all points with posterior probabilities of
% greater than 1-bias to the 'background' cluster (cluster_labels==1)
%%
if nargin < 5
    outlier_thresh = nan;
end
[comp_idx,~,P,~,M] = cluster(gmm_fit,features);
if ~isnan(outlier_thresh)
    min_mah_dists = sqrt(min(M,[],2));
    outliers = min_mah_dists > outlier_thresh;
%     comp_idx(outliers) = -1;
end

% uspks = comp_idx > 0;
% comp_idx(uspks) = cluster_labels(comp_idx(uspks));

back_comps = find(cluster_labels==1);
back_P = P(:,back_comps);
total_back_P = sum(back_P,2); %total prob each spike came from bckgnd dist
[~,best_back] = max(back_P,[],2); %most likely component of the bckgnd dist

%anything with enough probability of being from the background, set it to
%its best component
set_to_back = total_back_P >= (1-bias); 
comp_idx(set_to_back) = back_comps(best_back(set_to_back));

%reset outliers
comp_idx(outliers) = -1;

