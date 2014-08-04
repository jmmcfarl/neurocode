function [d,Dclust]  = gmm_dprime(G,clust_labels)
% [d]  = gmm_dprime(G, varargin)
% calcualte drpime between two Gaussians in gmdistribution fit

%%
if nargin < 2
    clust_labels = [];
end
if ~isobject(G)
    d = nan;
    Dclust = nan;
    return;
end
nc =size(G.mu,1);

%%
distance = mahal(G,G.mu);
Dmat = zeros(nc,nc);
for j = 1:nc
    for k = 1:j-1
        Dmat(j,k) = sqrt(2./((1./distance(j,k))+(1./distance(k,j))));
        Dmat(k,j) = Dmat(j,k);
    end
end
if nc == 2
    d = gmm_distance(G);
    Dclust = Dmat;
elseif isempty(clust_labels)
    d  = mean(Dmat(Dmat>0));
    Dclust = Dmat;
else
    weights = G.PComponents;
    
    unique_clabels = unique(clust_labels);
    for jj = 1:length(unique_clabels)
       cur_set = find(clust_labels == unique_clabels(jj));
       weights(cur_set) = weights(cur_set)/sum(weights(cur_set));
    end
    Dclust = zeros(length(unique_clabels));
    for jj = 1:(length(unique_clabels)-1)
        gauss_set1 = find(clust_labels == unique_clabels(jj));
        for kk = (jj+1):length(unique_clabels)
            gauss_set2 = find(clust_labels == unique_clabels(kk));
%             Dclust(jj,kk) = weights(gauss_set1)*Dmat(gauss_set1,gauss_set2)*weights(gauss_set2)';
            Dclust(jj,kk) = mean(reshape(Dmat(gauss_set1,gauss_set2),1,[]));
        end
    end
    d  = mean(Dclust(Dclust>0));
end

%%
