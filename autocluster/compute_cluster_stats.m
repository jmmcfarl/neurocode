function clusterDetails = compute_cluster_stats(clusterDetails,Spikes,spike_features,n_null_clusts)

if nargin < 4
    n_null_clusts = 0;
end

uids = clusterDetails.comp_idx > 0;
spike_clusts = int16(nan(size(clusterDetails.comp_idx)));
spike_clusts(uids) = (clusterDetails.cluster_labels(clusterDetails.comp_idx(uids)));
clusterDetails.spike_clusts = spike_clusts(:);
% n_SUs = length(unique(spike_clusts(uids))) - 1;
n_SUs = length(unique(clusterDetails.cluster_labels))-1;
n_spks = nan(1,n_SUs+1);
refractoriness = nan(n_SUs,2);
n_spks(1) = length(spike_clusts==1);
for ii = 1:n_SUs
   n_spks(ii+1) = sum(spike_clusts==ii+1); 
   cur_isis = diff(Spikes.times(spike_clusts==ii+1))*1e3;
   refractoriness(ii,1) = sum(cur_isis < 1)/length(cur_isis)*100;
   refractoriness(ii,2) = sum(cur_isis < 2)/length(cur_isis)*100;
end
clusterDetails.n_spks = n_spks;
clusterDetails.refract = refractoriness;

clusterDetails.Ncomps = length(clusterDetails.cluster_labels);

%compute two more measures of cluster quality
[clusterDetails.Lratios,clusterDetails.iso_dists] = compute_cluster_Lratio(spike_features,clusterDetails.gmm_fit,clusterDetails.comp_idx,clusterDetails.cluster_labels);

clusterDetails.dprime = gmm_dprime(clusterDetails.gmm_fit,clusterDetails.cluster_labels);
clusterDetails.LL = clusterDetails.gmm_fit.NlogL;

cluster_stats = get_cluster_stats(Spikes.V,spike_clusts);

%compute GMM model parameters projected into XY space
gmm_xyMeans = clusterDetails.gmm_fit.mu*clusterDetails.xy_projmat;
for ii = 1:size(gmm_xyMeans,1)
    gmm_xySigma(:,:,ii) = clusterDetails.xy_projmat' * squeeze(clusterDetails.gmm_fit.Sigma(:,:,ii)) * clusterDetails.xy_projmat;
end

%add in null cluster placeholders if needed
if n_null_clusts > 0
    cluster_stats.mean_spike = cat(2,cluster_stats.mean_spike,nan(size(cluster_stats.mean_spike,1),n_null_clusts));
    cluster_stats.std_spike = cat(2,cluster_stats.std_spike,nan(size(cluster_stats.std_spike,1),n_null_clusts));
    gmm_xyMeans = cat(1,gmm_xyMeans,nan(n_null_clusts,2));
    gmm_xySigma = cat(3,gmm_xySigma,nan(2,2,n_null_clusts));
end

clusterDetails.mean_spike = cluster_stats.mean_spike;
clusterDetails.std_spike = cluster_stats.std_spike;
clusterDetails.gmm_xyMeans = gmm_xyMeans;
clusterDetails.gmm_xySigma = gmm_xySigma;
