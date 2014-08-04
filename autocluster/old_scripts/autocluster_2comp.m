function [clusterDetails, spike_xy, spike_features] = autocluster_2comp(SpikeV,params)
% [clusterDetails, spike_xy, spike_features] = autocluster_2comp(SpikeV,params)
% tries different sets of features, as well as iterative template matching
% to find the best 2-cluster model. Also tries 3 component (2 cluster)
% model with template scores.

%%

if nargin < 2 || isempty(params)
    params = struct();
end
if ~isfield(params,'outlier_thresh')
    params.outlier_thresh = 5; %threshold mahalanobis distance to treat points as outliers
end
if ~isfield(params,'verbose')
    params.verbose = 0; %controls level of print detail
end
if ~isfield(params,'use_best_only')
    params.use_best_only = 0; %use only the best cluster waveform as template?
end
if ~isfield(params,'cluster_bias')
    params.cluster_bias = 0.5;  %posterior probability threshold for classifying as SU
end
if ~isfield(params,'reg_lambda')
    params.reg_lambda = 1e-5;
end

%% compute pcs
[N_spks,D,N_chs] = size(SpikeV);
if N_chs > 1
    AllV = reshape(SpikeV,N_spks,D*N_Chs);
else
    AllV = SpikeV;
end
C = cov(AllV);
[pc_coeffs, E] = eig(C);
[a,b] = max(diag(E));
Eval = diag(E);
if b == 1
    fprintf('Max Eigernvector First\n');
else
    pc_coeffs = fliplr(pc_coeffs); %put largest first;
    Eval = Eval(end:-1:1);
end

%arbitrary sign convention just so that its consistent when re-applying
for j = 1:size(pc_coeffs,1)
    if sum(pc_coeffs(:,j)) < 0
        pc_coeffs(:,j) = -pc_coeffs(:,j);
    end
end
pc_scores = AllV*pc_coeffs;

%% INITIAL CLUSTERING IN PC SPACE
n_pcs = 4;
n_comps = 2;
[GMM_obj{1}, distance(1),all_comp_idx(:,1),all_clust_labels{1},cluster_stats] = ...
    GMM_fit(SpikeV,pc_scores(:,1:n_pcs),n_comps,params,n_comps);
if isobject(GMM_obj{1})
    LL(1) = GMM_obj{1}.NlogL;
else
    LL(1) = nan;
end

%find peak and valley locations for the best cluster
cluster_means = cluster_stats.mean_spike;
[~,best_ch] = max(max(abs(cluster_means)));
spk_wvfrm = cluster_means(:,best_ch);
[~,minloc] = min(spk_wvfrm);
[~,maxloc] = max(spk_wvfrm);
first_peak = min([minloc maxloc]);
second_peak = max([minloc maxloc]);
if params.verbose > 0
    fprintf('PC d-prime: %.4f\n',distance(1));
end
% clf
% plot(pc_scores(:,1),pc_scores(:,2),'k.');
% hold on
% plot(pc_scores(all_comp_idx(:,1)==1,1),pc_scores(all_comp_idx(:,1)==1,2),'r.')
% shg
   

%% INITIAL CLUSTERING IN VOLTAGE SPACE
% use_tdims = [6 13 16 20];
use_tdims = [first_peak second_peak round((second_peak+first_peak)/2) second_peak + round((second_peak-first_peak)/2)];
use_tdims(use_tdims > D) = [];
n_comps = 2;
[GMM_obj{2}, distance(2),all_comp_idx(:,2),all_clust_labels{2}] = ...
    GMM_fit(SpikeV,SpikeV(:,use_tdims),n_comps,params,n_comps);
if isobject(GMM_obj{2})
    LL(2) = GMM_obj{2}.NlogL;
else
    LL(2) = nan;
end
if params.verbose > 0
    fprintf('Voltage d-prime: %.4f\n',distance(2));
end


%% INITIAL CLUSTERING IN ENERGY SPACE
spike_energy = sqrt(sum(SpikeV.^2,2));
spike_dt_energy = sqrt(sum(diff(SpikeV,1,2).^2,2));

n_comps = 2;
[GMM_obj{3}, distance(3),all_comp_idx(:,3),all_clust_labels{3}] = ...
    GMM_fit(SpikeV,[spike_energy spike_dt_energy],n_comps,params,n_comps);
if isobject(GMM_obj{3})
    LL(3) = GMM_obj{3}.NlogL;
else
    LL(3) = nan;
end
if params.verbose > 0
    fprintf('Energy d-prime: %.4f\n',distance(3));
end

%% INITIAL CLUSTERING IN PC-ENERGY SPACE
% 
% n_comps = 2;
% [GMM_obj{4}, distance(4),all_comp_idx(:,4),all_clust_labels{4}] = ...
%     GMM_fit(SpikeV,[pc_scores(:,1:2) spike_energy spike_dt_energy],n_comps,params,n_comps);
% if isobject(GMM_obj{4})
% LL(4) = GMM_obj{4}.NlogL;
% else
%     LL(4) = nan;
% end
% if params.verbose > 0
%     fprintf('PC-Energy d-prime: %.4f\n',distance(4));
% end

%% CHOOSE BEST INITIAL CLUSTERING
[d, best] = max(distance); %best clustering so far

init_comp_idx = all_comp_idx(:,best);
init_distance = distance(best);
init_cluster_labels = 1:n_comps;
[init_cluster_labels,cluster_stats] = relabel_clusters(SpikeV,init_comp_idx,init_cluster_labels);

if params.verbose > 0
   if best == 1
      fprintf('Using PC initial clustering\n'); 
   elseif best == 2
      fprintf('Using Voltage initial clustering\n'); 
   end
end

if best == 1
    clusterDetails.init_space = 'PC';
elseif best == 2
    clusterDetails.init_space = 'voltage';
elseif best == 3
    clusterDetails.init_space = 'energy';
end
clusterDetails.pc_vecs = pc_coeffs(:,1:n_pcs);
clusterDetails.tdims = use_tdims;

%% ITERATIVE TEMPLATE FEATURE CLUSTERING
min_comp_P = min(GMM_obj{best}.PComponents);
%only do template features if both components have decent spike rates
if min_comp_P > .005 
    [template_GMM,template_comp_idx,template_cluster_labels,cluster_stats,template_distance,used_it,template_scores,templates,template_params] = ...
        iterative_template_GMM(SpikeV,init_comp_idx,init_distance,params);
else
    used_it = nan;
    templates = nan;
    template_params = struct();
    template_distance = 0;
end
%if template clustering is better than the initial clustering, use it
if template_distance > init_distance
    if params.verbose > 0
       fprintf('Using template clustering, iteration %d\n',used_it); 
    end
    gmm_fit = template_GMM;
    spike_features = template_scores;
    cluster_labels = template_cluster_labels;
    comp_idx = template_comp_idx;
    clusterDetails.fin_space = 'template';
else
    gmm_fit = GMM_obj{best};
    cluster_labels = init_cluster_labels;
    comp_idx = init_comp_idx;    
    if best == 1
        spike_features = pc_scores(:,1:n_pcs);
        used_it = -1;
        clusterDetails.fin_space = 'PC';
    elseif best == 2
        spike_features = SpikeV(:,use_tdims);
        used_it = -2;
        clusterDetails.fin_space = 'voltage';
    elseif best == 3
        spike_features = [spike_energy spike_dt_energy];
        used_it = -3;
        clusterDetails.fin_space = 'energy';
    end
end
clusterDetails.LL = gmm_fit.NlogL;
clusterDetails.dprime = gmm_dprime(gmm_fit,cluster_labels);
clusterDetails.template_it = used_it;
clusterDetails.template_params = template_params;
clusterDetails.templates = templates;

%% IMPLEMENT BIASED CLUSTERING
if params.cluster_bias ~= 0.5
    if params.verbose > 0
       fprintf('Performing biased clustering\n'); 
    end
    comp_idx = biased_gmm_clustering(spike_features,gmm_fit,params.cluster_bias,cluster_labels,params.outlier_thresh);
end

%% STORE TO clusterDetails
[cluster_labels, cluster_stats] = relabel_clusters(SpikeV,comp_idx,cluster_labels);
clusterDetails.mean_spike = cluster_stats.mean_spike;
clusterDetails.std_spike = cluster_stats.std_spike;
clusterDetails.comp_idx = int16(comp_idx);
clusterDetails.cluster_labels = cluster_labels;

%% CREATE 2d PROJECTION OF FEATURES WITH FIRST DIMENSION BEST
[spike_xy,xyproj_mat] = Project_GMMfeatures(spike_features, gmm_fit,cluster_labels);
if mean(spike_xy(su_inds,1)) < 0
    spike_xy(:,1) = -spike_xy(:,1);
    xyproj_mat(:,1) = -xyproj_mat(:,1);
end

clusterDetails.xy_projmat = xyproj_mat;

%compute GMM model parameters projected into XY space
gmm_xyMeans = gmm_fit.mu*xyproj_mat;
for ii = 1:size(gmm_fit.Sigma,3)
    gmm_xySigma(:,:,ii) = xyproj_mat' * squeeze(gmm_fit.Sigma(:,:,ii)) * xyproj_mat;
end
clusterDetails.gmm_xyMeans = gmm_xyMeans;
clusterDetails.gmm_xySigma = gmm_xySigma;

%%
clusterDetails.gmm_fit = gmm_fit;
if isobject(gmm_fit)
    clusterDetails.failed = 0;
    clusterDetails.Ncomps = size(gmm_fit.mu,1);
else
    clusterDetails.failed = 1;
    clusterDetails.Ncomps = nan;
end
