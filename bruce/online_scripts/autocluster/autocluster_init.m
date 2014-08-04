function [clusterDetails, spike_xy, spike_features] = autocluster_init(Spikes,params)
% [clusterDetails, spike_xy, spike_features] = autocluster_init(Spikes,params)
% tries different sets of features, as well as iterative template matching
% to find the best 2-cluster model. Also tries 3 component (2 cluster)
% model with template scores.

%% SETS DEFAULT PARAMETERS NEEDED. Note, many of these are redundant with defaults set in detect_and_cluster_init. So, these are only used if this function is being used as a standalone.

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
if ~isfield(params,'max_back_comps')
    params.max_back_comps = 3;
end
if ~isfield(params,'n_pcs')
params.n_pcs = 4;
end
if ~isfield(params,'n_tdims')
params.n_tdims = 4;
end

%% compute pcs
[N_spks,D,N_chs] = size(Spikes.V);
if N_chs > 1
    AllV = reshape(Spikes.V,N_spks,D*N_chs);
else
    AllV = Spikes.V;
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
if ismember(1,params.try_features)
    n_comps = 2;
    pc_ks = nan(D,1);
    for ii = 1:D
        pc_ks(ii) = lillie_KSstat(pc_scores(:,ii));
    end
    [~,use_pcs] = sort(pc_ks,'descend');
    use_pcs = use_pcs(1:params.n_pcs);
    
    [GMM_obj{1}, distance(1),all_comp_idx(:,1),all_clust_labels{1},cluster_stats] = ...
        GMM_fit(Spikes.V,pc_scores(:,use_pcs),n_comps,params,n_comps);
    if isobject(GMM_obj{1})
        LL(1) = GMM_obj{1}.NlogL;
        
        new_GMM_obj = GMM_obj{1};
        new_comp_idx = all_comp_idx(:,1);
        new_clust_labels = all_clust_labels{1};
        cur_n_back_comps = 1;
        while cur_n_back_comps < params.max_back_comps
            cur_n_back_comps = cur_n_back_comps + 1;
            fprintf('Trying background split %d of %d\n',cur_n_back_comps,params.max_back_comps);
            
            [spike_xy,xyproj_mat] = Project_GMMfeatures(pc_scores(:,use_pcs), new_GMM_obj,new_clust_labels);
            %implement a sign convention that positive values of the first dimension
            %correspond to higher spike amplitudes
            bb = corr(spike_xy(:,1),abs(Spikes.trig_vals(:)));
            if bb < 0
                spike_xy(:,1) = -spike_xy(:,1);
                xyproj_mat(:,1) = -xyproj_mat(:,1);
            end
            [new_GMM_obj, new_distance,new_comp_idx, new_clust_labels] = ...
                try_backgnd_splitting(new_GMM_obj,Spikes.V,pc_scores(:,use_pcs),spike_xy,xyproj_mat,new_comp_idx,new_clust_labels,params);
           
            fprintf('Orig: %.3f New: %.3f \n',distance(1),new_distance);
            if new_distance > distance(1)
                GMM_obj{1} = new_GMM_obj;
                distance(1) = new_distance;
                all_comp_idx(:,1) = new_comp_idx;
                all_clust_labels{1} = new_clust_labels;
            end
        end
    else
        LL(1) = nan;
    end
    
    min_comp_P = min(GMM_obj{1}.PComponents);
    %only do template features if both components have decent spike rates
    if min_comp_P < params.min_Pcomp
        distance(1) = nan;
    end
    
    %find peak and valley locations for the best cluster
    cluster_means = reshape(cluster_stats.mean_spike,[D N_chs 2]);
    [peak_amps,peak_locs] = max(cluster_means);
    peak_amps = squeeze(peak_amps);
    if N_chs == 1
        [~,best_clust] = max(peak_amps);
        best_ch = 1;
    else
        peakloc = find(peak_amps == max(peak_amps(:)));
        [best_ch,best_clust] = ind2sub([N_chs 2],peakloc);
    end
    spk_wvfrm = squeeze(cluster_means(:,best_ch,best_clust));
    [~,minloc] = min(spk_wvfrm);
    [~,maxloc] = max(spk_wvfrm);
    first_peak = min([minloc maxloc]);
    second_peak = max([minloc maxloc]);
    if params.verbose > 0
        fprintf('PC d-prime: %.4f\n',distance(1));
    end
    clusterDetails.pc_vecs = pc_coeffs(:,use_pcs);
    clusterDetails.use_pcs = use_pcs;
else
    LL(1) = nan;
    distance(1) = nan;
end
%% INITIAL CLUSTERING IN VOLTAGE SPACE
if ismember(2,params.try_features)
    tdim_ks = nan(D*N_chs,1);
    for ii = 1:D*N_chs
        tdim_ks(ii) = lillie_KSstat(AllV(:,ii));
    end
    peak_locs = find(diff(sign(diff(tdim_ks))) < 0);
    if ~isempty(peak_locs)
        peak_locs = peak_locs + 1;
        peak_amps = tdim_ks(peak_locs);
        [~,use_peaks] = sort(peak_amps,'descend');
        peak_locs = peak_locs(use_peaks);
        if length(peak_locs) > params.n_tdims
            peak_locs = peak_locs(1:params.n_tdims);
        end
    else
        peak_locs = [first_peak second_peak round((second_peak+first_peak)/2) second_peak + round((second_peak-first_peak)/2)];
    end
    use_tdims = peak_locs;
    use_tdims(use_tdims > D*N_chs) = [];
    
    n_comps = 2;
    [GMM_obj{2}, distance(2),all_comp_idx(:,2),all_clust_labels{2}] = ...
        GMM_fit(Spikes.V,AllV(:,use_tdims),n_comps,params,n_comps);
    if isobject(GMM_obj{2})
        LL(2) = GMM_obj{2}.NlogL;
        
        new_GMM_obj = GMM_obj{2};
        new_comp_idx = all_comp_idx(:,2);
        new_clust_labels = all_clust_labels{2};
        cur_n_back_comps = 1;
        while cur_n_back_comps < params.max_back_comps
            cur_n_back_comps = cur_n_back_comps + 1;
            fprintf('Trying background split %d of %d\n',cur_n_back_comps,params.max_back_comps);
            
            [spike_xy,xyproj_mat] = Project_GMMfeatures(AllV(:,use_tdims), new_GMM_obj,new_clust_labels);
            %implement a sign convention that positive values of the first dimension
            %correspond to higher spike amplitudes
            bb = corr(spike_xy(:,1),abs(Spikes.trig_vals(:)));
            if bb < 0
                spike_xy(:,1) = -spike_xy(:,1);
                xyproj_mat(:,1) = -xyproj_mat(:,1);
            end
            [new_GMM_obj, new_distance,new_comp_idx, new_clust_labels] = ...
                try_backgnd_splitting(new_GMM_obj,Spikes.V,AllV(:,use_tdims),spike_xy,xyproj_mat,new_comp_idx,new_clust_labels,params);
            
            fprintf('Orig: %.3f New: %.3f \n',distance(2),new_distance);
            if new_distance > distance(2)
                GMM_obj{2} = new_GMM_obj;
                distance(2) = new_distance;
                all_comp_idx(:,2) = new_comp_idx;
                all_clust_labels{2} = new_clust_labels;
            end
        end
        
    else
        LL(2) = nan;
    end
        
    min_comp_P = min(GMM_obj{2}.PComponents);
    %only do template features if both components have decent spike rates
    if min_comp_P < params.min_Pcomp
        distance(2) = nan;
    end
    
    if params.verbose > 0
        fprintf('Voltage d-prime: %.4f\n',distance(2));
    end
clusterDetails.tdims = use_tdims;
else
    LL(2) = nan;
    distance(2) = nan;
end

%% INITIAL CLUSTERING IN ENERGY SPACE
if ismember(3,params.try_features)
    spike_energy = squeeze(sqrt(sum(Spikes.V.^2,2)));
    spike_dt_energy = squeeze(sqrt(sum(diff(Spikes.V,1,2).^2,2)));
    energy_features = [spike_energy spike_dt_energy];
    
    n_comps = 2;
    [GMM_obj{3}, distance(3),all_comp_idx(:,3),all_clust_labels{3}] = ...
        GMM_fit(Spikes.V,energy_features,n_comps,params,n_comps);
    if isobject(GMM_obj{3})
        LL(3) = GMM_obj{3}.NlogL;
        
        new_GMM_obj = GMM_obj{3};
        new_comp_idx = all_comp_idx(:,3);
        new_clust_labels = all_clust_labels{3};
        cur_n_back_comps = 1;
        while cur_n_back_comps < params.max_back_comps
            cur_n_back_comps = cur_n_back_comps + 1;
            fprintf('Trying background split %d of %d\n',cur_n_back_comps,params.max_back_comps);
            
            [spike_xy,xyproj_mat] = Project_GMMfeatures(energy_features, new_GMM_obj,new_clust_labels);
            %implement a sign convention that positive values of the first dimension
            %correspond to higher spike amplitudes
            bb = corr(spike_xy(:,1),abs(Spikes.trig_vals(:)));
            if bb < 0
                spike_xy(:,1) = -spike_xy(:,1);
                xyproj_mat(:,1) = -xyproj_mat(:,1);
            end
            [new_GMM_obj, new_distance,new_comp_idx, new_clust_labels] = ...
                try_backgnd_splitting(new_GMM_obj,Spikes.V,energy_features,spike_xy,xyproj_mat,new_comp_idx,new_clust_labels,params);
            
            fprintf('Orig: %.3f New: %.3f \n',distance(3),new_distance);
            if new_distance > distance(3)
                GMM_obj{3} = new_GMM_obj;
                distance(3) = new_distance;
                all_comp_idx(:,3) = new_comp_idx;
                all_clust_labels{3} = new_clust_labels;
            end            
        end
    else
        LL(3) = nan;
    end
    min_comp_P = min(GMM_obj{3}.PComponents);
    %only do template features if both components have decent spike rates
    if min_comp_P < params.min_Pcomp
        distance(3) = nan;
    end

    if params.verbose > 0
        fprintf('Energy d-prime: %.4f\n',distance(3));
    end
else
    LL(3) = nan;
    distance(3) = nan;
end
%% INITIAL CLUSTERING IN PC-ENERGY SPACE
% 
% n_comps = 2;
% [GMM_obj{4}, distance(4),all_comp_idx(:,4),all_clust_labels{4}] = ...
%     GMM_fit(Spikes.V,[pc_scores(:,1:2) spike_energy spike_dt_energy],n_comps,params,n_comps);
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
init_cluster_labels = all_clust_labels{best};
[init_cluster_labels,cluster_stats] = relabel_clusters(Spikes.V,init_comp_idx,init_cluster_labels);

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

%% ITERATIVE TEMPLATE FEATURE CLUSTERING
min_comp_P = min(GMM_obj{best}.PComponents);
%only do template features if both components have decent spike rates
if min_comp_P > .005 && ismember(4,params.try_features)
    [template_GMM,template_comp_idx,template_cluster_labels,cluster_stats,template_distance,used_it,template_scores,templates,template_params] = ...
        iterative_template_GMM(Spikes,init_comp_idx,init_cluster_labels,init_distance,params);
    min_comp_P = min(template_GMM.PComponents);
    %only do template features if both components have decent spike rates
    if min_comp_P < params.min_Pcomp
        template_distance = nan;
    end
    
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
        spike_features = pc_scores(:,use_pcs);
        used_it = -1;
        clusterDetails.fin_space = 'PC';
    elseif best == 2
        spike_features = AllV(:,use_tdims);
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
[cluster_labels, cluster_stats] = relabel_clusters(Spikes.V,comp_idx,cluster_labels);
clusterDetails.mean_spike = cluster_stats.mean_spike;
clusterDetails.std_spike = cluster_stats.std_spike;
clusterDetails.comp_idx = int16(comp_idx);
clusterDetails.cluster_labels = cluster_labels;

%% CREATE 2d PROJECTION OF FEATURES 
[spike_xy,xyproj_mat] = Project_GMMfeatures(spike_features, gmm_fit,cluster_labels);

%implement a sign convention that positive values of the first dimension
%correspond to higher spike amplitudes
bb = corr(spike_xy(:,1),abs(Spikes.trig_vals(:)));
if bb < 0
    spike_xy(:,1) = -spike_xy(:,1);
    xyproj_mat(:,1) = -xyproj_mat(:,1);
end

clusterDetails.xy_projmat = xyproj_mat;
clusterDetails.spike_xy = spike_xy;

%compute GMM model parameters projected into XY space. Useful for plotting
%gaussian components
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
