function [gmm_obj, gmm_distance, comp_idx, cluster_labels, cluster_stats, outliers] = ...
    GMM_fit(SpikeV, X, n_clusts, params, max_clusters,init_cluster_data,init_cluster_labels)

warning('off','stats:gmdistribution:FailedToConverge');

% if n_clusts > 3
%     error('not supported');
% end
if nargin < 4 || isempty(params)
    params = struct();
end
if nargin < 5
    max_clusters = 3;
end
if nargin < 6
    init_cluster_data = [];
end
if nargin < 7 
    init_cluster_labels = 1:n_clusts;
end
if ~isempty(init_cluster_data)
    use_init = true;
else
    use_init = false;
end
if ~isfield(params,'gmm_inits')
    params.gmm_inits = 50;
end
if ~isfield(params,'outlier_thresh')
    params.outlier_thresh = 7;
end
if ~isfield(params,'reg_lambda')
    params.reg_lambda = 1e-5;
end
if ~isfield(params,'max_iters')
    params.max_iters = 10;
end
if ~isfield(params,'TolFun')
    params.TolFun = 1e-4;
end
if ~isfield(params,'min_Pcomp')
    params.min_Pcomp = .0001;
end

fail = zeros(10,1);

if size(X,3) > 1
    X = reshape(X,size(X,1),size(X,2)*size(X,3));
end

if use_init %if an initial clustering is specified, just use this
    params.max_iters = 100;
    if isstruct(init_cluster_data)
        n_comps = size(init_cluster_data.mu,1);
        uids = 1:size(X,1);
        Gs{1} = fit_GMM(X(uids,:),n_comps,'Start',init_cluster_data);        
    else
        uids = init_cluster_data >= 1;
        n_comps = length(unique(init_cluster_data(uids)));
        if length(unique(init_cluster_data(uids))) == 1
            init_cluster_data = ceil(rand(size(init_cluster_data))*n_comps);
        end
        if length(unique(init_cluster_data(uids))) == 1
            Gs{1} = nan;
        else
            Gs{1} = fit_GMM(X(uids,:),n_comps,'Start',init_cluster_data(uids));
        end
    end
    ds(1) = gmm_dprime(Gs{1});
    clust_labels{1} = init_cluster_labels;
    if isobject(Gs{1})
        comp_idx = cluster(Gs{1},X);
%         ds(1) = gmm_sample_dprime(X,comp_idx,clust_labels{1});
    else
        comp_idx = nan;
        fail(1) = 1;
    end
    [L(:,1),iso_distance(:,1)] = compute_cluster_Lratio(X,Gs{1},comp_idx,clust_labels{1});
else
    %% TRY FITTING WITH RANDOM INITS
    temp_Gs = cell(params.gmm_inits,1);
    temp_ds = nan(params.gmm_inits,1);
    temp_min_pcomp = nan(params.gmm_inits,1);
    for ii = 1:params.gmm_inits
        if params.reg_lambda > 0
            temp_Gs{ii} = fit_GMM(X,n_clusts,'Regularize',params.reg_lambda,'Options',statset('MaxIter',params.max_iters,'TolFun',params.TolFun));
        else
            temp_Gs{ii} = fit_GMM(X,n_clusts,'Options',statset('MaxIter',params.max_iters,'TolFun',params.TolFun));
        end
        temp_ds(ii) = gmm_dprime(temp_Gs{ii});
        if ~isnan(temp_ds(ii))
            temp_min_pcomp(ii) = min(temp_Gs{ii}.PComponents);
        end
    end
    temp_ds(temp_min_pcomp < params.min_Pcomp) = nan;
    [ds(1),best] = max(temp_ds);
    if ~isnan(ds(1))
        Gs{1} = temp_Gs{best};
        clust_labels{1} = 1:n_clusts;
        comp_idx = cluster(Gs{1},X);
%          ds(1) = gmm_sample_dprime(X,comp_idx,clust_labels{1});
       [L(:,1),iso_distance(:,1)] = compute_cluster_Lratio(X,Gs{1},comp_idx,clust_labels{1});
        fail(1) = 0;
    else
        fail(1) = 1;
        Gs{1} = nan;
        ds(1) = nan;
        L(:,1) = nan;
        iso_distance(:,1) = nan;
    end
    %% TRY FITTING WITH PREDEFINED INIT
    %sets the initial component means to be separated along the first PC. If 3
    %clusters seperates the 3rd along the 2nd PC. Initializes the component
    %covariances to all be some fraction of the full covariance of X
    
    if size(X,2) > 1
        C = cov(X);
        [E,V] = eig(C);
        pc = E(:,end); %first PC
        pcb = E(:,end-1); %second pc
        pc = pc./max(abs(pc)); %normalize to have max 1
        pcb = pcb./max(abs(pcb));
        for j = 1:size(X,2)
            S.mu(1,j) = mean(X(:,j)) + pc(j);
            S.mu(2,j) = mean(X(:,j)) - pc(j);
            if n_clusts == 3
                S.mu(3,j) = mean(X(:,j)) + pcb(j);
            end
        end
        for j = 1:n_clusts
            S.Sigma(:,:,j) = C./sqrt(2);
        end
        
        Gs{2} = fit_GMM(X,n_clusts,'Start',S,'Regularize',params.reg_lambda);
        ds(2) = gmm_dprime(Gs{2});
        clust_labels{2} = 1:n_clusts;
        if isobject(Gs{2})
            comp_idx = cluster(Gs{2},X);
%         ds(2) = gmm_sample_dprime(X,comp_idx,clust_labels{2});
        else
            comp_idx = nan;
            fail(2) = 1;
        end
        [L(:,2),iso_distance(:,2)] = compute_cluster_Lratio(X,Gs{2},comp_idx,clust_labels{2});
    end
    
    %% TRY INITIALIZING WITH K-MEANS
    kmeans_idx = kmeans(X,n_clusts);
    Gs{3} = fit_GMM(X,n_clusts,'Regularize',params.reg_lambda,'Start',kmeans_idx);
    ds(3) = gmm_dprime(Gs{3});
    clust_labels{3} = 1:n_clusts;
    if isobject(Gs{3})
        comp_idx = cluster(Gs{3},X);
%         ds(3) = gmm_sample_dprime(X,comp_idx,clust_labels{3});
    else
        comp_idx = nan;
        fail(3) = 1;
    end
    [L(:,3),iso_distance(:,3)] = compute_cluster_Lratio(X,Gs{3},comp_idx,clust_labels{3});
    
end
%% PICK BEST FIT
if all(fail == 1)
    fprintf('No GMM fits succeeded! Aborting.\n');
    gmm_obj = [];
    gmm_distance = nan;
    comp_idx = ones(size(X,1),1);
    cluster_labels = 1;
    cluster_stats = get_cluster_stats(SpikeV,comp_idx);
    outliers = [];
    return;
end
[gmm_distance,best] = max(ds);
gmm_obj = Gs{best};
cluster_labels = clust_labels{best};
%% CHECK FOR OUTLIERS
if ~isnan(params.outlier_thresh) && isobject(gmm_obj)
    use_nclusts = size(gmm_obj.mu,1);
    [idx,nlogl,P,logpdf,M] = cluster(gmm_obj,X);
    min_mah_dists = sqrt(min(M,[],2));
    outliers = find(min_mah_dists > params.outlier_thresh);
    use_idx = setdiff(1:size(X,1),outliers);
    S.mu = gmm_obj.mu; S.Sigma = gmm_obj.Sigma;
    new_gmm_obj = fit_GMM(X(use_idx,:),use_nclusts,'Start',S,'Regularize',params.reg_lambda);
    if isobject(new_gmm_obj)
        gmm_obj = new_gmm_obj;
    end
    gmm_distance = gmm_dprime(gmm_obj,cluster_labels);
    comp_idx = cluster(gmm_obj,X);
    comp_idx(outliers) = -1;
%     gmm_distance = gmm_sample_dprime(X,comp_idx,cluster_labels);
else
    if isobject(gmm_obj)
        comp_idx = cluster(gmm_obj,X);
    else
        comp_idx = nan;
    end
    outliers = [];
end

if isobject(gmm_obj)
    %make sure smallest Pcomp is bigger than minPcomp
    if min(gmm_obj.PComponents) < params.min_Pcomp
        gmm_distance = nan;
    end
    
    [cluster_labels,cluster_stats] = relabel_clusters(SpikeV,comp_idx,cluster_labels);
else
    
    gmm_distance = nan;
    cluster_labels = nan;
    cluster_stats = nan;
end
% if isobject(gmm_obj)
%     [Lratio,iso_distance] = compute_cluster_Lratio(X,gmm_obj,comp_idx,cluster_labels);
%     gmm_distance = mean(iso_distance);
% end