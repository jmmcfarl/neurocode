function [GMM_obj,comp_idx,clust_labels,cluster_stats,distance,used_it,template_scores,templates,template_params] = ...
    iterative_template_GMM(Spikes,init_comp_idx,init_cluster_labels,init_distance,params)
% [GMM_obj,comp_idx,distance] = iterative_template_GMM(Spikes,init_comp_idx,t_derivs,use_ms,verbose)
% NOTE, this function must start from an initial clustering with 2
% components (representing two clusters)

%%
if nargin < 4
    init_distance = Inf;
end
if nargin < 5
    params = struct();
end
if ~isfield(params,'n_used_templates')
    params.n_used_templates = 4;
end
if ~isfield(params,'max_iter')
    params.max_iter = 10;
end
if ~isfield(params,'deriv_types')
   params.deriv_types = [0 1]; %default use both 0th and 1st derivative templates
end
if ~isfield(params,'clust_eps')
    params.clust_eps = 0.01;
end
if ~isfield(params,'use_ms')
    params.use_ms = 1; %default use mean-subtraction
end
if ~isfield(params,'outlier_thresh')
    params.outlier_thresh = nan;
end
if ~isfield(params,'verbose')
    params.verbose = 0;
end
if ~isfield(params,'use_best_only')
    params.use_best_only = 0;
end
if any(params.deriv_types > 1)
    error('Only first derivative supported');
end
warning('off','stats:gmdistribution:FailedToConverge');
warning('off','stats:gmdistribution:MaxIterations');

[N_spks,D,N_chs] = size(Spikes.V);
n_comps = length(unique(init_comp_idx(init_comp_idx > 0)));
% if n_comps > 2
%     error('Must start with 2-component model for iterative template clustering');
% end

n_clusters = 2;

%% DO initial template GMM 
if params.verbose > 0
    fprintf('Fitting GMM to initial template scores\n');
end
cluster_assignments = zeros(size(init_comp_idx));
uspks = init_comp_idx > 0; %non-outlier spikes
unique_cids = unique(init_comp_idx(uspks)); %unique Gaussian components
n_clusters = nanmax(init_cluster_labels);
for ii = 1:n_clusters
    component_set = find(init_cluster_labels == ii);
    cluster_assignments(ismember(init_comp_idx,component_set)) = ii;
end

cluster_stats = get_cluster_stats(Spikes.V,cluster_assignments); %initial cluster stats
cluster_means = reshape(cluster_stats.mean_spike,[D size(cluster_stats.mean_spike,2) N_chs]);
[peak_amps,peak_locs] = max(cluster_means);
peak_amps = squeeze(peak_amps);
if N_chs == 1
    [~,best_clust] = max(peak_amps);
    best_ch = 1;
else
    peakloc = find(peak_amps == max(peak_amps(:)));
    [best_clust,best_ch] = ind2sub([2 N_chs],peakloc);
end

if params.use_best_only
   cluster_means = cluster_means(:,best_clust,:); 
   n_used_clusters = 1;
else
    n_used_clusters = n_clusters;
end

%create the templates
n_templates = length(params.deriv_types)*N_chs*n_used_clusters; %one template for each derivative type, channel, and cluster
templates = zeros(D,n_templates);
channels = zeros(n_templates,1);
t_derivs = zeros(n_templates,1);
cnt = 1;
if any(params.deriv_types == 0)
    for cc = 1:N_chs
        templates(:,cnt:(cnt+n_used_clusters-1)) = squeeze(cluster_means(:,:,cc));
        channels(cnt:(cnt+n_used_clusters-1)) = cc;
        t_derivs(cnt:(cnt+n_used_clusters-1)) = 0;
        cnt = cnt + n_used_clusters;
    end
end
if any(params.deriv_types == 1)
    for cc = 1:N_chs
        templates(2:end,cnt:(cnt+n_used_clusters-1)) = squeeze(diff(cluster_means(:,:,cc)));
        channels(cnt:(cnt+n_used_clusters-1)) = cc;
        t_derivs(cnt:(cnt+n_used_clusters-1)) = 1;
        cnt = cnt + n_used_clusters;
    end
end

template_params = struct('channels',channels,'t_derivs',t_derivs,'mean_sub',...
    params.use_ms,'use_best_only',params.use_best_only);

[template_scores] = get_template_scores(Spikes.V,templates,template_params);
template_ks = nan(n_templates,1);
for ii = 1:n_templates
   template_ks(ii) = lillie_KSstat(template_scores(:,ii)); 
end
[~,ks_ord] = sort(template_ks,'descend');
use_template_set = ks_ord(1:params.n_used_templates);
template_params.channels = template_params.channels(use_template_set);
template_params.t_derivs = template_params.t_derivs(use_template_set);
template_scores = template_scores(:,use_template_set);
templates = templates(:,use_template_set);


%fit GMM with two clusters and two components. Also try adding a third component
[GMM_obj{1}, distance(1),all_comp_idx(:,1), all_clust_labels{1}, cluster_stats] = GMM_fit(Spikes.V,template_scores,n_clusters,params,n_clusters+1);

if isobject(GMM_obj{1})
        
    new_GMM_obj = GMM_obj{1};
    new_comp_idx = all_comp_idx(:,1);
    new_clust_labels = all_clust_labels{1};
    cur_n_back_comps = 1;
    while cur_n_back_comps < params.max_back_comps
        cur_n_back_comps = cur_n_back_comps + 1;
        fprintf('Trying background split %d of %d\n',cur_n_back_comps,params.max_back_comps);
        
        [spike_xy,xyproj_mat] = Project_GMMfeatures(template_scores, new_GMM_obj,new_clust_labels);
        %implement a sign convention that positive values of the first dimension
        %correspond to higher spike amplitudes
        bb = corr(spike_xy(:,1),abs(Spikes.trig_vals(:)));
        if bb < 0
            spike_xy(:,1) = -spike_xy(:,1);
            xyproj_mat(:,1) = -xyproj_mat(:,1);
        end
        
        [new_GMM_obj, new_distance,new_comp_idx, new_clust_labels] = ...
            try_backgnd_splitting(new_GMM_obj,Spikes.V,template_scores,spike_xy,xyproj_mat,new_comp_idx,new_clust_labels,params);
                
        fprintf('Orig: %.3f New: %.3f \n',distance(1),new_distance);
        if new_distance > distance(1)
            GMM_obj{1} = new_GMM_obj;
            distance(1) = new_distance;
            all_comp_idx(:,1) = new_comp_idx;
            all_clust_labels{1} = new_clust_labels;
        end
    end
end

%how much did the component assignments change
comp_change(1) = sum(all_comp_idx(:,1) ~= init_comp_idx)/N_spks;
if params.verbose > 1
    fprintf('Iteration 1, comp_change: %.4f   GMM d-prime: %.4f\n',comp_change(end),distance(end));
end

%%
it_cnt = 2;
if distance(1) > init_distance
    better_clust = true;
else
    fprintf('Template clustering worse than initial, keeping initial\n');
    comp_idx = init_comp_idx(:,1);
    clust_labels = all_clust_labels{1};
    GMM_obj = GMM_obj{1};
    used_it = 0;
    return;
end
%while the change in component assignments is large enough, and the cluster
%quality is improving, keep iterating
while comp_change(end) > params.clust_eps && better_clust
    
    %compute new templates
    prev_templates = templates;
    cluster_means = reshape(cluster_stats.mean_spike,[D size(cluster_stats.mean_spike,2) N_chs]);
    [peak_amps,peak_locs] = max(cluster_means);
    peak_amps = squeeze(peak_amps);
    if N_chs == 1
        [~,best_clust] = max(peak_amps);
        best_ch = 1;
    else
        peakloc = find(peak_amps == max(peak_amps(:)));
        [best_clust,best_ch] = ind2sub([2 N_chs],peakloc);
    end
    
    if params.use_best_only
        cluster_means = cluster_means(:,best_clust,:);
        n_used_clusters = 1;
    else
        n_used_clusters = n_clusters;
    end
    
    cnt = 1;
    if any(params.deriv_types == 0)
        for cc = 1:N_chs
            templates(:,cnt:(cnt+n_used_clusters-1)) = squeeze(cluster_means(:,:,cc));
            cnt = cnt + n_used_clusters;
        end
    end
    if any(params.deriv_types == 1)
        for cc = 1:N_chs
            templates(2:end,cnt:(cnt+n_used_clusters-1)) = squeeze(diff(cluster_means(:,:,cc)));
            cnt = cnt + n_used_clusters;
        end
    end
    templates = templates(:,use_template_set);
    prev_template_scores = template_scores;
    [template_scores] = get_template_scores(Spikes.V,templates,template_params);
    
    %now fit GMM using the previous component assignments as a starting
    %point
    [GMM_obj{it_cnt}, distance(it_cnt),all_comp_idx(:,it_cnt), all_clust_labels{it_cnt}, cluster_stats] = ...
        GMM_fit(Spikes.V,template_scores,n_clusters,params,0,all_comp_idx(:,it_cnt-1),all_clust_labels{it_cnt-1});
        
    comp_change(it_cnt) = sum(all_comp_idx(:,it_cnt) ~= all_comp_idx(:,it_cnt-1))/N_spks;
    if distance(it_cnt) <= distance(it_cnt - 1) || isnan(distance(it_cnt))
        better_clust = false;
    end
    if params.verbose > 1
        fprintf('Iteration %d, comp_change: %.4f   GMM d-prime: %.4f\n',it_cnt,comp_change(end),distance(end));
    end
    it_cnt = it_cnt + 1;
end

%%
if better_clust
    GMM_obj = GMM_obj{end};
    comp_idx = all_comp_idx(:,end);
    clust_labels = all_clust_labels{end};
    distance = distance(end);
    used_it = it_cnt;
else
    GMM_obj = GMM_obj{end-1};
    template_scores = prev_template_scores;
    comp_idx = all_comp_idx(:,end-1);
    clust_labels = all_clust_labels{end-1};
    distance = distance(end-1);
    templates = prev_templates;
    used_it = it_cnt-1;
end
