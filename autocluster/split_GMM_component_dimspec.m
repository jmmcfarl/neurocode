function [] = split_GMM_component_dimspec(block_num,probe_num,use_features,use_2d)

%%
if nargin < 4 || isempty(use_2d)
    use_2d = false
end

global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData raw_block_nums

fprintf('Loading block %d Clusters\n',block_num);
cur_clust_data = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
load(cur_clust_data,'Clusters');

if Expt_name(1) == 'G'
    loadedData = [data_dir sprintf('/Expt%d.p%dFullV.mat',raw_block_nums(block_num),probe_num)];
else
    sfile_name = [data_dir sprintf('/Expt%dFullV.mat',raw_block_nums(block_num))];
    if Vloaded ~= raw_block_nums(block_num)
        fprintf('Loading data file %s\n',sfile_name);
        [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
        Vloaded = raw_block_nums(block_num);
    end
end

%%
cluster = Clusters{probe_num};
new_cluster = cluster;

[cluster,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,cluster,[],1);
N_spks = size(spike_xy,1);

if isnan(use_features)
    init_features = spike_xy;
    gmm_mus = cluster.gmm_xyMeans;
    gmm_Sigs = cluster.gmm_xySigma;    
else
    init_features = spike_features(:,use_features);
    gmm_mus = cluster.gmm_fit.mu(:,use_features);
    gmm_Sigs = cluster.gmm_fit.Sigma(use_features,use_features,:);
end

N_sus = length(unique(cluster.cluster_labels)) - 1;
spk_inds = find(cluster.spike_clusts > 0);

clear h leg_labels
f1 = figure();
subplot(2,1,1)
hold on
plot(init_features(spk_inds,1),init_features(spk_inds,2),'k.');
cmap = cluster_cmap(length(cluster.cluster_labels));
for ii = 1:length(cluster.cluster_labels)
    h(ii) = plot_gaussian_2d(gmm_mus(ii,:)',squeeze(gmm_Sigs(:,:,ii)),[2],cmap(ii,:),2);
    leg_labels{ii} = sprintf('Component %d',ii);
end
legend(h,leg_labels);
axis tight
xl = xlim(); yl = ylim();

subplot(2,1,2);hold on
[handles, details] = DensityPlot_jmm(init_features(:,1),init_features(:,2),'sqrtsc','ynormal','sd',[1 1]);
set(gca,'ytick',[]);
if cluster.template_it == -1
    title('Used PCs','fontsize',12);
elseif cluster.template_it == -2
    title('Used voltage proj','fontsize',12);
elseif cluster.template_it == -3
    title('Used energy proj','fontsize',12);
else
    title('Used template projection','fontsize',12);
end
xlim(xl); ylim(yl);
fp = get(gcf,'Position'); fp(4) = fp(4) + 600; fp(3) = fp(3) + 100;
set(gcf,'Position',fp);

target_comp = input('Which component do you want to split?');
if ~isempty(target_comp)
    target_inds = find(ismember(cluster.comp_idx,target_comp));
    subplot(2,1,1)
    fprintf('Create line with two points\n');
    [x,y] = ginput(2);
    sep_vec = [x(2)-x(1); y(2)-y(1)];
    orth_vec = [0 1;-1 0]*sep_vec;
    proj_data = init_features*orth_vec;
    proj_thresh = [x(1) y(1)]*orth_vec;
    above_thresh = target_inds(proj_data(target_inds) > proj_thresh);
    
    N_comps = cluster.Ncomps + 1;
    init_comp_idx = cluster.comp_idx;
    init_comp_idx(above_thresh) = max(init_comp_idx) + 1;
    init_cluster_labels = cluster.cluster_labels;
    init_cluster_labels(end+1) = max(init_cluster_labels) + 1;
    
    if use_2d
        spike_features = init_features;
    end
    
    [new_cluster.gmm_fit, new_cluster.dprime, new_cluster.comp_idx, new_cluster.cluster_labels, new_cluster.cluster_stats, outliers] = ...
        GMM_fit(Spikes.V, spike_features, N_comps, cluster.params, N_comps,init_comp_idx,init_cluster_labels);
    if ~isobject(new_cluster.gmm_fit)
        fprintf('GMM fitting failed, aborting...\n');
        if ishandle(f1)
            close(f1);
        end
        return;
    end
    if use_2d
        [new_cluster.spike_xy,new_cluster.xy_projmat] = Project_GMMfeatures(spike_features, new_cluster.gmm_fit, new_cluster.cluster_labels);
    end
    new_cluster.gmm_xyMeans = new_cluster.gmm_fit.mu*new_cluster.xy_projmat;
    for ii = 1:size(new_cluster.gmm_fit.Sigma,3)
        new_cluster.gmm_xySigma(:,:,ii) = new_cluster.xy_projmat' * squeeze(new_cluster.gmm_fit.Sigma(:,:,ii)) * new_cluster.xy_projmat;
    end
end

spike_labels = zeros(size(new_cluster.comp_idx));
uids = find(new_cluster.comp_idx > 0);
spike_labels(uids) = new_cluster.cluster_labels(new_cluster.comp_idx(uids));
% spike_labels = new_cluster.spike_clusts;
spk_inds = find(spike_labels >= 1);
N_sus = nanmax(new_cluster.cluster_labels) - 1;
cmap = jet(N_sus);

if use_2d
    gmm_mus = new_cluster.gmm_fit.mu;
    gmm_Sigs = new_cluster.gmm_fit.Sigma;
else
    gmm_mus = new_cluster.gmm_fit.mu(:,use_features);
    gmm_Sigs = new_cluster.gmm_fit.Sigma(use_features,use_features,:);
end

subplot(2,1,1);hold off;
plot(init_features(spk_inds,1),init_features(spk_inds,2),'k.');
hold on
cmap = cluster_cmap(length(new_cluster.cluster_labels));
clear h leg_labels
for ii = 1:length(new_cluster.cluster_labels)
    h(ii) = plot_gaussian_2d(gmm_mus(ii,:)',squeeze(gmm_Sigs(:,:,ii)),[2],cmap(ii,:),2);
    leg_labels{ii} = sprintf('Comp %d',ii);
end
legend(h,leg_labels);
axis tight

new_labels = input('What are the cluster labels for the components (input as vector)?\n');
while length(new_labels) ~= size(new_cluster.gmm_xyMeans,1)
    fprintf('Must assign a label to each component!\n');
    new_labels = input('What are the cluster labels for the components (input as vector)?\n');
end
new_cluster.cluster_labels = new_labels;
spike_labels(uids) = new_cluster.cluster_labels(new_cluster.comp_idx(uids));
mu_inds = find(spike_labels == 1);
N_sus = length(unique(new_cluster.cluster_labels)) - 1;
cmap = cluster_cmap(N_sus);
clear h leg_labels;

hold off
h(1) = plot(init_features(mu_inds,1),init_features(mu_inds,2),'k.');
hold on
leg_labels{1} = sprintf('Cluster %d',1);
hold on
for ii = 1:N_sus
    h(ii+1) = plot(init_features(spike_labels == ii + 1,1),init_features(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
    leg_labels{ii+1} = sprintf('Cluster %d',ii+1);
end
cmap = jet(length(new_cluster.cluster_labels));
for ii = 1:length(new_cluster.cluster_labels)
    plot_gaussian_2d(gmm_mus(ii,:)',squeeze(gmm_Sigs(:,:,ii)),[2],'r',2);
end
legend(h,leg_labels);
axis tight

keep = input('Use new cluster (y/n)?','s');
if strcmpi(keep,'y')
    %     spike_stats = get_cluster_stats(Spikes.V,spike_labels);
    %     new_cluster.mean_spike = spike_stats.mean_spike;
    %     new_cluster.std_spike = spike_stats.std_spike;
    %     new_cluster.LL = new_cluster.gmm_fit.NlogL;
    %     [new_cluster.Lratios,new_cluster.iso_dists] = compute_cluster_Lratio(spike_features,new_cluster.gmm_fit,new_cluster.comp_idx,new_cluster.cluster_labels);
    [new_cluster,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,new_cluster,[],0,Spikes);
    Clusters{probe_num} = new_cluster;
    fprintf('Saving cluster details\n');
    save(cur_clust_data,'Clusters');
    
    if block_num == new_cluster.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,new_cluster.params);
        pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        fillPage(gcf,'papersize',[14 8]);
        print(pname,'-dpng');
        
        fprintf('Saving to RefClusters\n');
        rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
        load(rclust_dat_name,'RefClusters');
        RefClusters{probe_num} = new_cluster;
        save(rclust_dat_name,'RefClusters');
    else
        resp = input('Save to RefClusters (y/n)?','s');
        if strcmpi(resp,'y')
            fprintf('Saving to RefClusters\n');
            rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
            load(rclust_dat_name,'RefClusters');
            RefClusters{probe_num} = new_cluster;
            RefClusters{probe_num}.base_block = block_num;
            save(rclust_dat_name,'RefClusters');
            
            sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,new_cluster.params);
            fillPage(gcf,'papersize',[14 8]);
            pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
            print(pname,'-dpng');
         close(sum_fig);
       end
    end
end

if ishandle(f1)
    close(f1);
end

