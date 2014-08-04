function [] = delete_GMM_component(block_num,probe_num)

%%
global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData

fprintf('Loading block %d Clusters\n',block_num);
cur_clust_data = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
load(cur_clust_data,'Clusters');

if Expt_name(1) == 'G'
    loadedData = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,probe_num)];
else
    sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
    if Vloaded ~= block_num
        fprintf('Loading data file %s\n',sfile_name);
        [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
        Vloaded = block_num;
    end
end

%%
cur_cluster = Clusters{probe_num};
new_cluster = cur_cluster;

[cur_cluster,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,cur_cluster,[],1);
N_spks = size(spike_xy,1);
% spike_labels = zeros(size(cur_cluster.comp_idx));
% uids = find(cur_cluster.comp_idx > 0);
% spike_labels(uids) = cur_cluster.cluster_labels(cur_cluster.comp_idx(uids));

N_sus = length(unique(cur_cluster.cluster_labels)) - 1;
spk_inds = find(cur_cluster.spike_clusts > 0);
% mu_inds = find(spike_labels == 1);
% cmap = cluster_cmap(N_sus);

clear h leg_labels
f1 = figure();
subplot(2,1,1)
hold on
plot(spike_xy(spk_inds,1),spike_xy(spk_inds,2),'k.');
cmap = cluster_cmap(length(cur_cluster.cluster_labels));
for ii = 1:length(cur_cluster.cluster_labels)
    h(ii) = plot_gaussian_2d(cur_cluster.gmm_xyMeans(ii,:)',squeeze(cur_cluster.gmm_xySigma(:,:,ii)),[2],cmap(ii,:),2);
    leg_labels{ii} = sprintf('Component %d',ii);
end
legend(h,leg_labels);
axis tight
xl = xlim(); yl = ylim();

subplot(2,1,2);hold on
[handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
set(gca,'ytick',[]);
if cur_cluster.template_it == -1
    title('Used PCs','fontsize',12);
elseif cur_cluster.template_it == -2
    title('Used voltage proj','fontsize',12);
elseif cur_cluster.template_it == -3
    title('Used energy proj','fontsize',12);
else
    title('Used template projection','fontsize',12);
end
xlim(xl); ylim(yl);
fp = get(gcf,'Position'); fp(4) = fp(4) + 600; fp(3) = fp(3) + 100;
set(gcf,'Position',fp);

target_comp = input('Which component do you want to delete?');
if ~isempty(target_comp)
    [idx,nlogl,P] = cluster(cur_cluster.gmm_fit,spike_features);
    P(:,target_comp) = 0;
    [~,new_comp_ids] = max(P,[],2);
    [~,~,new_comp_ids] = unique(new_comp_ids);
    
    new_cluster_labels = cur_cluster.cluster_labels;
    new_cluster_labels(target_comp) = [];
    N_comps = length(new_cluster_labels);

    [new_cluster.gmm_fit, new_cluster.dprime, new_cluster.comp_idx, new_cluster.cluster_labels, new_cluster.cluster_stats, outliers] = ...
        GMM_fit(Spikes.V, spike_features, N_comps, cur_cluster.params, N_comps,new_comp_ids,new_cluster_labels);
    if ~isobject(new_cluster.gmm_fit)
        fprintf('GMM fitting failed, aborting...\n');
        if ishandle(f1)
            close(f1);
        end
        return;
    end
    new_cluster.gmm_xyMeans = new_cluster.gmm_fit.mu*cur_cluster.xy_projmat;
    for ii = 1:size(new_cluster.gmm_fit.Sigma,3)
        new_cluster.gmm_xySigma(:,:,ii) = cur_cluster.xy_projmat' * squeeze(new_cluster.gmm_fit.Sigma(:,:,ii)) * cur_cluster.xy_projmat;
    end
end

spike_labels = zeros(size(new_cluster.comp_idx));
uids = find(new_cluster.comp_idx > 0);
spike_labels(uids) = new_cluster.cluster_labels(new_cluster.comp_idx(uids));
% spike_labels = new_cluster.spike_clusts;
spk_inds = find(spike_labels >= 1);
N_sus = nanmax(new_cluster.cluster_labels) - 1;
cmap = jet(N_sus);

subplot(2,1,1);hold off;
plot(spike_xy(spk_inds,1),spike_xy(spk_inds,2),'k.');
hold on
cmap = cluster_cmap(length(new_cluster.cluster_labels));
clear h leg_labels
for ii = 1:length(new_cluster.cluster_labels)
    h(ii) = plot_gaussian_2d(new_cluster.gmm_xyMeans(ii,:)',squeeze(new_cluster.gmm_xySigma(:,:,ii)),[2],cmap(ii,:),2);
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
h(1) = plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
hold on
leg_labels{1} = sprintf('Cluster %d',1);
hold on
for ii = 1:N_sus
    h(ii+1) = plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
    leg_labels{ii+1} = sprintf('Cluster %d',ii+1);
end
cmap = jet(length(new_cluster.cluster_labels));
for ii = 1:length(new_cluster.cluster_labels)
    plot_gaussian_2d(new_cluster.gmm_xyMeans(ii,:)',squeeze(new_cluster.gmm_xySigma(:,:,ii)),[2],'r',2);
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

