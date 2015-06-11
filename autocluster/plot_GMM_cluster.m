function [f1] = plot_GMM_cluster(block_num,probe_num,precomp_spike_data)

if nargin < 3
    precomp_spike_data = [];
end

%%
global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData raw_block_nums

fprintf('Loading block %d Clusters\n',block_num);
cur_clust_data = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
load(cur_clust_data,'Clusters');

cluster = Clusters{probe_num};
new_cluster = cluster;

if ~isempty(precomp_spike_data)
    loadedData = [];
    Vloaded = nan;
    Spikes = load_spike_data(precomp_spike_data);
    [cluster,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,cluster,[],1,Spikes);
else
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
    
    [cluster,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,cluster,[],1);
end

%%
N_spks = size(spike_xy,1);
spike_labels = zeros(size(cluster.comp_idx));
uids = find(cluster.comp_idx > 0);
spike_labels(uids) = cluster.cluster_labels(cluster.comp_idx(uids));

% N_sus = length(unique(cluster.cluster_labels)) - 1;
N_sus = nanmax(cluster.cluster_labels) - 1;
spk_inds = find(cluster.spike_clusts > 0);
mu_inds = find(spike_labels == 1);
cmap = cluster_cmap(N_sus);


clear h leg_labels
f1 = figure();
subplot(2,1,1);
h(1) = plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
hold on
leg_labels{1} = sprintf('Cluster %d',1);
hold on
leg_cnt = 2;
for ii = 1:N_sus
    if sum(spike_labels == ii + 1) > 0
    h(ii+1) = plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
    leg_labels{leg_cnt} = sprintf('Cluster %d',ii+1);
    leg_cnt = leg_cnt + 1;
    end
end
cmap = jet(length(cluster.cluster_labels));
for ii = 1:length(cluster.cluster_labels)
    plot_gaussian_2d(cluster.gmm_xyMeans(ii,:)',squeeze(cluster.gmm_xySigma(:,:,ii)),[2],'r',2);
end
legend(h,leg_labels);
axis tight
xl = xlim(); yl = ylim();

subplot(2,1,2);hold on
[handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
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


