function [full_scatter_fig] = regenerate_allblock_xyscatters_v2(ref_block,probe_nums,target_blocks,clust_nums)

global base_save_dir full_save_dir


fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name,'RefClusters');

full_scatter_fig = figure('visible','off');
n_blocks = length(target_blocks);
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);

for bb = 1:length(target_blocks)
    fprintf('Block %d\n',target_blocks(bb));
    
    probe_num = probe_nums(target_blocks(bb));
    clust_num = clust_nums(target_blocks(bb));
    
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',target_blocks(bb))];
    load(cur_dat_name,'Clusters');
    
    spike_xy = Clusters{probe_num}.spike_xy;
    N_spks = size(spike_xy,1);
    spike_labels = Clusters{probe_num}.spike_clusts;
    N_sus = length(unique(spike_labels(spike_labels > 1)));
    spike_labels(~ismember(spike_labels,clust_num)) = 1;
    clust_spikes = find(spike_labels == clust_num);
    
    mu_inds = find(spike_labels == 1);
    out_inds = find(spike_labels == -1);
    
    subplot(n_cols,n_rows,bb);hold off
    plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.'); hold on
    plot(spike_xy(out_inds,1),spike_xy(out_inds,2),'r.');
    plot(spike_xy(clust_spikes,1),spike_xy(clust_spikes,2),'.');
    set(gca,'xtick',[],'ytick',[]);axis tight
    for ii = 1:length(Clusters{probe_num}.cluster_labels)
        if ~isnan(Clusters{probe_num}.gmm_xyMeans(ii,1))
        h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_xyMeans(ii,:)',squeeze(Clusters{probe_num}.gmm_xySigma(:,:,ii)),[2],'r',1);
        end
    end
    if target_blocks(bb) == ref_block
        title(['Block #',int2str(target_blocks(bb))],'Color','r');
    else
        title(['Block #',int2str(target_blocks(bb))],'Color','k');
    end
    
end

probe_num = probe_nums(ref_block);
%scale axes
set(0,'CurrentFigure',full_scatter_fig);
base_block_ind = find(target_blocks == ref_block);
if ~isempty(base_block_ind)
    subplot(n_cols,n_rows,base_block_ind);
    xl = xlim(); yl = ylim();
else
    xl = minmax(RefClusters{probe_num}.spike_xy(:,1));
    yl = minmax(RefClusters{probe_num}.spike_xy(:,2));
end
for bb = 1:length(target_blocks)
    subplot(n_cols,n_rows,bb);
    xlim(xl); ylim(yl);
end

