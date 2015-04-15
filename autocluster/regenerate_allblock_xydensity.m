function [] = regenerate_allblock_xydensity(probe_num,target_blocks)

global base_save_dir full_save_dir

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name,'RefClusters');
cur_base_block = RefClusters{probe_num}.base_block;

pfname_de = [full_save_dir sprintf('/Probe%d_fullclust_dens.fig',probe_num)];
full_dense_fig = figure('visible','off');
n_blocks = length(target_blocks);
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);

ref_xy = RefClusters{probe_num}.spike_xy;
ref_spike_ids = RefClusters{probe_num}.spike_clusts;
uids = ref_spike_ids > 0;
base_block_loc = find(target_blocks == cur_base_block);
if ~isempty(base_block_loc)
    subplot(n_cols,n_rows,base_block_loc);
    [handles, details] = DensityPlot_jmm(ref_xy(uids,1),ref_xy(uids,2),'sqrtsc','ynormal','sd',[1 1]);
    dens_xrange = minmax(details.x);
    dens_yrange = minmax(details.y);
else
    dens_xrange = minmax(ref_xy(uids,1));
    dens_yrange = minmax(ref_xy(uids,2));
end

for bb = 1:length(target_blocks)
    fprintf('Block %d\n',target_blocks(bb));
    
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',target_blocks(bb))];
    load(cur_dat_name,'Clusters');
    
    spike_xy = Clusters{probe_num}.spike_xy;
    spike_clusts = Clusters{probe_num}.spike_clusts;
    uids = spike_clusts > 0;
    subplot(n_cols,n_rows,bb);hold off
    [handles, details] = DensityPlot_jmm(spike_xy(uids,1),spike_xy(uids,2),'sqrtsc','ynormal','sd',[1 1],'xrange',dens_xrange,'yrange',dens_yrange);
    hold on
    set(gca,'xtick',[],'ytick',[]);axis tight
    for ii = 1:length(Clusters{probe_num}.cluster_labels)
        if ~isnan(Clusters{probe_num}.gmm_xyMeans(ii,1))
            h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_xyMeans(ii,:)',squeeze(Clusters{probe_num}.gmm_xySigma(:,:,ii)),[2],'r',0.5);
        end
    end
    if target_blocks(bb) == cur_base_block
        title(['Block #',int2str(target_blocks(bb))],'Color','r');
    else
        title(['Block #',int2str(target_blocks(bb))],'Color','k');
    end
    
end

%scale axes
set(0,'CurrentFigure',full_dense_fig);
if ~isempty(base_block_loc)
    subplot(n_cols,n_rows,base_block_loc);
    xl = xlim(); yl = ylim();
else
    xl = dens_xrange;
    yl = dens_yrange;
end
for bb = 1:length(target_blocks)
    subplot(n_cols,n_rows,bb);
    xlim(xl); ylim(yl);
end

% saveas(full_dense_fig,pfname_de);

fillPage(full_dense_fig,'papersize',[14 14]);
fname = [full_save_dir sprintf('/Probe%d_fullclust_dens',probe_num)];
print(full_dense_fig,fname,'-dpng');
close(full_dense_fig);

