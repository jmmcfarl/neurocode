function [] = regenerate_allblock_xydensity(probe_num,target_blocks)

global base_save_dir full_save_dir 

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name,'RefClusters');
cur_base_block = RefClusters{probe_num}.base_block;

pfname_de = [full_save_dir sprintf('/Probe%d_fullclust_dens.fig',probe_num)];
full_dense_fig = figure('visible','off');
n_blocks = max(target_blocks);
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);

ref_xy = RefClusters{probe_num}.spike_xy;
ref_spike_ids = RefClusters{probe_num}.spike_clusts;
uids = ref_spike_ids > 0;
subplot(n_cols,n_rows,cur_base_block);
[handles, details] = DensityPlot_jmm(ref_xy(uids,1),ref_xy(uids,2),'sqrtsc','ynormal','sd',[1 1]);
dens_xrange = minmax(details.x);
dens_yrange = minmax(details.y);

for bb = target_blocks
    fprintf('Block %d\n',bb);
    
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',bb)];
    load(cur_dat_name,'Clusters');
    
    spike_xy = Clusters{probe_num}.spike_xy;
    spike_clusts = Clusters{probe_num}.spike_clusts;
    uids = spike_clusts > 0;
    subplot(n_cols,n_rows,bb);hold off
    [handles, details] = DensityPlot_jmm(spike_xy(uids,1),spike_xy(uids,2),'sqrtsc','ynormal','sd',[1 1],'xrange',dens_xrange,'yrange',dens_yrange);
    hold on
    set(gca,'xtick',[],'ytick',[]);axis tight
    for ii = 1:length(Clusters{probe_num}.cluster_labels)
        h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_xyMeans(ii,:)',squeeze(Clusters{probe_num}.gmm_xySigma(:,:,ii)),[2],'r',0.5);
    end
    if bb == cur_base_block
        title(['Block #',int2str(bb)],'Color','r');
    else
        title(['Block #',int2str(bb)],'Color','k');
    end
        
end

%scale axes
set(0,'CurrentFigure',full_dense_fig);
subplot(n_cols,n_rows,cur_base_block);
xl = xlim(); yl = ylim();
for bb = target_blocks
    subplot(n_cols,n_rows,bb);
    xlim(xl); ylim(yl);
end

saveas(full_dense_fig,pfname_de);

fillPage(full_dense_fig,'papersize',[14 14]);
fname = [full_save_dir sprintf('/Probe%d_fullclust_dens',probe_num)];
print(full_dense_fig,fname,'-dpng');
close(full_dense_fig);

