function [] = regenerate_allprobe_xydensity(target_probes,block_num)

global base_save_dir full_save_dir 

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name,'RefClusters');
% cur_base_block = RefClusters{probe_num}.base_block;

pfname_de = [full_save_dir sprintf('/Block%d_fullclust_dens.fig',block_num)];
full_dense_fig = figure('visible','off');
n_probes = max(target_probes);
n_cols = ceil(sqrt(n_probes)); n_rows = ceil(n_probes/n_cols);

% ref_xy = RefClusters{probe_num}.spike_xy;
% ref_spike_ids = RefClusters{probe_num}.spike_clusts;
% uids = ref_spike_ids > 0;
% subplot(n_cols,n_rows,cur_base_block);
% [handles, details] = DensityPlot_jmm(ref_xy(uids,1),ref_xy(uids,2),'sqrtsc','ynormal','sd',[1 1]);
% dens_xrange = minmax(details.x);
% dens_yrange = minmax(details.y);

cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
load(cur_dat_name,'Clusters');

for probe_num = target_probes
    fprintf('Probe %d\n',probe_num);
    
    spike_xy = Clusters{probe_num}.spike_xy;
    spike_clusts = Clusters{probe_num}.spike_clusts;
    uids = spike_clusts > 0;
    if ~isempty(uids)
    subplot(n_cols,n_rows,probe_num);hold off
%     [handles, details] = DensityPlot_jmm(spike_xy(uids,1),spike_xy(uids,2),'sqrtsc','ynormal','sd',[1 1],'xrange',dens_xrange,'yrange',dens_yrange);
    [handles, details] = DensityPlot_jmm(spike_xy(uids,1),spike_xy(uids,2),'sqrtsc','ynormal','sd',[1 1]);
    hold on
    set(gca,'xtick',[],'ytick',[]);axis tight
    for ii = 1:length(Clusters{probe_num}.cluster_labels)
         if ~isnan(Clusters{probe_num}.gmm_xyMeans(ii,1))
       h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_xyMeans(ii,:)',squeeze(Clusters{probe_num}.gmm_xySigma(:,:,ii)),[2],'r',0.5);
         end
    end
%     if bb == cur_base_block
%         title(['Block #',int2str(bb)],'Color','r');
%     else
        title(['Probe #',int2str(probe_num)],'Color','k');
%     end
    end
end

% %scale axes
% set(0,'CurrentFigure',full_dense_fig);
% subplot(n_cols,n_rows,cur_base_block);
% xl = xlim(); yl = ylim();
% for bb = target_blocks
%     subplot(n_cols,n_rows,bb);
%     xlim(xl); ylim(yl);
% end

% saveas(full_dense_fig,pfname_de);

fillPage(full_dense_fig,'papersize',[14 14]);
fname = [full_save_dir sprintf('/Block%d_fullclust_dens',block_num)];
print(full_dense_fig,fname,'-dpng');
close(full_dense_fig);

