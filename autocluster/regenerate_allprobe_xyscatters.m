function [full_scatter_fig] = regenerate_allprobe_xyscatters(target_probes,block_num)

global base_save_dir full_save_dir 

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name,'RefClusters');
% cur_base_block = RefClusters{probe_num}.base_block;

pfname_sc = [full_save_dir sprintf('/Block%d_fullclust_scatter.fig',block_num)];
full_scatter_fig = figure('visible','off');
n_probes = max(target_probes);
n_cols = ceil(sqrt(n_probes)); n_rows = ceil(n_probes/n_cols);

for probe_num = target_probes
    fprintf('Probe %d\n',probe_num);
    
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
    load(cur_dat_name,'Clusters');
    
    spike_xy = Clusters{probe_num}.spike_xy;
    N_spks = size(spike_xy,1);
    spike_labels = Clusters{probe_num}.spike_clusts;
    N_sus = length(unique(spike_labels(spike_labels > 1)));
    use_SUs = (1:N_sus)+1;
    
    mu_inds = find(spike_labels == 1);
    out_inds = find(spike_labels == -1);
    cmap = cluster_cmap(length(use_SUs));
    
    subplot(n_cols,n_rows,probe_num);hold off
    plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.'); hold on
    plot(spike_xy(out_inds,1),spike_xy(out_inds,2),'r.');
    for ii = 1:length(use_SUs)
        plot(spike_xy(spike_labels == use_SUs(ii),1),spike_xy(spike_labels == use_SUs(ii),2),'.','color',cmap(ii,:));
    end
    set(gca,'xtick',[],'ytick',[]);axis tight
    for ii = 1:length(Clusters{probe_num}.cluster_labels)
                if ~isnan(Clusters{probe_num}.gmm_xyMeans(ii,1))
        h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_xyMeans(ii,:)',squeeze(Clusters{probe_num}.gmm_xySigma(:,:,ii)),[2],'r',1);
                end
    end
%     if bb == cur_base_block
%         title(['Block #',int2str(bb)],'Color','r');
%     else
        title(['Probe #',int2str(probe_num)],'Color','k');
%     end
        
end

% %scale axes
% set(0,'CurrentFigure',full_scatter_fig);
% subplot(n_cols,n_rows,cur_base_block);
% xl = xlim(); yl = ylim();
% for bb = target_blocks
%     subplot(n_cols,n_rows,bb);
%     xlim(xl); ylim(yl);
% end

% saveas(full_scatter_fig,pfname_sc);

fillPage(full_scatter_fig,'papersize',[14 14]);
% if to_print
    fname = [full_save_dir sprintf('/Block%d_fullclust_scatter',block_num)];
    print(full_scatter_fig,fname,'-dpng');
    close(full_scatter_fig);
    full_scatter_fig = nan;
% end
