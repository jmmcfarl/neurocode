function [full_scatter_fig] = regenerate_allblock_xyscatters(probe_num,target_blocks,target_clusts,to_print)

global base_save_dir full_save_dir 

if nargin < 3 || isempty(target_clusts)
    target_clusts = Inf;
end
if nargin < 4 || isempty(to_print)
    to_print = 1;
end

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name,'RefClusters');
cur_base_block = RefClusters{probe_num}.base_block;

pfname_sc = [full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
full_scatter_fig = figure('visible','off');
n_blocks = max(target_blocks);
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);

for bb = target_blocks
    fprintf('Block %d\n',bb);
    
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',bb)];
    load(cur_dat_name,'Clusters');
    
    spike_xy = Clusters{probe_num}.spike_xy;
    N_spks = size(spike_xy,1);
    spike_labels = Clusters{probe_num}.spike_clusts;
    N_sus = length(unique(spike_labels(spike_labels > 1)));
    if isinf(target_clusts)
        use_SUs = (1:N_sus)+1;
    else
        use_SUs = target_clusts; 
        spike_labels(~ismember(spike_labels,target_clusts)) = 1;
    end
    mu_inds = find(spike_labels == 1);
    out_inds = find(spike_labels == -1);
    cmap = cluster_cmap(length(use_SUs));
    
    subplot(n_cols,n_rows,bb);hold off
    plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.'); hold on
    plot(spike_xy(out_inds,1),spike_xy(out_inds,2),'r.');
    for ii = 1:length(use_SUs)
        plot(spike_xy(spike_labels == use_SUs(ii),1),spike_xy(spike_labels == use_SUs(ii),2),'.','color',cmap(ii,:));
    end
    set(gca,'xtick',[],'ytick',[]);axis tight
    for ii = 1:length(Clusters{probe_num}.cluster_labels)
        h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_xyMeans(ii,:)',squeeze(Clusters{probe_num}.gmm_xySigma(:,:,ii)),[2],'r',1);
    end
    if bb == cur_base_block
        title(['Block #',int2str(bb)],'Color','r');
    else
        title(['Block #',int2str(bb)],'Color','k');
    end
        
end

%scale axes
set(0,'CurrentFigure',full_scatter_fig);
subplot(n_cols,n_rows,cur_base_block);
xl = xlim(); yl = ylim();
for bb = target_blocks
    subplot(n_cols,n_rows,bb);
    xlim(xl); ylim(yl);
end

saveas(full_scatter_fig,pfname_sc);

fillPage(full_scatter_fig,'papersize',[14 14]);
if to_print
    fname = [full_save_dir sprintf('/Probe%d_fullclust_scatter',probe_num)];
    print(full_scatter_fig,fname,'-dpng');
    close(full_scatter_fig);
    full_scatter_fig = nan;
end
