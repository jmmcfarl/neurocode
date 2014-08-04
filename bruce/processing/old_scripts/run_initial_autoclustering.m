clear all
close all
addpath('~/James_scripts/autocluster/');
Expt_name = 'G095';

%%
data_loc = '/media/NTlab_data1/Data/bruce/';
data_dir = [data_loc Expt_name];
base_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
junk_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/junk_init'];
good_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/good_init'];

if ~exist(base_save_dir,'dir');
    mkdir(base_save_dir);
end
if ~exist(junk_save_dir,'dir');
    mkdir(junk_save_dir);
end
if ~exist(good_save_dir,'dir');
    mkdir(good_save_dir);
end

data_dir2 = ['~/Data/bruce/' Expt_name];

%%
cd(data_dir2);
if Expt_name(1) == 'G';
    load(sprintf('jbe%sExpts.mat',Expt_name));
    n_probes = 96;
elseif Expt_name(1) == 'M'
    
end
n_blocks = length(Expts);
block_durs = nan(n_blocks,1);
for ii = 1:n_blocks
    expt_durs(ii) = (Expts{ii}.Header.End - Expts{ii}.Header.Start)/1e4;
end

%choose an initial clustering block that is near the middle and is long
% look_range = round(n_blocks/3):round(2*n_blocks/3); %look in the middle third of blocks
% [~,longest] = max(expt_durs(look_range));
% base_block = look_range(longest);
base_block = round(n_blocks/2);
fprintf('Initial clustering for all probes for block %d of %d\n',base_block,n_blocks);

%%
full_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',base_block)];
if exist(full_dat_name,'file');
    fprintf('Loading clusters for block %d\n',base_block);
    load(full_dat_name);
else
    fprintf('Initiating new cluster block\n');
    Clusters = cell(n_probes,1);
end

clust_params.gmm_inits = 20;
clust_params.reg_lambda = 0;
clust_params.min_Pcomp = 0.005;
all_probe_fig = figure('visible','off');
n_cols = ceil(sqrt(n_probes)); n_rows = ceil(n_probes/n_cols);
for probe_num = 1:n_probes
% for probe_num = 25
    fprintf('\nClustering probe %d of %d\n',probe_num,n_probes);
    
    sfile_name = [data_dir sprintf('/Expt%d.p%dFullV.mat',base_block,probe_num)];
    
    [cluster_details,spike_features,sum_fig] = detect_and_cluster_2comp(sfile_name,clust_params);
    cluster_details.base_block = base_block;
    Clusters{probe_num} = cluster_details;
    
    fprintf('Saving cluster details\n');
    save(full_dat_name,'Clusters');
    
    spike_xy = spike_features*cluster_details.xy_projmat;
    N_spks = size(spike_xy,1);
    su_inds = find(ismember(cluster_details.cluster_idx,find(cluster_details.cluster_labels == 2)));
    mu_inds = setdiff(1:N_spks,su_inds);
    
    %save unit cluster plot
    fillPage(sum_fig,'papersize',[14 8]);
    if cluster_details.dprime > 2 && cluster_details.su_avg_rate >= 0.25
        putative_good(probe_num) = 1;
        pname = [good_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,base_block)];
    else
        putative_good(probe_num) = 0;
        pname = [junk_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,base_block)];
    end
    print(sum_fig,pname,'-dpng');
    close(sum_fig);
    
    %add to all probe plot
    set(0,'CurrentFigure',all_probe_fig);
    subplot(n_cols,n_rows,probe_num);hold on
    plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.','markersize',1);
    plot(spike_xy(su_inds,1),spike_xy(su_inds,2),'r.','markersize',1);
    set(gca,'xtick',[],'ytick',[]);axis tight
    for ii = 1:length(cluster_details.cluster_labels)
        if cluster_details.cluster_labels(ii) == 1
            h1 = plot_gaussian_2d(cluster_details.gmm_xyMeans(ii,:)',squeeze(cluster_details.gmm_xySigma(:,:,ii)),[1 2],'g',1);
        else
            h1 = plot_gaussian_2d(cluster_details.gmm_xyMeans(ii,:)',squeeze(cluster_details.gmm_xySigma(:,:,ii)),[1 2],'b',1);
        end
    end
    if putative_good(probe_num) == 1
        title(sprintf('P%d',probe_num),'color','r');
    else
        title(sprintf('P%d',probe_num),'color','k');
    end
    
end

good_probes = find(putative_good);

fillPage(all_probe_fig,'papersize',[14 14]);
pname = [base_save_dir '/Allprobe_baseblock_scatter'];
print(all_probe_fig,pname,'-dpng');
close(all_probe_fig);

figure
subplot(2,1,1)
plot(cellfun(@(x) x.dprime,Clusters),'o-')
ylabel('Dprime','fontsize',12);
subplot(2,1,2)
plot(cellfun(@(x) x.Lratio,Clusters),'o-')
ylabel('Lratio','fontsize',12);

%% DO ANY RETRIGGERING
retrig_ch = 18;
retrig_rate = 250;
retrig_sign = -1;
reapply = 1;

load(full_dat_name);

fprintf('Retriggering probe %d\n',retrig_ch);

sfile_name = [data_dir sprintf('/Expt%d.p%dFullV.mat',base_block,retrig_ch)];
cur_clust_params = Clusters{retrig_ch}.params;
cur_clust_params.target_rate = retrig_rate;
cur_clust_params.thresh_sign = retrig_sign;
% cur_clust_params.min_Pcomp = 0.0001;
cur_clust_params.summary_plot = 2; %make summary plot visible

if reapply == 1
    [cluster_details,spike_features,sum_fig] = detect_and_cluster_2comp(sfile_name,cur_clust_params);
else
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,Clusters{retrig_ch},cur_clust_params);
    sum_fig = create_summary_cluster_fig(cluster_details,Spikes,spike_xy,cur_clust_params);
    clear Spikes
end

resp = input('Keep new clustering (y/n)?','s');
if strcmpi(resp,'Y')
    fprintf('Saving cluster details\n');
    Clusters{retrig_ch} = cluster_details;
    Clusters{retrig_ch}.params = cur_clust_params;
    save(full_dat_name,'Clusters');

    fillPage(gcf,'papersize',[14 8]);
    pname = [good_save_dir sprintf('/Probe%d_Block%d_initclust',retrig_ch,base_block)];
    print(pname,'-dpng');
    close(sum_fig);
else
    fprintf('Keeping original clustering\n');
    if ishandle(sum_fig)
    close(sum_fig);
    end
end

%% SPLIT COMPONENTS
close all

probe_num = 18;

full_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',base_block)];
load(full_dat_name);

[new_cluster,used_new] = split_GMM_cluster(sfile_name,Clusters{probe_num});
if used_new == 1
    Clusters{probe_num} = new_cluster;
    fprintf('Saving cluster details\n');
    save(full_dat_name,'Clusters');
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,Clusters{probe_num},[],1);
    sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
    pname = [good_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,base_block)];
    print(pname,'-dpng');
    close(sum_fig);
end
%% NOW APPLY THIS CLUSTERING TO ALL EXPTS FOR ALL PROBES
junk_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/junk_full'];
good_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/good_full'];
if ~exist(junk_save_dir,'dir');
    mkdir(junk_save_dir);
end
if ~exist(good_save_dir,'dir');
    mkdir(good_save_dir);
end 
    set(gca,'xtick',[],'ytick',[]);axis tight
    xl = xlim(); yl = ylim();
    for ii = 1:length(base_Cluster.cluster_labels)
        if base_Cluster.cluster_labels(ii) == 1
            h1 = plot_gaussian_2d(base_Cluster.gmm_xyMeans(ii,:)',squeeze(base_Cluster.gmm_xySigma(:,:,ii)),[1 2],'g',1);
        else
            h1 = plot_gaussian_2d(base_Cluster.gmm_xyMeans(ii,:)',squeeze(base_Cluster.gmm_xySigma(:,:,ii)),[1 2],'b',1);
        end
    end
    title(['Block #',int2str(base_block)],'Color','r');
    
    set(0,'CurrentFigure',f2);
    subplot(n_cols,n_rows,base_block);hold on
    [handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
    dens_xrange = minmax(details.x);
    dens_yrange = minmax(details.y);
    set(gca,'xtick',[],'ytick',[]);
    xlim(dens_xrange); ylim(dens_yrange);
    title(['Block #',int2str(base_block)],'Color','r');
    
    probe_iso_quality(base_block,probe_num,:) = [base_Cluster.dprime base_Cluster.LL base_Cluster.Lratio base_Cluster.iso_dist];
    probe_clust_means(base_block,1) = mean(spike_xy(mu_inds,1));
    probe_clust_means(base_block,2) = mean(spike_xy(su_inds,1));
    probe_clust_stds(base_block,1) = std(spike_xy(mu_inds,1));
    probe_clust_stds(base_block,2) = std(spike_xy(su_inds,1));
    
    for bb = target_blocks
        cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',bb)];
        if exist(cur_dat_name,'file');
            fprintf('Loading clusters for block %d\n',bb);
            load(cur_dat_name);
        else
            fprintf('Initiating new Cluster array for block %d\n',bb);
            Clusters = cell(n_probes,1);
        end
        
        fprintf('Probe %d: Applying clusters to block %d of %d\n',probe_num,bb,n_blocks);
        sfile_name = [data_dir sprintf('/Expt%d.p%dFullV.mat',bb,probe_num)];
        [new_Cluster,spike_features,spike_xy] = apply_clustering(sfile_name,base_Cluster);
        N_spks = size(spike_xy,1);
        
        set(0,'CurrentFigure',f2);
        subplot(n_cols,n_rows,bb); hold on
        [handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1],'xrange',dens_xrange,'yrange',dens_yrange);
        set(gca,'xtick',[],'ytick',[]);
        xlim(dens_xrange); ylim(dens_yrange);
        title(['Block #',int2str(bb)],'Color','k');
       
        if new_Cluster.failed == 1
            probe_clust_means(bb,:) = nan;
            probe_clust_stds(bb,:) = nan;
            
            set(0,'CurrentFigure',f1);
            subplot(n_cols,n_rows,bb); hold on
            plot(spike_xy(:,1),spike_xy(:,2),'k.','markersize',1);
            set(gca,'xtick',[],'ytick',[]);
            xlim(xl); ylim(yl);
            title(['Block #',int2str(bb)],'Color','k');
            
             probe_iso_quality(bb,probe_num,:) = [nan nan nan nan];           
        else
            su_inds = find(ismember(new_Cluster.cluster_idx,find(new_Cluster.cluster_labels == 2)));
            mu_inds = setdiff(1:N_spks,su_inds);
            
            probe_clust_means(bb,1) = mean(spike_xy(mu_inds,1));
            probe_clust_means(bb,2) = mean(spike_xy(su_inds,1));
            probe_clust_stds(bb,1) = std(spike_xy(mu_inds,1));
            probe_clust_stds(bb,2) = std(spike_xy(su_inds,1));
            
            set(0,'CurrentFigure',f1);
            subplot(n_cols,n_rows,bb); hold on
            plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.','markersize',1);
            plot(spike_xy(su_inds,1),spike_xy(su_inds,2),'r.','markersize',1);
            set(gca,'xtick',[],'ytick',[]);
            xlim(xl); ylim(yl);
            for ii = 1:length(new_Cluster.cluster_labels)
                if new_Cluster.cluster_labels(ii) == 1
                    h1 = plot_gaussian_2d(new_Cluster.gmm_xyMeans(ii,:)',squeeze(new_Cluster.gmm_xySigma(:,:,ii)),[1 2],'g',1);
                else
                    h1 = plot_gaussian_2d(new_Cluster.gmm_xyMeans(ii,:)',squeeze(new_Cluster.gmm_xySigma(:,:,ii)),[1 2],'b',1);
                end
            end
            title(['Block #',int2str(bb)],'Color','k');
                        
            probe_iso_quality(bb,probe_num,:) = [new_Cluster.dprime new_Cluster.LL new_Cluster.Lratio new_Cluster.iso_dist];
        end
        fprintf('Saving clusters for block %d\n',bb);
        Clusters{probe_num} = new_Cluster;
        save(cur_dat_name,'Clusters');
    end
    
    fillPage(f1,'papersize',[15 15]);
    if ismember(probe_num,good_probes)
        pname = [good_save_dir sprintf('/Probe%d_fullclust_scatter',probe_num)];
    else
        pname = [junk_save_dir sprintf('/Probe%d_fullclust_scatter',probe_num)];
    end
    print(f1,pname,'-dpng');
    close(f1);
    
    fillPage(f2,'papersize',[15 15]);
    if ismember(probe_num,good_probes)
        pname = [good_save_dir sprintf('/Probe%d_fullclust_dens',probe_num)];
    else
        pname = [junk_save_dir sprintf('/Probe%d_fullclust_dens',probe_num)];
    end
    print(f2,pname,'-dpng');
    close(f2);
    
    f3 = figure('visible','off');
    subplot(3,1,1)
    ax=plotyy(1:n_blocks,squeeze(probe_iso_quality(:,probe_num,1)),1:n_blocks,squeeze(probe_iso_quality(:,probe_num,4)));
    xlim([0 n_blocks+1]);
    xlabel('Block number','fontsize',12);
    set(get(ax(1),'Ylabel'),'String','Dprime');
    set(get(ax(2),'Ylabel'),'String','Iso distance');
    subplot(3,1,2)
    plot(1:n_blocks,squeeze(probe_iso_quality(:,probe_num,3)),'o-');
    xlim([0 n_blocks+1]);
    xlabel('Block number','fontsize',12);
    ylabel('Lratio','fontsize',12);
    subplot(3,1,3)
    errorbar(1:n_blocks,probe_clust_means(:,1),probe_clust_stds(:,1),'ko-');
    hold on
    errorbar(1:n_blocks,probe_clust_means(:,2),probe_clust_stds(:,2),'ro-');
    xlabel('Block number','fontsize',12);
    ylabel('Feature projection','fontsize',12);
    xlim([0 n_blocks+1]);
    fillPage(f3,'papersize',[5 10]);
    
    if ismember(probe_num,good_probes)
        pname = [good_save_dir sprintf('/Probe%d_fullclust_quality',probe_num)];
    else
        pname = [junk_save_dir sprintf('/Probe%d_fullclust_quality',probe_num)];
    end
    print(f3,pname,'-dpng');
    close(f3);
    
end


