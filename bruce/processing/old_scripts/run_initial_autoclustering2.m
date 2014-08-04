clear all
close all
addpath('~/James_scripts/autocluster/');
Expt_name = 'M270';

%%
% data_loc = '/media/NTlab_data1/Data/bruce/';
data_loc = '/home/james/Data/bruce/';
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
    load(sprintf('lem%sExpts.mat',Expt_name));
    n_probes = 24;
    
end
n_blocks = length(Expts);
block_durs = nan(n_blocks,1);
for ii = 1:n_blocks
    if ~isempty(Expts{ii})
        expt_durs(ii) = (Expts{ii}.Header.End - Expts{ii}.Header.Start)/1e4;
    end
end

%choose an initial clustering block that is near the middle and is long
look_range = round(n_blocks/3):round(2*n_blocks/3); %look in the middle third of blocks
[~,longest] = max(expt_durs(look_range));
base_block = look_range(longest);
% base_block = round(n_blocks/2);
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
clust_params.use_ms = 0;
clust_params.target_rate = 'median';

Vloaded = nan;

all_probe_fig = figure('visible','off');
n_cols = ceil(sqrt(n_probes)); n_rows = ceil(n_probes/n_cols);
for probe_num = 1:n_probes
% for probe_num = 22
    fprintf('\nClustering probe %d of %d\n',probe_num,n_probes);
    
    if Expt_name(1) == 'G'
        sfile_name = [data_dir sprintf('/Expt%d.p%dFullV.mat',base_block,probe_num)];
        [cluster_details,spike_features,sum_fig] =detect_and_cluster_2comp(sfile_name,clust_params);
    else
        sfile_name = [data_dir sprintf('/Expt%dFullV.mat',base_block)];
        use_chs = [probe_num-1 probe_num probe_num + 1];
        use_chs(use_chs < 1 | use_chs > n_probes) = [];
        if Vloaded ~= base_block
            fprintf('Loading data file %s\n',sfile_name);
            [sfile.V,sfile.Vtime,sfile.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
            Vloaded = base_block;
        end
        [cluster_details,spike_features,sum_fig] = detect_and_cluster_2comp(sfile,clust_params,use_chs);
    end
    
    
    cluster_details.base_block = base_block;
    Clusters{probe_num} = cluster_details;
    
    fprintf('Saving cluster details\n');
    save(full_dat_name,'Clusters');
    
    spike_xy = spike_features*cluster_details.xy_projmat;
    N_spks = size(spike_xy,1);
    su_inds = find(ismember(cluster_details.comp_idx,find(cluster_details.cluster_labels == 2)));
    mu_inds = setdiff(1:N_spks,su_inds);
    
    %save unit cluster plot
    fillPage(sum_fig,'papersize',[14 8]);
    if cluster_details.dprime > 2 && cluster_details.n_spks(2)/cluster_details.recDur >= 0.25
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
plot(cellfun(@(x) x.Lratios,Clusters),'o-')
ylabel('Lratio','fontsize',12);

%% DO ANY RETRIGGERING
retrig_ch = 22;
retrig_rate = 250;
retrig_sign = -1;
reapply = 0;
block_num = base_block;

cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',block_num)];
load(cur_dat_name);
fprintf('Retriggering probe %d\n',retrig_ch);

if Expt_name(1) == 'G'
    sfile = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,retrig_ch)];
else
    sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
    if Vloaded ~= block_num
        fprintf('Loading data file %s\n',sfile_name);
        [sfile.V,sfile.Vtime,sfile.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
        Vloaded = block_num;
    end
end

cur_clust_params = Clusters{retrig_ch}.params;
cur_clust_params.target_rate = retrig_rate;
cur_clust_params.thresh_sign = retrig_sign;
% cur_clust_params.min_Pcomp = 0.0001;
cur_clust_params.summary_plot = 2; %make summary plot visible

if reapply == 0
    use_chs = [retrig_ch-1 retrig_ch retrig_ch + 1];
    use_chs(use_chs < 1 | use_chs > n_probes) = [];
    [cluster_details,spike_features,sum_fig] = detect_and_cluster_2comp(sfile,cur_clust_params,use_chs);
    cluster_details.base_block = block_num;
else
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile,Clusters{retrig_ch},cur_clust_params);
    sum_fig = create_summary_cluster_fig(cluster_details,Spikes,spike_xy,cur_clust_params);
    clear Spikes
end

resp = input('Keep new clustering (y/n)?','s');
if strcmpi(resp,'Y')
    fprintf('Saving cluster details\n');
    Clusters{retrig_ch} = cluster_details;
    Clusters{retrig_ch}.params = cur_clust_params;
    save(cur_dat_name,'Clusters');
    
    fillPage(gcf,'papersize',[14 8]);
    pname = [good_save_dir sprintf('/Probe%d_Block%d_initclust',retrig_ch,block_num)];
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

probe_num = 22;
block_num = base_block;

cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',block_num)];
load(cur_dat_name);

if Expt_name(1) == 'G'
    sfile = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,retrig_ch)];
else
    sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
    if Vloaded ~= block_num
        fprintf('Loading data file %s\n',sfile_name);
        [sfile.V,sfile.Vtime,sfile.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
        Vloaded = block_num;
    end
end

[new_cluster,used_new] = split_GMM_cluster(sfile,Clusters{probe_num});
if used_new == 1
    Clusters{probe_num} = new_cluster;
    fprintf('Saving cluster details\n');
    save(cur_dat_name,'Clusters');
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = [good_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        print(pname,'-dpng');
        close(sum_fig);
    end
end

%% MERGE COMPONENTS
close all

probe_num = 22;
block_num = base_block;
% block_num = 12;

cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',block_num)];
load(cur_dat_name);

if Expt_name(1) == 'G'
    sfile = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,retrig_ch)];
else
    sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
    if Vloaded ~= block_num
        fprintf('Loading data file %s\n',sfile_name);
        [sfile.V,sfile.Vtime,sfile.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
        Vloaded = block_num;
    end
end
[new_cluster,used_new] = merge_GMM_cluster(sfile,Clusters{probe_num});
if used_new == 1
    Clusters{probe_num} = new_cluster;
    fprintf('Saving cluster details\n');
    save(cur_dat_name,'Clusters');
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = [good_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        print(pname,'-dpng');
        close(sum_fig);
    end
end

%% REASSIGN COMPONENTS
close all

probe_num = 22;
block_num = base_block;
% block_num = 12;

cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',block_num)];
load(cur_dat_name);

if Expt_name(1) == 'G'
    sfile = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,retrig_ch)];
else
    sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
    if Vloaded ~= block_num
        fprintf('Loading data file %s\n',sfile_name);
        [sfile.V,sfile.Vtime,sfile.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
        Vloaded = block_num;
    end
end
[new_cluster,used_new] = reassign_GMM_components(sfile,Clusters{probe_num});
if used_new == 1
    Clusters{probe_num} = new_cluster;
    fprintf('Saving cluster details\n');
    save(cur_dat_name,'Clusters');
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = [good_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        print(pname,'-dpng');
        close(sum_fig);
    end
end

%% ADD NEW CLUSTER
close all

probe_num = 22;
block_num = base_block;

cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',block_num)];
load(cur_dat_name);

if Expt_name(1) == 'G'
    sfile = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,retrig_ch)];
else
    sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
    if Vloaded ~= block_num
        fprintf('Loading data file %s\n',sfile_name);
        [sfile.V,sfile.Vtime,sfile.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
        Vloaded = block_num;
    end
end
[new_cluster,used_new] = add_GMM_cluster(sfile,Clusters{probe_num});
if used_new == 1
    Clusters{probe_num} = new_cluster;
    fprintf('Saving cluster details\n');
    save(cur_dat_name,'Clusters');
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = [good_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        print(pname,'-dpng');
        close(sum_fig);
    end
end

%% RELABEL CLUSTERS
close all

probe_num = 44;
target_blocks = [22 23 24 25 28];
target_labels = [1 3 2];

%look for existing full scatter figure and open if it exists
if ismember(probe_num,good_probes)
    pfname = [good_full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
else
    pfname = [junk_full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
end
if exist(pfname,'file') && ~ishandle(full_scatter_fig)
    open(pfname); full_scatter_fig = gcf(); set(full_scatter_fig,'visible','off');
else
    fprintf('Need to generate full cluster scatter first\n');
end
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);
subplot_handles = get(full_scatter_fig,'Children');
set(0,'CurrentFigure',full_scatter_fig);
xl = xlim(); yl = ylim();

for bb = target_blocks
    
    cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',bb)];
    if exist(cur_dat_name,'file')
        fprintf('Loading clusters for block %d\n',bb);
        load(cur_dat_name);
    else
        fprintf('Initiating new Cluster array for block %d\n',bb);
        Clusters = cell(n_probes,1);
    end
    
    fprintf('Probe %d: Applying clusters to block %d\n',probe_num,bb);
    sfile_name = [data_dir sprintf('/Expt%d.p%dFullV.mat',bb,probe_num)];
    [new_Cluster,spike_features,spike_xy] = apply_clustering(sfile_name,Clusters{probe_num},[],1);
    N_spks = size(spike_xy,1);
    
    %make switch of cluster labels
    cur_labels = new_Cluster.cluster_labels;
    new_labels = cur_labels;
    for ii = 1:length(target_labels)
        new_labels(cur_labels == ii) = target_labels(ii);
    end
    new_Cluster.cluster_labels = new_labels;
    Clusters{probe_num} = new_Cluster;
    save(cur_dat_name);
    
    spike_labels = zeros(size(new_Cluster.comp_idx));
    uids = find(new_Cluster.comp_idx > 0);
    spike_labels(uids) = new_Cluster.cluster_labels(new_Cluster.comp_idx(uids));
    mu_inds = find(spike_labels == 1);
    
    probe_clust_means(bb,1) = mean(spike_xy(mu_inds,1));
    probe_clust_stds(bb,1) = std(spike_xy(mu_inds,1));
    for ii = 1:N_sus
        probe_clust_means(bb,ii+1) = mean(spike_xy(spike_labels==ii+1,1));
        probe_clust_stds(bb,ii+1) = std(spike_xy(spike_labels==ii+1,1));
    end
    
    subplot(n_cols,n_rows,bb); hold off
    plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.','markersize',1);
    hold on
    for ii = 1:N_sus
        plot(spike_xy(spike_labels==ii+1,1),spike_xy(spike_labels==ii+1,2),'.','color',cmap(ii,:));
    end
    set(gca,'xtick',[],'ytick',[]);
    xlim(xl); ylim(yl);
    for ii = 1:length(new_Cluster.cluster_labels)
        h1 = plot_gaussian_2d(new_Cluster.gmm_xyMeans(ii,:)',squeeze(new_Cluster.gmm_xySigma(:,:,ii)),[1 2],'r',1);
    end
    title(['Block #',int2str(bb)],'Color','k');
    
end
pname = [good_full_save_dir sprintf('/Probe%d_fullclust_scatter',probe_num)];
print(full_scatter_fig,pname,'-dpng');
saveas(full_scatter_fig,pfname);
close(full_scatter_fig);

%% REMOVE CLUSTER
close all

probe_num = 44;
target_blocks = [22 23 24 25 28];
target_cluster = [2];

%look for existing full scatter figure and open if it exists
if ismember(probe_num,good_probes)
    pfname = [good_full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
else
    pfname = [junk_full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
end
if exist(pfname,'file') && ~ishandle(full_scatter_fig)
    open(pfname); full_scatter_fig = gcf(); set(full_scatter_fig,'visible','off');
else
    fprintf('Need to generate full cluster scatter first\n');
end
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);
subplot_handles = get(full_scatter_fig,'Children');
set(0,'CurrentFigure',full_scatter_fig);
xl = xlim(); yl = ylim();

for bb = target_blocks
    
    cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',bb)];
    if exist(cur_dat_name,'file')
        fprintf('Loading clusters for block %d\n',bb);
        load(cur_dat_name);
    else
        fprintf('Initiating new Cluster array for block %d\n',bb);
        Clusters = cell(n_probes,1);
    end
    
    fprintf('Probe %d: Applying clusters to block %d\n',probe_num,bb);
    sfile_name = [data_dir sprintf('/Expt%d.p%dFullV.mat',bb,probe_num)];
    [new_Cluster,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,Clusters{probe_num},[],1);
    N_spks = size(spike_xy,1);
    
    gauss_means = new_Cluster.gmm_fit.mu;
    gauss_Sigmas = new_Cluster.gmm_fit.Sigma;
    bad_components = find(new_Cluster.cluster_labels ==target_cluster);
    gauss_means(bad_components,:) = [];
    gauss_Sigmas(:,:,bad_components) = [];
    S.mu = gauss_means;
    S.Sigma = gauss_Sigmas;
    init_labels = new_Cluster.cluster_labels;
    init_labels(bad_components) = [];
    [gmm_obj, gmm_distance, clust_ids, cluster_labels, cluster_stats, outliers] = ...
        GMM_fit(Spikes.V, spike_features, size(S.mu,1), new_Cluster.params, [],S,init_labels);
    new_Cluster.gmm_fit = gmm_obj;
    new_Cluster.cluster_labels = init_labels;
    new_Cluster.comp_idx = clust_ids;
    Clusters{probe_num} = new_Cluster;
    save(cur_dat_name);
    
    spike_labels = zeros(size(new_Cluster.comp_idx));
    uids = find(new_Cluster.comp_idx > 0);
    spike_labels(uids) = new_Cluster.cluster_labels(new_Cluster.comp_idx(uids));
    mu_inds = find(spike_labels == 1);
    
    probe_clust_means(bb,1) = mean(spike_xy(mu_inds,1));
    probe_clust_stds(bb,1) = std(spike_xy(mu_inds,1));
    for ii = 1:N_sus
        probe_clust_means(bb,ii+1) = mean(spike_xy(spike_labels==ii+1,1));
        probe_clust_stds(bb,ii+1) = std(spike_xy(spike_labels==ii+1,1));
    end
    
    subplot(n_cols,n_rows,bb); hold off
    plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.','markersize',1);
    hold on
    for ii = 1:N_sus
        plot(spike_xy(spike_labels==ii+1,1),spike_xy(spike_labels==ii+1,2),'.','color',cmap(ii,:));
    end
    set(gca,'xtick',[],'ytick',[]);
    xlim(xl); ylim(yl);
    for ii = 1:length(new_Cluster.cluster_labels)
        h1 = plot_gaussian_2d(new_Cluster.gmm_xyMeans(ii,:)',squeeze(new_Cluster.gmm_xySigma(:,:,ii)),[1 2],'r',1);
    end
    title(['Block #',int2str(bb)],'Color','k');
    
end
pname = [good_full_save_dir sprintf('/Probe%d_fullclust_scatter',probe_num)];
print(full_scatter_fig,pname,'-dpng');
saveas(full_scatter_fig,pfname);
close(full_scatter_fig);

%% NOW APPLY THIS CLUSTERING TO ALL EXPTS FOR ALL PROBES
junk_full_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/junk_full'];
good_full_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/good_full'];
if ~exist(junk_full_save_dir,'dir');
    mkdir(junk_full_save_dir);
end
if ~exist(good_full_save_dir,'dir');
    mkdir(good_full_save_dir);
end


n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);

target_blocks = 1:n_blocks;

for probe_num = 1:n_probes
    % for probe_num = 44
    
    close all
    
    full_dens_fig = figure('visible','off');
    %look for existing full scatter figure and open if it exists
    if ismember(probe_num,good_probes)
        pfname = [good_full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
    else
        pfname = [junk_full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
    end
    %     if exist(pfname,'file') && ~ishandle(full_scatter_fig)
    %         open(pfname); full_scatter_fig = gcf(); set(full_scatter_fig,'visible','off');
    %     else
    full_scatter_fig = figure('visible','off');
    %     end
    
    fprintf('Loading clusters for block %d\n',base_block);
    load(full_dat_name);
    
    base_Cluster = Clusters{probe_num};
    if Expt_name(1) == 'G'
        sfile_name = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,probe_num)];
    else
        sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
    end
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,base_Cluster,[],1);
    N_spks = size(spike_xy,1);
    spike_labels = zeros(size(cluster_details.comp_idx));
    uids = find(cluster_details.comp_idx > 0);
    spike_labels(uids) = cluster_details.cluster_labels(cluster_details.comp_idx(uids));
    mu_inds = find(spike_labels == 1);
    N_sus = length(unique(spike_labels(uids))) - 1;
    cmap = jet(N_sus);
    hold on
    set(0,'CurrentFigure',full_scatter_fig);
    subplot(n_cols,n_rows,base_block);hold on
    plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.');
    for ii = 1:N_sus
        plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
    end
    set(gca,'xtick',[],'ytick',[]);axis tight
    xl = xlim(); yl = ylim();
    for ii = 1:length(base_Cluster.cluster_labels)
        %         if base_Cluster.cluster_labels(ii) == 1
        h1 = plot_gaussian_2d(base_Cluster.gmm_xyMeans(ii,:)',squeeze(base_Cluster.gmm_xySigma(:,:,ii)),[1 2],'r',1);
        %         else
        %             h1 = plot_gaussian_2d(base_Cluster.gmm_xyMeans(ii,:)',squeeze(base_Cluster.gmm_xySigma(:,:,ii)),[1 2],'b',1);
        %         end
    end
    title(['Block #',int2str(base_block)],'Color','r');
    
    set(0,'CurrentFigure',full_dens_fig);
    subplot(n_cols,n_rows,base_block);hold on
    [handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
    dens_xrange = minmax(details.x);
    dens_yrange = minmax(details.y);
    set(gca,'xtick',[],'ytick',[]);
    xlim(dens_xrange); ylim(dens_yrange);
    title(['Block #',int2str(base_block)],'Color','r');
    
    probe_iso_quality(base_block,probe_num,:) = [base_Cluster.dprime base_Cluster.LL base_Cluster.Lratios(1) base_Cluster.iso_dists(1)];
    probe_clust_means = nan(n_blocks,N_sus+1);
    probe_clust_stds = nan(n_blocks,N_sus+1);
    probe_clust_means(base_block,1) = mean(spike_xy(mu_inds,1));
    probe_clust_stds(base_block,1) = std(spike_xy(mu_inds,1));
    for ii = 1:N_sus
        probe_clust_means(base_block,ii+1) = mean(spike_xy(spike_labels==ii+1,1));
        probe_clust_stds(base_block,ii+1) = std(spike_xy(spike_labels==ii+1,1));
    end
    
    for bb = target_blocks
        cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',bb)];
        if exist(cur_dat_name,'file')
            fprintf('Loading clusters for block %d\n',bb);
            load(cur_dat_name);
        else
            fprintf('Initiating new Cluster array for block %d\n',bb);
            Clusters = cell(n_probes,1);
        end
        
        fprintf('Probe %d: Applying clusters to block %d of %d\n',probe_num,bb,n_blocks);
        if Expt_name(1) == 'G'
            sfile_name = [data_dir sprintf('/Expt%d.p%dFullV.mat',bb,probe_num)];
        else
            sfile_name = [data_dir sprintf('/Expt%dFullV.mat',bb)];
        end
        if exist(sfile_name,'file')
            [new_Cluster,spike_features,spike_xy] = apply_clustering(sfile_name,base_Cluster);
            N_spks = size(spike_xy,1);
            
            set(0,'CurrentFigure',full_dens_fig);
            subplot(n_cols,n_rows,bb); hold on
            [handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1],'xrange',dens_xrange,'yrange',dens_yrange);
            set(gca,'xtick',[],'ytick',[]);
            xlim(dens_xrange); ylim(dens_yrange);
            title(['Block #',int2str(bb)],'Color','k');
            
            if new_Cluster.failed == 1
                probe_clust_means(bb,:) = nan;
                probe_clust_stds(bb,:) = nan;
                
                set(0,'CurrentFigure',full_scatter_fig);
                subplot(n_cols,n_rows,bb); hold on
                plot(spike_xy(:,1),spike_xy(:,2),'k.','markersize',1);
                set(gca,'xtick',[],'ytick',[]);
                xlim(xl); ylim(yl);
                title(['Block #',int2str(bb)],'Color','k');
                
                probe_iso_quality(bb,probe_num,:) = [nan nan nan nan];
            else
                spike_labels = zeros(size(new_Cluster.comp_idx));
                uids = find(new_Cluster.comp_idx > 0);
                spike_labels(uids) = new_Cluster.cluster_labels(new_Cluster.comp_idx(uids));
                mu_inds = find(spike_labels == 1);
                
                probe_clust_means(bb,1) = mean(spike_xy(mu_inds,1));
                probe_clust_stds(bb,1) = std(spike_xy(mu_inds,1));
                for ii = 1:N_sus
                    probe_clust_means(bb,ii+1) = mean(spike_xy(spike_labels==ii+1,1));
                    probe_clust_stds(bb,ii+1) = std(spike_xy(spike_labels==ii+1,1));
                end
                
                set(0,'CurrentFigure',full_scatter_fig);
                subplot(n_cols,n_rows,bb); hold on
                plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.','markersize',1);
                for ii = 1:N_sus
                    plot(spike_xy(spike_labels==ii+1,1),spike_xy(spike_labels==ii+1,2),'.','color',cmap(ii,:));
                end
                set(gca,'xtick',[],'ytick',[]);
                xlim(xl); ylim(yl);
                for ii = 1:length(new_Cluster.cluster_labels)
                    %                 if new_Cluster.cluster_labels(ii) == 1
                    h1 = plot_gaussian_2d(new_Cluster.gmm_xyMeans(ii,:)',squeeze(new_Cluster.gmm_xySigma(:,:,ii)),[1 2],'r',1);
                    %                 else
                    %                     h1 = plot_gaussian_2d(new_Cluster.gmm_xyMeans(ii,:)',squeeze(new_Cluster.gmm_xySigma(:,:,ii)),[1 2],'b',1);
                    %                 end
                end
                title(['Block #',int2str(bb)],'Color','k');
                
                probe_iso_quality(bb,probe_num,:) = [new_Cluster.dprime new_Cluster.LL new_Cluster.Lratio new_Cluster.iso_dist];
            end
            fprintf('Saving clusters for block %d\n',bb);
            Clusters{probe_num} = new_Cluster;
            save(cur_dat_name,'Clusters');
        else
            probe_iso_quality(bb,probe_num,:) = nan;
        end
    end
    
    fillPage(full_scatter_fig,'papersize',[15 15]);
    if ismember(probe_num,good_probes)
        pname = [good_full_save_dir sprintf('/Probe%d_fullclust_scatter',probe_num)];
        pfname = [good_full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
    else
        pname = [junk_full_save_dir sprintf('/Probe%d_fullclust_scatter',probe_num)];
        pfname = [junk_full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
    end
    print(full_scatter_fig,pname,'-dpng');
    saveas(full_scatter_fig,pfname)
    close(full_scatter_fig);
    
    fillPage(full_dens_fig,'papersize',[15 15]);
    if ismember(probe_num,good_probes)
        pname = [good_full_save_dir sprintf('/Probe%d_fullclust_dens',probe_num)];
    else
        pname = [junk_full_save_dir sprintf('/Probe%d_fullclust_dens',probe_num)];
    end
    print(full_dens_fig,pname,'-dpng');
    close(full_dens_fig);
    
    f3 = figure('visible','off');
    subplot(3,1,1)
    ax=plotyy(1:n_blocks,squeeze(probe_iso_quality(:,probe_num,1)),1:n_blocks,squeeze(probe_iso_quality(:,probe_num,4)));
    xlim(ax(1),[0 n_blocks+1]); xlim(ax(2),[0 n_blocks+1]);
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
        pname = [good_full_save_dir sprintf('/Probe%d_fullclust_quality',probe_num)];
    else
        pname = [junk_full_save_dir sprintf('/Probe%d_fullclust_quality',probe_num)];
    end
    print(f3,pname,'-dpng');
    close(f3);
    
end


