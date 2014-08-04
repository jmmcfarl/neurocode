clear all
close all
addpath('~/James_scripts/autocluster/');
Expt_name = 'G095';

%%
data_loc = '/media/NTlab_data1/Data/bruce/';
% data_loc = '/home/james/Data/bruce/';
data_dir = [data_loc Expt_name];
base_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
init_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/init'];

if ~exist(base_save_dir,'dir');
    mkdir(base_save_dir);
end
if ~exist(init_save_dir,'dir');
    mkdir(init_save_dir);
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
Vloaded = nan;

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
% clust_params.target_rate = 'median';
clust_params.target_rate = 50;

all_probe_fig = figure('visible','off');
n_cols = ceil(sqrt(n_probes)); n_rows = ceil(n_probes/n_cols);
for probe_num = 1:n_probes
    %     for probe_num = 10
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
    putative_good(probe_num) = 1;
    pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,base_block)];
    print(sum_fig,pname,'-dpng');
    close(sum_fig);
    
    %add to all probe plot
    set(0,'CurrentFigure',all_probe_fig);
    subplot(n_cols,n_rows,probe_num);hold on
    plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.','markersize',1);
    plot(spike_xy(su_inds,1),spike_xy(su_inds,2),'r.','markersize',1);
    set(gca,'xtick',[],'ytick',[]);axis tight
    for ii = 1:length(cluster_details.cluster_labels)
        %         if cluster_details.cluster_labels(ii) == 1
        %             h1 = plot_gaussian_2d(cluster_details.gmm_xyMeans(ii,:)',squeeze(cluster_details.gmm_xySigma(:,:,ii)),[2],'g',1);
        %         else
        h1 = plot_gaussian_2d(cluster_details.gmm_xyMeans(ii,:)',squeeze(cluster_details.gmm_xySigma(:,:,ii)),[2],'b',1);
        %         end
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
retrig_ch = 76;
retrig_rate = 50;
retrig_sign = -1;
reapply = 0;
block_num = base_block;
% block_num = 4;

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

cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',base_block)];
load(cur_dat_name);

fprintf('Retriggering probe %d\n',retrig_ch);
cur_clust_params = Clusters{retrig_ch}.params;
cur_clust_params.target_rate = retrig_rate;
cur_clust_params.thresh_sign = retrig_sign;
% cur_clust_params.min_Pcomp = 0.0001;
cur_clust_params.summary_plot = 2; %make summary plot visible
cur_clust_params.try_features = [1 2 3 4];
cur_clust_params.spk_pts = [-12:27];

cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',block_num)];
load(cur_dat_name);

if reapply == 0
    if n_probes == 24
        use_chs = [retrig_ch-1 retrig_ch retrig_ch + 1];
        use_chs(use_chs < 1 | use_chs > n_probes) = [];
    else
        use_chs = [];
    end
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
    save(cur_dat_name,'Clusters');
    
    fillPage(gcf,'papersize',[14 8]);
    pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',retrig_ch,block_num)];
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

probe_num = 59;
block_num = base_block;
% block_num = 4;

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
        pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        fillPage(gcf,'papersize',[14 8]);
        print(pname,'-dpng');
        close(sum_fig);
    end
end

%% MERGE COMPONENTS
close all

probe_num = 13;
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
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        fillPage(gcf,'papersize',[14 8]);
        print(pname,'-dpng');
        close(sum_fig);
    end
end

%% REASSIGN COMPONENTS
close all

probe_num = 44;
block_num = base_block;
% block_num = 12;

cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',block_num)];
load(cur_dat_name);

if Expt_name(1) == 'G'
    sfile = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,probe_num)];
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
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        fillPage(gcf,'papersize',[14 8]);
        print(pname,'-dpng');
        close(sum_fig);
    end
end

%% ADD NEW CLUSTER
close all

probe_num = 44;
block_num = base_block;
% block_num = 4;

cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',block_num)];
load(cur_dat_name);

if Expt_name(1) == 'G'
    sfile = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,probe_num)];
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
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(cluster_details,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        fillPage(gcf,'papersize',[14 8]);
        print(pname,'-dpng');
        close(sum_fig);
    end
end
%% REASSIGN COMPONENTS
close all

probe_num = 24;
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
full_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/full'];
if ~exist(full_save_dir,'dir');
    mkdir(full_save_dir);
end

n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);

target_blocks = [base_block setdiff(1:n_blocks,base_block)];
% target_blocks(target_blocks == 16) = [];

target_probes = [7 9 18 33 44 49 55 58 59 74 76 89];

close all

fprintf('Loading clusters for block %d\n',base_block);
full_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',base_block)];
load(full_dat_name);
base_Clusters = Clusters;

all_clust_means = cell(n_probes,1);
all_clust_stds = cell(n_probes,1);
for bb = target_blocks
    if Expt_name(1) == 'M'
        sfile_name = [data_dir sprintf('/Expt%dFullV.mat',bb)];
        if Vloaded ~= bb
            fprintf('Loading data file %s\n',sfile_name);
            [sfile.V,sfile.Vtime,sfile.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
            Vloaded = bb;
        end
    end
    cur_dat_name = [base_save_dir sprintf('/Block%d_clusters.mat',bb)];
    if exist(cur_dat_name,'file')
        fprintf('Loading clusters for block %d\n',bb);
        load(cur_dat_name);
    else
        fprintf('Initializing clusters for block %d\n',bb);
        Clusters = cell(n_probes,1);
    end
    for probe_num = target_probes
        fprintf('Applying clustering for probe %d\n',probe_num);
        if Expt_name(1) == 'G'
            sfile = [data_dir sprintf('/Expt%d.p%dFullV.mat',bb,probe_num)];
        end
        %look for existing full scatter figure and open if it exists
        pfname_sc = [full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
        if exist(pfname_sc,'file') && ~ishandle(full_scatter_fig)
            full_scatter_fig = open(pfname_sc); set(full_scatter_fig,'visible','off');
        else
            full_scatter_fig = figure('visible','off');
            %                 full_scatter_fig = figure();
        end
        %look for existing full scatter figure and open if it exists
        pfname_de = [full_save_dir sprintf('/Probe%d_fullclust_dens.fig',probe_num)];
        if exist(pfname_de,'file') && ~ishandle(full_dens_fig)
            full_dens_fig = open(pfname_de); set(full_dens_fig,'visible','off');
        else
            full_dens_fig = figure('visible','off');
            %                 full_dens_fig = figure();
        end
        
        if bb == base_block
            [new_cluster,spike_features,spike_xy,Spikes] = apply_clustering(sfile,base_Clusters{probe_num},[],1);
        else
            [new_cluster,spike_features,spike_xy,Spikes] = apply_clustering(sfile,base_Clusters{probe_num});
        end
        
        if ~new_cluster.failed
            N_spks = size(spike_xy,1);
            spike_labels = zeros(size(new_cluster.comp_idx));
            uids = find(new_cluster.comp_idx > 0);
            spike_labels(uids) = new_cluster.cluster_labels(new_cluster.comp_idx(uids));
            mu_inds = find(spike_labels == 1);
            N_sus = length(unique(spike_labels(uids))) - 1;
            cmap = cluster_cmap(N_sus);
            
            set(0,'CurrentFigure',full_scatter_fig);
            %         figure(full_scatter_fig);
            subplot(n_cols,n_rows,bb);hold off
            plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.'); hold on
            for ii = 1:N_sus
                plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
            end
            set(gca,'xtick',[],'ytick',[]);axis tight
            if bb == base_block
                xl(probe_num,:) = xlim(); yl(probe_num,:) = ylim();
            else
                xlim(xl(probe_num,:)); ylim(yl(probe_num,:));
            end
            for ii = 1:length(new_cluster.cluster_labels)
                h1 = plot_gaussian_2d(new_cluster.gmm_xyMeans(ii,:)',squeeze(new_cluster.gmm_xySigma(:,:,ii)),[2],'r',1);
            end
            if bb == base_block
                title(['Block #',int2str(bb)],'Color','r');
            else
                title(['Block #',int2str(bb)],'Color','k');
            end
            saveas(full_scatter_fig,pfname_sc);
            close(full_scatter_fig);
            
            all_clust_means{probe_num}(bb,1) = mean(spike_xy(mu_inds,1));
            all_clust_stds{probe_num}(bb,1) = std(spike_xy(mu_inds,1));
            for ii = 1:N_sus
                all_clust_means{probe_num}(bb,ii+1) = mean(spike_xy(spike_labels==ii+1,1));
                all_clust_stds{probe_num}(bb,ii+1) = std(spike_xy(spike_labels==ii+1,1));
            end
        else
            close(full_scatter_fig);
        end
        
        set(0,'CurrentFigure',full_dens_fig);
        %         figure(full_dens_fig);
        subplot(n_cols,n_rows,bb);hold on
        [handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
        if bb == base_block
            dens_xrange(probe_num,:) = minmax(details.x);
            dens_yrange(probe_num,:) = minmax(details.y);
        end
        set(gca,'xtick',[],'ytick',[]);
        xlim(dens_xrange(probe_num,:)); ylim(dens_yrange(probe_num,:));
        if bb == base_block
            title(['Block #',int2str(bb)],'Color','r');
        else
            title(['Block #',int2str(bb)],'Color','k');
        end
        saveas(full_dens_fig,pfname_de);
        close(full_dens_fig);
        
        probe_iso_quality(bb,probe_num,:) = [new_cluster.dprime new_cluster.LL new_cluster.Lratios(1) new_cluster.iso_dists(1)];
        
        Clusters{probe_num} = new_cluster;
        
    end
    
    fprintf('Saving clusters for block %d\n',bb);
    save(cur_dat_name,'Clusters');
end

%%
close all
for probe_num = target_probes
    %look for existing full scatter figure and open if it exists
    pfname_sc = [full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
    if exist(pfname_sc,'file') %&& ~ishandle(full_scatter_fig)
        full_scatter_fig = open(pfname_sc); set(full_scatter_fig,'visible','on');
    end
    %look for existing full scatter figure and open if it exists
    pfname_de = [full_save_dir sprintf('/Probe%d_fullclust_dens.fig',probe_num)];
    if exist(pfname_de,'file') %&& ~ishandle(full_dens_fig)
        full_dens_fig = open(pfname_de); set(full_dens_fig,'visible','on');
    end
    
    distFig('r1');
    
    
    pause
    close all
end
