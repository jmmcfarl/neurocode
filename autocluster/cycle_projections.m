function [] = cycle_projections(block_num,probe_num)

global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData raw_block_nums

fprintf('Loading block %d Clusters\n',block_num);
cur_clust_data = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
load(cur_clust_data,'Clusters');

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
params = Clusters{probe_num}.params;

[~,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,Clusters{probe_num},Clusters{probe_num}.params,1);

cluster_labels = Clusters{probe_num}.cluster_labels;
%%
[N_spks,N_samps,N_chs] = size(Spikes.V);
Fs = 3e4;

spike_labels = Clusters{probe_num}.spike_clusts;
uids = find(Clusters{probe_num}.comp_idx > 0);
mu_inds = find(spike_labels == 1);
out_inds = find(Clusters{probe_num}.comp_idx == -1);
N_sus = length(unique(spike_labels(uids))) - 1;
for ii = 1:N_sus
    su_inds{ii} = find(spike_labels == ii + 1);
end
cmap = cluster_cmap(N_sus);

figure
% DENSITY PLOT IN XY SPACE
subplot(2,1,1);
[handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
hold on
for ii = 1:length(cluster_labels)
    if cluster_labels(ii) == 1
        h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_xyMeans(ii,:)',squeeze(Clusters{probe_num}.gmm_xySigma(:,:,ii)),[2],'k',2);
    else
        h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_xyMeans(ii,:)',squeeze(Clusters{probe_num}.gmm_xySigma(:,:,ii)),[2],cmap(cluster_labels(ii)-1,:),2);
    end
end
set(gca,'ytick',[]);
if Clusters{probe_num}.template_it == -1
    title('Used PCs','fontsize',12);
elseif Clusters{probe_num}.template_it == -2
    title('Used voltage proj','fontsize',12);
elseif Clusters{probe_num}.template_it == -3
    title('Used energy proj','fontsize',12);
else
    title('Used template projection','fontsize',12);
end
xrange = minmax(details.x);
yrange = minmax(details.y);

% SCATTERPLOT IN XY SPACE
subplot(2,1,2); hold on
plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.','markersize',4);
for ii = 1:N_sus
    plot(spike_xy(spike_labels==ii+1,1),spike_xy(spike_labels==ii+1,2),'.','color',cmap(ii,:),'markersize',8);
end
plot(spike_xy(out_inds,1),spike_xy(out_inds,2),'r*','markersize',2);
for ii = 1:length(Clusters{probe_num}.cluster_labels)
    h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_xyMeans(ii,:)',squeeze(Clusters{probe_num}.gmm_xySigma(:,:,ii)),[2],'r',2);
end
xlim(xrange); ylim(yrange);

pause

%%
N_features = size(spike_features,2);
N_projs = N_features*(N_features-1)/2;
for kk = 1:N_features-1
for jj = (kk+1):N_features
    fprintf('Feature %d vs Feature %d\n',kk,jj);
    clf
    
    % DENSITY PLOT IN XY SPACE
    subplot(2,1,1);
    [handles, details] = DensityPlot_jmm(spike_features(:,kk),spike_features(:,jj),'sqrtsc','ynormal','sd',[1 1]);
    hold on
    for ii = 1:length(cluster_labels)
        if cluster_labels(ii) == 1
            h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_fit.mu(ii,[kk jj])',squeeze(Clusters{probe_num}.gmm_fit.Sigma([kk jj],[kk jj],ii)),[2],'k',2);
        else
            h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_fit.mu(ii,[kk jj])',squeeze(Clusters{probe_num}.gmm_fit.Sigma([kk jj],[kk jj],ii)),[2],cmap(cluster_labels(ii)-1,:),2);
        end
    end
    set(gca,'ytick',[]);
    if Clusters{probe_num}.template_it == -1
        title('Used PCs','fontsize',12);
    elseif Clusters{probe_num}.template_it == -2
        title('Used voltage proj','fontsize',12);
    elseif Clusters{probe_num}.template_it == -3
        title('Used energy proj','fontsize',12);
    else
        title('Used template projection','fontsize',12);
    end
    xrange = minmax(details.x);
    yrange = minmax(details.y);
    
    % SCATTERPLOT IN XY SPACE
    subplot(2,1,2); hold on
    plot(spike_features(mu_inds,kk),spike_features(mu_inds,jj),'k.','markersize',4);
    for ii = 1:N_sus
        plot(spike_features(spike_labels==ii+1,kk),spike_features(spike_labels==ii+1,jj),'.','color',cmap(ii,:),'markersize',8);
    end
    plot(spike_features(out_inds,kk),spike_features(out_inds,jj),'r*','markersize',2);
    for ii = 1:length(Clusters{probe_num}.cluster_labels)
        h1 = plot_gaussian_2d(Clusters{probe_num}.gmm_fit.mu(ii,[kk jj])',squeeze(Clusters{probe_num}.gmm_fit.Sigma([kk jj],[kk jj],ii)),[2],'r',2);
    end
    xlim(xrange); ylim(yrange);
    
    pause
    
end
end