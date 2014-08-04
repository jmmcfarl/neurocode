function [] = autocluster_excludeSU(block_num,probe_num,exclude_SUs,clust_params)

global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData

fprintf('Loading block %d Clusters\n',block_num);
cur_clust_data = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
load(cur_clust_data,'Clusters');

if Expt_name(1) == 'G'
    loadedData = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,probe_num)];
else
    sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
    if Vloaded ~= block_num
        fprintf('Loading data file %s\n',sfile_name);
        [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
        Vloaded = block_num;
    end
end
params = Clusters{probe_num}.params;

if nargin >= 4 
    cur_fields = fieldnames(clust_params);
    for ii = 1:length(cur_fields)
       params = setfield(params,cur_fields{ii},getfield(clust_params,cur_fields{ii}));
    end
end

exclude_SU_inds = Clusters{probe_num}.spk_inds(ismember(Clusters{probe_num}.spike_clusts, exclude_SUs));

[clusterDetails,spike_features,sum_fig] = detect_and_cluster_init_excludespks(loadedData,params,Clusters{probe_num}.use_chs,exclude_SU_inds);

figure(sum_fig);
keep = input('Use new cluster (y/n)?','s');
if strcmpi(keep,'y')
    pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust_exclude%d',probe_num,block_num,exclude_SUs)];
    fillPage(gcf,'papersize',[14 8]);
    print(pname,'-dpng');
    close(sum_fig);
    
    appended_features = spike_features*clusterDetails.xy_projmat(:,1);
    
    [~,orig_spike_features,spike_xy,Spikes] = apply_clustering(loadedData,Clusters{probe_num},Clusters{probe_num}.params,1);
    
    full_spike_features = [orig_spike_features appended_features];
    
    NComps = Clusters{probe_num}.Ncomps;
    [GMM_obj, distance,all_comp_idx,all_clust_labels,cluster_stats] = ...
        GMM_fit(Spikes.V,full_spike_features,NComps,Clusters{probe_num}.params,NComps,Clusters{probe_num}.comp_idx,Clusters{probe_num}.cluster_labels);
    
    newCluster = Clusters{probe_num};
    newCluster.gmm_fit = GMM_obj;
    
    [cluster_labels, cluster_stats] = relabel_clusters(Spikes.V,all_comp_idx,all_clust_labels);
    newCluster.cluster_labels = cluster_labels;
    newCluster.mean_spike = cluster_stats.mean_spike;
    newCluster.std_spike = cluster_stats.std_spike;
    
    [spike_xy,xyproj_mat] = Project_GMMfeatures(full_spike_features, GMM_obj,cluster_labels);
    newCluster.spike_xy = spike_xy;
    newCluster.xy_projmat = xyproj_mat;
    %compute GMM model parameters projected into XY space
    gmm_xyMeans = newCluster.gmm_fit.mu*newCluster.xy_projmat;
    for ii = 1:size(gmm_xyMeans,1)
        gmm_xySigma(:,:,ii) = newCluster.xy_projmat' * squeeze(newCluster.gmm_fit.Sigma(:,:,ii)) * newCluster.xy_projmat;
    end
    newCluster.gmm_xyMeans = gmm_xyMeans;
    newCluster.gmm_xySigma = gmm_xySigma;
    
    sum_fig = create_summary_cluster_fig(newCluster,Spikes,spike_xy,newCluster.params)
    
[~,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,Clusters{probe_num},Clusters{probe_num}.params,1);
    
    new_cluster_labels = clusterDetails.cluster_labels;
    new_comp_idx = Clusters{probe_num}.comp_idx;
    max_new_SUs = max(new_cluster_labels);
    new_refracts = clusterDetails.refract;
    new_isos = clusterDetails.iso_dists;
    new_Lrats = clusterDetails.Lratios;
    new_EX_SUs = repmat({[exclude_SUs]},1,max_new_SUs-1);
    for ii = 1:length(exclude_SUs)
        cur_comps = find(Clusters{probe_num}.cluster_labels==exclude_SUs(ii));
        for jj = 1:length(cur_comps)
            new_comp_idx(new_comp_idx == cur_comps(jj)) = length(new_cluster_labels)+1;
            new_cluster_labels = [new_cluster_labels max_new_SUs + 1];
        end
        new_refracts = [new_refracts; Clusters{probe_num}.refract(exclude_SUs(ii)-1,:)];
        new_isos = [new_isos Clusters{probe_num}.iso_dists(exclude_SUs(ii)-1)];
        new_Lrats = [new_Lrats Clusters{probe_num}.Lratios(exclude_SUs(ii)-1)];
        new_EX_SUs{max_new_SUs} = [];
        max_new_SUs = max_new_SUs + 1;
    end
    uids = new_comp_idx > 0;
    new_spike_clusts = int16(nan(size(new_comp_idx)));
    new_spike_clusts(uids) = (new_cluster_labels(new_comp_idx(uids)));
    new_spike_clusts = new_spike_clusts(:);
    
    newCluster = Clusters{probe_num};
    newCluster.cluster_labels = new_cluster_labels;
    newCluster.comp_idx = new_comp_idx;
    newCluster.spike_clusts = new_spike_clusts;
    newCluster.EX_SUs = new_EX_SUs;
    newCluster.refract = new_refracts;
    newCluster.iso_dists = new_isos;
    newCluster.Lratios = new_Lrats;    
    
    newCluster = compute_cluster_stats(newCluster,Spikes,spike_features);
    
    Clusters{probe_num} = newCluster;
    fprintf('Saving cluster details\n');
    save(cur_clust_data,'Clusters');
    
    if block_num == newCluster.base_block
        sum_fig = create_summary_cluster_fig(newCluster,Spikes,spike_xy,newCluster.params);
        pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        fillPage(gcf,'papersize',[14 8]);
        print(pname,'-dpng');
        
        fprintf('Saving to RefClusters\n');
        rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
        load(rclust_dat_name,'RefClusters');
        RefClusters{probe_num} = newCluster;
        save(rclust_dat_name,'RefClusters');
    else
        resp = input('Save to RefClusters (y/n)?','s');
        if strcmpi(resp,'y')
            fprintf('Saving to RefClusters\n');
            rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
            load(rclust_dat_name,'RefClusters');
            RefClusters{probe_num} = newCluster;
            RefClusters{probe_num}.base_block = block_num;
            save(rclust_dat_name,'RefClusters');
            
            sum_fig = create_summary_cluster_fig(newCluster,Spikes,spike_xy,new_cluster.params);
            fillPage(gcf,'papersize',[14 8]);
            pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
            print(pname,'-dpng');
            close(sum_fig);
        end
    end
end

