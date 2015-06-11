function [] = autocluster_excludeSU(block_num,probe_num,exclude_SUs,clust_params)

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

if nargin >= 4
    cur_fields = fieldnames(clust_params);
    for ii = 1:length(cur_fields)
        params = setfield(params,cur_fields{ii},getfield(clust_params,cur_fields{ii}));
    end
end

exclude_SU_inds = Clusters{probe_num}.spk_inds(ismember(Clusters{probe_num}.spike_clusts, exclude_SUs));

[clusterDetails,spike_features,sum_fig,Spikes] = detect_and_cluster_init_excludespks(loadedData,params,Clusters{probe_num}.use_chs,exclude_SU_inds);

figure(sum_fig);
keep = input('Use these feature (y/n)?','s');
close(sum_fig);
if strcmpi(keep,'y')
    %     pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust_exclude%d',probe_num,block_num,exclude_SUs)];
    %     fillPage(gcf,'papersize',[14 8]);
    %     print(pname,'-dpng');
    %     close(sum_fig);
    
    %     appended_features = spike_features*clusterDetails.xy_projmat(:,1);
    %
    %     [~,orig_spike_features,spike_xy,Spikes] = apply_clustering(loadedData,Clusters{probe_num},Clusters{probe_num}.params,1);
    %
    %     full_spike_features = [orig_spike_features appended_features];
    %
    %     NComps = Clusters{probe_num}.Ncomps;
    %     [GMM_obj, distance,all_comp_idx,all_clust_labels,cluster_stats] = ...
    %         GMM_fit(Spikes.V,full_spike_features,NComps,Clusters{probe_num}.params,NComps,Clusters{probe_num}.comp_idx,Clusters{probe_num}.cluster_labels);
    %
    %     newCluster = Clusters{probe_num};
    %     newCluster.gmm_fit = GMM_obj;
    %
    %     [cluster_labels, cluster_stats] = relabel_clusters(Spikes.V,all_comp_idx,all_clust_labels);
    %     newCluster.cluster_labels = cluster_labels;
    %     newCluster.mean_spike = cluster_stats.mean_spike;
    %     newCluster.std_spike = cluster_stats.std_spike;
    %
    %     [spike_xy,xyproj_mat] = Project_GMMfeatures(full_spike_features, GMM_obj,cluster_labels);
    %     newCluster.spike_xy = spike_xy;
    %     newCluster.xy_projmat = xyproj_mat;
    %     %compute GMM model parameters projected into XY space
    %     gmm_xyMeans = newCluster.gmm_fit.mu*newCluster.xy_projmat;
    %     for ii = 1:size(gmm_xyMeans,1)
    %         gmm_xySigma(:,:,ii) = newCluster.xy_projmat' * squeeze(newCluster.gmm_fit.Sigma(:,:,ii)) * newCluster.xy_projmat;
    %     end
    %     newCluster.gmm_xyMeans = gmm_xyMeans;
    %     newCluster.gmm_xySigma = gmm_xySigma;
    %
    %     sum_fig = create_summary_cluster_fig(newCluster,Spikes,spike_xy,newCluster.params)
    %
    %     [~,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,Clusters{probe_num},Clusters{probe_num}.params,1);
    %
    %     new_cluster_labels = clusterDetails.cluster_labels;
    %     new_comp_idx = Clusters{probe_num}.comp_idx;
    %     max_new_SUs = max(new_cluster_labels);
    %     new_refracts = clusterDetails.refract;
    %     new_isos = clusterDetails.iso_dists;
    %     new_Lrats = clusterDetails.Lratios;
    %     new_EX_SUs = repmat({[exclude_SUs]},1,max_new_SUs-1);
    %     for ii = 1:length(exclude_SUs)
    %         cur_comps = find(Clusters{probe_num}.cluster_labels==exclude_SUs(ii));
    %         for jj = 1:length(cur_comps)
    %             new_comp_idx(new_comp_idx == cur_comps(jj)) = length(new_cluster_labels)+1;
    %             new_cluster_labels = [new_cluster_labels max_new_SUs + 1];
    %         end
    %         new_refracts = [new_refracts; Clusters{probe_num}.refract(exclude_SUs(ii)-1,:)];
    %         new_isos = [new_isos; Clusters{probe_num}.iso_dists(exclude_SUs(ii)-1)];
    %         new_Lrats = [new_Lrats; Clusters{probe_num}.Lratios(exclude_SUs(ii)-1)];
    %         new_EX_SUs{max_new_SUs} = [];
    %         max_new_SUs = max_new_SUs + 1;
    %     end
    %     uids = new_comp_idx > 0;
    %     new_spike_clusts = int16(nan(size(new_comp_idx)));
    %     new_spike_clusts(uids) = (new_cluster_labels(new_comp_idx(uids)));
    %     new_spike_clusts = new_spike_clusts(:);
    %
    %     newCluster = Clusters{probe_num};
    %     newCluster.cluster_labels = new_cluster_labels;
    %     newCluster.comp_idx = new_comp_idx;
    %     newCluster.spike_clusts = new_spike_clusts;
    %     newCluster.EX_SUs = new_EX_SUs;
    %     newCluster.refract = new_refracts;
    %     newCluster.iso_dists = new_isos;
    %     newCluster.Lratios = new_Lrats;
    
    orig_params = Clusters{probe_num}.params;
    params.outlier_thresh = orig_params.outlier_thresh; %use original outlier detection threshold
    [~,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,clusterDetails,params,1);
    
    n_comps = length(clusterDetails.cluster_labels) + length(exclude_SUs);
    
    [GMM_obj, distance,all_comp_idx,all_clust_labels,cluster_stats] = ...
        GMM_fit(Spikes.V,spike_features,n_comps,params,n_comps);
    newCluster = clusterDetails;
    newCluster.params = params;
    newCluster.gmm_fit = GMM_obj;
    newCluster.comp_idx = all_comp_idx;
    newCluster.cluster_labels = all_clust_labels;
    newCluster = compute_cluster_stats(newCluster,Spikes,spike_features);
    
    % all_spike_features =
    [Lratios,iso_distance] = compute_cluster_Lratio(spike_features,newCluster.gmm_fit,newCluster.comp_idx,newCluster.cluster_labels);
    newCluster.Lratios = Lratios;
    newCluster.iso_dists = iso_distance;
    
    params.summary_plot = 2;
    sum_fig = create_summary_cluster_fig(newCluster,Spikes,spike_xy,params);
    
    resp = input('Keep new clustering (y/n, d for cluster dump)?','s');
    if strcmpi(resp,'Y')
        fprintf('Saving cluster details\n');
            rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
            load(rclust_dat_name);
        newCluster.base_block = RefClusters{probe_num}.base_block;
        Clusters{probe_num} = newCluster;
%         Clusters{probe_num}.base_block = RefCl
        save(cur_clust_data,'Clusters');
        
        if block_num == RefClusters{probe_num}.base_block
            fprintf('Saving to RefClusters\n');
            RefClusters{probe_num} = newCluster;
            save(rclust_dat_name,'RefClusters');
            
            fillPage(gcf,'papersize',[14 8]);
            pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
            print(pname,'-dpng');
        else
            resp = input('Save to RefClusters (y/n)?','s');
            if strcmpi(resp,'y')
                fprintf('Saving to RefClusters\n');
                rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
                load(rclust_dat_name);
                RefClusters{probe_num} = newCluster;
                RefClusters{probe_num}.base_block = block_num;
                save(rclust_dat_name,'RefClusters');
                
                fillPage(gcf,'papersize',[14 8]);
                pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
                print(pname,'-dpng');
            end
        end
        close(sum_fig);
        
    elseif strcmpi(resp,'D')
        fprintf('Dumping cluster output for inspection\n');
        Cdump = newCluster;
        if ishandle(sum_fig)
            close(sum_fig);
        end
    else
        fprintf('Keeping original clustering\n');
        if ishandle(sum_fig)
            close(sum_fig);
        end
    end
end
