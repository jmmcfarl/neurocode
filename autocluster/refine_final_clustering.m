all_target_blocks = 1:40;
fprintf('Loading Ref_Clusters\n');
rclust_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_name);

%% JUST PLOT CLUSTER
block_num = 4;
probe_num = 15;
spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,block_num)];
plot_GMM_cluster(block_num,probe_num,spk_data_name);

%% SPLIT COMPONENTS
block_num = 28;
probe_num = 18;
spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,block_num)];
split_GMM_component(block_num,probe_num,spk_data_name);

%% DELETE COMPONENT
block_num = 32;
probe_num = 20;
spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,block_num)];
delete_GMM_component(block_num,probe_num,spk_data_name);

%% EXCLUSION CLUSTERING
block_num = 35;
probe_num = 10;
exclude_SUs = [2];

clear clust_params
clust_params.outlier_thresh = 5;

autocluster_excludeSU(block_num,probe_num,exclude_SUs,clust_params);

%% CYCLE THROUGH PROJECTIONS
block_num = 5;
probe_num = 19;
cycle_projections(block_num,probe_num);

%% SPLIT COMPONENTS
block_num = 14;
probe_num = 17;
use_2d = false;
use_proj = [2 3];
split_GMM_component_dimspec(block_num,probe_num,use_proj,use_2d);

%% CYCLE THROUGH CLUSTER PLOTS and do split/delete operations
close all

probe_num = 15;
block_num = RefClusters{probe_num}.base_block;
spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,block_num)];
f1 = plot_GMM_cluster(block_num,probe_num,spk_data_name);

target_blocks = setdiff(5:40,block_num);
for bb = 1:length(target_blocks)
    block_num = target_blocks(bb);
    spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,block_num)];
    f2 = plot_GMM_cluster(block_num,probe_num,spk_data_name);
    
    resp = input('Delete (d), Split/Relable (s)?','s');
    if strcmpi(resp,'d')
        delete_GMM_component(block_num,probe_num,spk_data_name);
    elseif strcmpi(resp,'s')
        split_GMM_component(block_num,probe_num,spk_data_name);
    end
    close(f2);
end

regenerate_allblock_xydensity(probe_num,all_target_blocks);
regenerate_allblock_xyscatters(probe_num,all_target_blocks);

%% apply new refcluster to target blocks using saved spike data
probe_number = 18;
base_block = RefClusters{probe_num}.base_block;
target_blocks = setdiff(all_target_blocks,base_block);

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name);

for bb = 1:length(target_blocks)
    spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,target_blocks(bb))];
    Spikes = load_spike_data(spk_data_name);
    
        
    %load existing clusters for this block
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',target_blocks(bb))];
    if exist(cur_dat_name,'file') && ~force_new_clusters
        fprintf('Loading clusters for block %d\n',target_blocks(bb));
        load(cur_dat_name);
    else
        fprintf('Initializing clusters for block %d\n',target_blocks(bb));
        Clusters = cell(n_probes,1);
    end

    try
    [new_cluster,spike_features,spike_xy] = ...
        apply_clustering(loadedData,RefClusters{probe_num},[],0,Spikes);
        new_cluster.failed = false;
    catch
        new_cluster.failed = true;
        fprintf('Couldnt cluster probe %d block %d!\n',probe_num,target_blocks(bb));
    end
    Clusters{probe_num} = new_cluster;
    
    fprintf('Saving clusters for block %d\n',target_blocks(bb));
    save(cur_dat_name,'Clusters');

end

regenerate_allblock_xydensity(probe_num,all_target_blocks);
regenerate_allblock_xyscatters(probe_num,all_target_blocks);

%% check cluster alignment
probe_num = 18;
target_blocks = 1:40;
check_cluster_alignment(probe_num,target_blocks);
regenerate_allblock_xyscatters(probe_num,all_target_blocks);
