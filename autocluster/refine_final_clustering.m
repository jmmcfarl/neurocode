
%% JUST PLOT CLUSTER
block_num = 11;
probe_num = 23;
spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,block_num)];
plot_GMM_cluster(block_num,probe_num,spk_data_name);

%% SPLIT COMPONENTS
block_num = 49;
probe_num = 23;
spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,block_num)];
split_GMM_component(block_num,probe_num,spk_data_name);

%% DELETE COMPONENT
block_num = 32;
probe_num = 20;
spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,block_num)];
delete_GMM_component(block_num,probe_num,spk_data_name);

%% EXCLUSION CLUSTERING
block_num = 14;
probe_num = 17;
exclude_SUs = [2 3];

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
