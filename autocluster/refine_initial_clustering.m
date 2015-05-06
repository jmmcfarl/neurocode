
%% RETRIGGERING
block_num = 40;
probe_num = 24;
trig_rate = 300;
trig_sign = -1;
reapply = 0;
clear add_params
add_params.try_features = [1 2 4];
% add_params.outlier_thresh = 5;

fprintf('Loading Ref_Clusters\n');
rclust_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_name);

Cdump = retrigger_and_cluster(RefClusters,block_num,probe_num,trig_rate,trig_sign,reapply,add_params);

%% SPLIT COMPONENTS
block_num = 15;
probe_num = 8;
split_GMM_component(block_num,probe_num);

%% DELETE COMPONENT
block_num = 15;
probe_num = 5;
delete_GMM_component(block_num,probe_num);

%% EXCLUSION CLUSTERING
block_num = 38;
probe_num = 17;
exclude_SUs = [2];

clear clust_params
clust_params.outlier_thresh = 5;
clust_params.trig_rate = 150;
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
