fprintf('Loading Ref_Clusters\n');
rclust_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_name);

%% RETRIGGERING
block_num = 55;
probe_num = 86;
trig_rate = 100;
trig_sign = -1;
reapply = 0;
retrigger_and_cluster(RefClusters,block_num,probe_num,trig_rate,trig_sign,reapply);

%% SPLIT COMPONENTS
block_num = 68;
probe_num = 17;
split_GMM_component(block_num,probe_num);

%% DELETE COMPONENT
block_num = 45;
probe_num = 23;
delete_GMM_component(block_num,probe_num);
