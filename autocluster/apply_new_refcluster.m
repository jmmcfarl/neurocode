function apply_new_refcluster(target_blocks,target_probe)

global base_save_dir spkdata_dir Expt_name

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name);

for bb = target_blocks
    %load existing clusters for this block
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',bb)];
    if exist(cur_dat_name,'file')
        fprintf('Loading clusters for block %d\n',bb);
        load(cur_dat_name);
    else
        error('No existing cluster file');
    end
    
    fprintf('Applying clustering for block %d, probe %d\n',bb,target_probe);
    spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',target_probe,bb)];
    Spikes = load_spike_data(spk_data_name);
    try
        [new_cluster,spike_features,spike_xy,Spikes] = apply_clustering([],RefClusters{target_probe},[],0,Spikes);
    catch
        new_cluster.failed = 1;
        fprintf('Couldnt cluster probe %d block %d!\n',target_probe,bb);
    end
    
    Clusters{target_probe} = new_cluster;
    fprintf('Saving clusters for block %d\n',bb);
    save(cur_dat_name,'Clusters');

end

%%
% regenerate_allblock_xydensity(target_probe,target_blocks);
% regenerate_allblock_xyscatters(target_probe,target_blocks);

