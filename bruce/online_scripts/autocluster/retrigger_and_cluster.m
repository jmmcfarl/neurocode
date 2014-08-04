function retrigger_and_cluster(RefClusters,block_num,probe_num,trig_rate,trig_sign,reapply)

global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData

if nargin < 5 || isempty(trig_sign)
    trig_sign = -1;
end
if nargin < 6 || isempty(reapply)
    reapply = 0;
end

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

fprintf('Retriggering probe %d\n',probe_num);
cur_clust_params = RefClusters{probe_num}.params;
cur_clust_params.summary_plot = 2; %make summary plot visible
cur_clust_params.target_rate = trig_rate;
cur_clust_params.thresh_sign = trig_sign;

cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
load(cur_dat_name,'Clusters');

if reapply == 0
    if n_probes == 24
        use_chs = [probe_num-1 probe_num probe_num + 1];
        use_chs(use_chs < 1 | use_chs > n_probes) = [];
    else
        use_chs = [];
    end
    [cluster_details,spike_features,sum_fig] = detect_and_cluster_init(loadedData,cur_clust_params,use_chs);
else
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,RefClusters{probe_num},cur_clust_params,-1);
    sum_fig = create_summary_cluster_fig(cluster_details,Spikes,spike_xy,cur_clust_params);
    clear Spikes
end

cluster_details.base_block = RefClusters{probe_num}.base_block;

resp = input('Keep new clustering (y/n)?','s');
if strcmpi(resp,'Y')
    fprintf('Saving cluster details\n');
    Clusters{probe_num} = cluster_details;
    save(cur_dat_name,'Clusters');
    
    if block_num == RefClusters{probe_num}.base_block
        fprintf('Saving to RefClusters\n');
        rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
        RefClusters{probe_num} = cluster_details;
        save(rclust_dat_name,'RefClusters');
        
        fillPage(gcf,'papersize',[14 8]);
        pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
        print(pname,'-dpng');
    else
        resp = input('Save to RefClusters (y/n)?','s');
        if strcmpi(resp,'y')
            fprintf('Saving to RefClusters\n');
            rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
            RefClusters{probe_num} = cluster_details;
            RefClusters{probe_num}.base_block = block_num;
            save(rclust_dat_name,'RefClusters');
            
            fillPage(gcf,'papersize',[14 8]);
            pname = [init_save_dir sprintf('/Probe%d_Block%d_initclust',probe_num,block_num)];
            print(pname,'-dpng');
        end
    end
    close(sum_fig);
    
else
    fprintf('Keeping original clustering\n');
    if ishandle(sum_fig)
        close(sum_fig);
    end
end
