clear all
close all
cd ~/James_scripts/autocluster/
Expt_name = 'G095';

%%
clust_params.gmm_inits = 20;
clust_params.reg_lambda = 0;
clust_params.min_Pcomp = 0.005;

probe_num = 44;
base_block = 14;
dat_name = sprintf('Expt%d.p%d_Clusters.mat',base_block,probe_num);
sfile_name = sprintf('Expt%d.p%dFullV.mat',base_block,probe_num);

[cluster_details,spike_features,sum_fig] = detect_and_cluster_2comp(sfile_name,clust_params);
cluster_details.base_block = base_block;
Clusters{probe_num} = cluster_details;

fprintf('Saving cluster details\n');
save(dat_name,'Clusters');

spike_xy = spike_features*cluster_details.xy_projmat;
N_spks = size(spike_xy,1);
su_inds = find(ismember(cluster_details.comp_idx,find(cluster_details.cluster_labels == 2)));
mu_inds = setdiff(1:N_spks,su_inds);

%save unit cluster plot
fillPage(sum_fig,'papersize',[14 8]);
putative_good(probe_num) = 1;
pname = sprintf('Probe%d_Block%d_initclust',probe_num,base_block);
print(sum_fig,pname,'-dpng');
close(sum_fig);

%% DO ANY RETRIGGERING
retrig_ch = 44;
retrig_rate = 50; %retriggering at higher rates than 50 makes this one mess up it seems...
retrig_sign = -1;
reapply = 1; %tries to reapply existing model params
block_num = base_block;

load(dat_name);
fprintf('Retriggering probe %d\n',retrig_ch);
cluster_details = Clusters{retrig_ch};
cur_clust_params = cluster_details.params;
cur_clust_params.target_rate = retrig_rate;
cur_clust_params.thresh_sign = retrig_sign;
% cur_clust_params.min_Pcomp = 0.0001;
cur_clust_params.summary_plot = 2; %make summary plot visible

if reapply == 0
    [cluster_details,spike_features,sum_fig] = detect_and_cluster_2comp(cluster_details.rawV_file,cur_clust_params);
else
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(cluster_details.rawV_file,cluster_details,cur_clust_params);
    sum_fig = create_summary_cluster_fig(cluster_details,Spikes,spike_xy,cur_clust_params);
    clear Spikes
end

resp = input('Keep new clustering (y/n)?','s');
if strcmpi(resp,'Y')
    fprintf('Saving cluster details\n');
    Clusters{retrig_ch} = cluster_details;
    Clusters{retrig_ch}.params = cur_clust_params;
    save(dat_name,'Clusters');

    fillPage(gcf,'papersize',[14 8]);
    pname = sprintf('Probe%d_Block%d_initclust',retrig_ch,block_num);
    print(pname,'-dpng');
    close(sum_fig);
else
    fprintf('Keeping original clustering\n');
    if ishandle(sum_fig)
    close(sum_fig);
    end
end

%% ADD NEW CLUSTER
close all

probe_num = 44;
block_num = base_block;

load(dat_name);

[new_cluster,used_new] = add_GMM_cluster(Clusters{probe_num}.rawV_file,Clusters{probe_num});
if used_new == 1
    Clusters{probe_num} = new_cluster;
    fprintf('Saving cluster details\n');
    save(dat_name,'Clusters');
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = sprintf('Probe%d_Block%d_initclust',probe_num,block_num);
        print(pname,'-dpng');
        close(sum_fig);
    end
end

%% SPLIT COMPONENTS
close all

probe_num = 44;
block_num = base_block;

load(dat_name);

[new_cluster,used_new] = split_GMM_cluster(Clusters{probe_num}.rawV_file,Clusters{probe_num});
if used_new == 1
    Clusters{probe_num} = new_cluster;
    fprintf('Saving cluster details\n');
    save(dat_name,'Clusters');
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = sprintf('Probe%d_Block%d_initclust',probe_num,block_num);
        print(pname,'-dpng');
        close(sum_fig);
    end
end

%% MERGE COMPONENTS
close all

probe_num = 44;
block_num = base_block;

load(dat_name);

[new_cluster,used_new] = merge_GMM_cluster(Clusters{probe_num}.rawV_file,Clusters{probe_num});
if used_new == 1
    Clusters{probe_num} = new_cluster;
    fprintf('Saving cluster details\n');
    save(cur_dat_name,'Clusters');
    [cluster_details,spike_features,spike_xy,Spikes] = apply_clustering(sfile_name,Clusters{probe_num},[],1);
    if block_num == cluster_details.base_block
        sum_fig = create_summary_cluster_fig(new_cluster,Spikes,spike_xy,Clusters{probe_num}.params);
        pname = sprintf('Probe%d_Block%d_initclust',probe_num,block_num);
        print(pname,'-dpng');
        close(sum_fig);
    end
end


