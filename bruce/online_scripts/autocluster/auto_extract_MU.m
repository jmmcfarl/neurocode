clear all
close all
addpath('~/James_scripts/autocluster/');

global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData
Expt_name = 'M0050';

data_loc = '/media/NTlab_data3/Data/bruce/';
% data_loc = '/home/james/Data/bruce/';

base_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
init_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/init'];

if ~exist(base_save_dir,'dir');
    mkdir(base_save_dir);
end
if ~exist(init_save_dir,'dir');
    mkdir(init_save_dir);
end

%location of FullV files
data_dir = [data_loc Expt_name];

%location of Expts.mat files
% data_dir2 = ['~/Data/bruce/' Expt_name];
data_dir2 = ['/media/NTlab_data3/Data/bruce/' Expt_name];

Vloaded = nan;

cd(data_dir2);
% if Expt_name(1) == 'G';
%     load(sprintf('jbe%sExpts.mat',Expt_name));
%     n_probes = 96;
% elseif Expt_name(1) == 'M'
%     load(sprintf('lem%sExpts.mat',Expt_name));
    load(sprintf('jbe%sExpts.mat',Expt_name));
    n_probes = 24;
% end

%%
poss_base_blocks = [2]; %set of blocks to try fitting initial models on
target_probes = 1:n_probes;

%% SET CLUSTERING PARAMETERS
clear clust_params
clust_params.use_ms = 0; %mean-subtraction on templates
clust_params.target_rate = 50; %target rate for spike detection

%% PERFORM INITIAL CLUSTERING
for bb = 1:length(poss_base_blocks) %loop over initial set of blocks
    cur_base_block = poss_base_blocks(bb);
    full_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',cur_base_block)];
    %     second_dat_name = [base_save_dir sprintf('/Block%d_initClusters.mat',cur_base_block)];
    if exist(full_dat_name,'file');
        fprintf('Loading clusters for block %d\n',cur_base_block);
        load(full_dat_name,'Clusters');
    else
        fprintf('Initiating new cluster block\n');
        Clusters = cell(n_probes,1);
    end
    
    all_probe_fig = figure('visible','off');
    n_cols = ceil(sqrt(n_probes)); n_rows = ceil(n_probes/n_cols);
    for probe_num = target_probes
        fprintf('\nClustering probe %d of %d\n',probe_num,n_probes);
        
        if Expt_name(1) == 'G' %for UTAH array data Load in Voltage signasl for each probe
            loadedData = [data_dir sprintf('/Expt%d.p%dFullV.mat',cur_base_block,probe_num)];
            use_chs = [];
            cluster_details = detect_MU(loadedData,clust_params,use_chs);
        else %for Laminar probe data load in all voltage signals for a given block
            sfile_name = [data_dir sprintf('/Expt%dFullV.mat',cur_base_block)];
            use_chs = [probe_num];
            if Vloaded ~= cur_base_block
                fprintf('Loading data file %s\n',sfile_name);
                [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
                Vloaded = cur_base_block;
                cluster_details = detect_MU(loadedData,clust_params,probe_num);
            end
        end
        cluster_details.base_block = cur_base_block;
        Clusters{probe_num} = cluster_details;
        fprintf('Saving cluster details\n');
        save(full_dat_name,'Clusters');
        %         save(second_dat_name,'Clusters');
    end
end

%% save refclusters
RefClusters = Clusters;
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
save(rclust_dat_name,'RefClusters');

%%
target_blocks = [1:48];
force_new_clusters = false;

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name);

all_clust_means = cell(n_probes,1);
all_clust_stds = cell(n_probes,1);
for bb = target_blocks
    %for LP load all Voltage signals for this block
    if Expt_name(1) == 'M'
        sfile_name = [data_dir sprintf('/Expt%dFullV.mat',bb)];
        if Vloaded ~= bb
            fprintf('Loading data file %s\n',sfile_name);
            [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
            Vloaded = bb;
        end
    end
    
    %load existing clusters for this block
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',bb)];
    if exist(cur_dat_name,'file') && ~force_new_clusters
        fprintf('Loading clusters for block %d\n',bb);
        load(cur_dat_name);
    else
        fprintf('Initializing clusters for block %d\n',bb);
        Clusters = cell(n_probes,1);
    end
    
    for probe_num = target_probes
        fprintf('Applying clustering for probe %d\n',probe_num);
        if Expt_name(1) == 'G'
            loadedData = [data_dir sprintf('/Expt%d.p%dFullV.mat',bb,probe_num)];
        end
        
        cur_base_block = RefClusters{probe_num}.base_block;
        cur_trig_thresh = RefClusters{probe_num}.trig_thresh;
        if Expt_name(1) == 'G'
        new_cluster = detect_MU(loadedData,clust_params,[],cur_trig_thresh);   
        else
        new_cluster = detect_MU(loadedData,clust_params,probe_num,cur_trig_thresh);   
        end
        Clusters{probe_num} = new_cluster;
        
    end
    
    fprintf('Saving clusters for block %d\n',bb);
    save(cur_dat_name,'Clusters');
end
