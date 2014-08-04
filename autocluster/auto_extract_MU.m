clear all
close all
addpath('~/James_scripts/autocluster/');

global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData
Expt_name = 'M266';

% data_loc = '/media/NTlab_data1/Data/bruce/';
data_loc = '/home/james/Data/bruce/';

base_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering2'];
init_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering2/init'];

if ~exist(base_save_dir,'dir');
    mkdir(base_save_dir);
end
if ~exist(init_save_dir,'dir');
    mkdir(init_save_dir);
end

%location of FullV files
data_dir = [data_loc Expt_name];

%location of Expts.mat files
data_dir2 = ['~/Data/bruce/' Expt_name];

Vloaded = nan;

cd(data_dir2);
if Expt_name(1) == 'G';
    load(sprintf('jbe%sExpts.mat',Expt_name));
    n_probes = 96;
elseif Expt_name(1) == 'M'
    load(sprintf('lem%sExpts.mat',Expt_name));
    n_probes = 24;    
end

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
        else %for Laminar probe data load in all voltage signals for a given block
            sfile_name = [data_dir sprintf('/Expt%dFullV.mat',cur_base_block)];
            use_chs = [probe_num];
            if Vloaded ~= cur_base_block
                fprintf('Loading data file %s\n',sfile_name);
                [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
                Vloaded = cur_base_block;
            end
        end
        cluster_details = detect_MU(loadedData,clust_params,probe_num);
        cluster_details.base_block = cur_base_block;
        Clusters{probe_num} = cluster_details;
        fprintf('Saving cluster details\n');
        save(full_dat_name,'Clusters');
%         save(second_dat_name,'Clusters');
    end    
end

