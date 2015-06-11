clear all
close all
addpath('~/James_scripts/autocluster/');

global data_dir base_save_dir init_save_dir spkdata_dir Expt_name monk_name rec_type Vloaded n_probes loadedData raw_block_nums
Expt_name = 'M320';
monk_name = 'lem';
rec_type = 'LP';

rec_number = 1;

% block_set = [1:100];

Expt_num = str2num(Expt_name(2:end));

data_loc = '/media/NTlab_data3/Data/bruce/';

%location of Expts.mat files
data_dir2 = [data_loc Expt_name];

spkdata_dir = [data_loc Expt_name '/spikes/'];

base_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
if rec_number > 1 %if you're splitting the recording into multiple separate chunks for clustering
   base_save_dir = [base_save_dir sprintf('/rec%d',rec_number)]; 
end
init_save_dir = [base_save_dir '/init'];

if ~exist(base_save_dir,'dir');
    mkdir(base_save_dir);
end
if ~exist(init_save_dir,'dir');
    mkdir(init_save_dir);
end
if ~exist(spkdata_dir,'dir')
    mkdir(spkdata_dir);
end

%location of FullV files
data_dir = [data_loc Expt_name];

Vloaded = nan;
cd(data_dir2);
if strcmp(Expt_name,'G029')
    load('G029Expts.mat');
else
    load(sprintf('%s%sExpts.mat',monk_name,Expt_name));
end


if strcmp(rec_type,'UA')
    n_probes = 96;
elseif strcmp(rec_type,'LP')
    n_probes = 24;
end

target_probes = 1:n_probes;
force_new_clusters = false; %if you want to ov

elen = cellfun(@(x) length(x),Expts);
target_blocks = find(elen > 0);
if isfield(Expts{1}.Header,'exptno')
    raw_block_nums = cellfun(@(X) X.Header.exptno,Expts,'uniformoutput',1); %block numbering for EM/LFP data sometimes isnt aligned with Expts struct
else
    raw_block_nums = 1:max(target_blocks);
end

%don't apply to blocks where we dont have the FullV data
missing_Vdata = [];
for bb = target_blocks
    if strcmp(rec_type,'LP')
        check_name = [data_dir sprintf('/Expt%dFullV.mat',raw_block_nums(bb))];
    elseif strcmp(rec_type,'UA')
        check_name = [data_dir sprintf('/Expt%d.p1FullV.mat',raw_block_nums(bb))];
    end
    if ~exist(check_name,'file')
        missing_Vdata = [missing_Vdata bb];
    end
end
target_blocks(ismember(target_blocks,missing_Vdata)) = [];


n_blocks = max(target_blocks);
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);

global full_save_dir

full_save_dir = [base_save_dir '/full'];
if ~exist(full_save_dir,'dir');
    mkdir(full_save_dir);
end

if strcmp(Expt_name,'M281')
    target_blocks(target_blocks == 13) = [];
end
if strcmp(Expt_name,'M289')
    target_blocks(target_blocks == 14 | target_blocks == 32) = [];
end

if exist('block_set','var')
target_blocks(~ismember(target_blocks,block_set)) = [];
end
%%
% ADD CLUSTER PERMUTATION CHECKING INTO THE APPLY REF CLUSTERS FUNCTION!

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name);

if ~exist('full_scatter_fig','var')
    full_scatter_fig = nan;
end
% if ~exist('full_dens_fig','var')
%     full_dens_fig = nan;
% end

all_clust_means = cell(n_probes,1);
all_clust_stds = cell(n_probes,1);
for bb = target_blocks
    %for LP load all Voltage signals for this block
    if strcmp(rec_type,'LP')
        sfile_name = [data_dir sprintf('/Expt%dFullV.mat',raw_block_nums(bb))];
        if Vloaded ~= raw_block_nums(bb)
            fprintf('Loading data file %s\n',sfile_name);
            [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
            Vloaded = raw_block_nums(bb);
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
        if strcmp(rec_type,'UA')
            loadedData = [data_dir sprintf('/Expt%d.p%dFullV.mat',raw_block_nums(bb),probe_num)];
        end
        
        %look for existing full scatter figure and open if it exists
        pfname_sc = [full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
        if exist(pfname_sc,'file') && ~ishandle(full_scatter_fig)
            full_scatter_fig = open(pfname_sc); set(full_scatter_fig,'visible','off');
        else
            full_scatter_fig = figure('visible','off');
        end
        %         %look for existing full density figure and open if it exists
        %         pfname_de = [full_save_dir sprintf('/Probe%d_fullclust_dens.fig',probe_num)];
        %         if exist(pfname_de,'file') && ~ishandle(full_dens_fig)
        %             full_dens_fig = open(pfname_de); set(full_dens_fig,'visible','off');
        %         else
        %             full_dens_fig = figure('visible','off');
        %         end
        
        cur_base_block = RefClusters{probe_num}.base_block;
        try
            %     if bb == cur_base_block
            %         [new_cluster,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,RefClusters{probe_num},[],1);
            %     else
            [new_cluster,spike_features,spike_xy,Spikes] = apply_clustering(loadedData,RefClusters{probe_num});
            %     end
          
            spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',probe_num,bb)];
            store_spike_data(Spikes,spk_data_name,1);
            
        catch
            new_cluster.failed = 1;
            fprintf('Couldnt cluster probe %d block %d!\n',probe_num,bb);
        end
        
        
        if ~new_cluster.failed
            N_spks = size(spike_xy,1);
            spike_labels = new_cluster.spike_clusts;
            mu_inds = find(spike_labels == 1);
            out_inds = find(spike_labels == -1);
            N_sus = length(unique(spike_labels(spike_labels > 1)));
            cmap = cluster_cmap(N_sus);
            
            set(0,'CurrentFigure',full_scatter_fig);
            subplot(n_cols,n_rows,bb);hold off
            plot(spike_xy(mu_inds,1),spike_xy(mu_inds,2),'k.'); hold on
            plot(spike_xy(out_inds,1),spike_xy(out_inds,2),'r.');
            for ii = 1:N_sus
                plot(spike_xy(spike_labels == ii + 1,1),spike_xy(spike_labels == ii + 1,2),'.','color',cmap(ii,:));
            end
            set(gca,'xtick',[],'ytick',[]);axis tight
            for ii = 1:length(new_cluster.cluster_labels)
                if ~isnan(new_cluster.gmm_xyMeans(ii,1))
                h1 = plot_gaussian_2d(new_cluster.gmm_xyMeans(ii,:)',squeeze(new_cluster.gmm_xySigma(:,:,ii)),[2],'r',1);
                end
            end
            if bb == cur_base_block
                title(['Block #',int2str(bb)],'Color','r');
            else
                title(['Block #',int2str(bb)],'Color','k');
            end
            saveas(full_scatter_fig,pfname_sc);
            close(full_scatter_fig);
            
            all_clust_means{probe_num}(bb,1) = mean(spike_xy(mu_inds,1));
            all_clust_stds{probe_num}(bb,1) = std(spike_xy(mu_inds,1));
            for ii = 1:N_sus
                all_clust_means{probe_num}(bb,ii+1) = mean(spike_xy(spike_labels==ii+1,1));
                all_clust_stds{probe_num}(bb,ii+1) = std(spike_xy(spike_labels==ii+1,1));
            end
        else
            close(full_scatter_fig);
        end
        
        %         set(0,'CurrentFigure',full_dens_fig);
        %         subplot(n_cols,n_rows,bb);hold on
        %         [handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
        %         if bb == cur_base_block
        %             dens_xrange(probe_num,:) = minmax(details.x);
        %             dens_yrange(probe_num,:) = minmax(details.y);
        %         end
        %         set(gca,'xtick',[],'ytick',[]);
        %         if bb == cur_base_block
        %             title(['Block #',int2str(bb)],'Color','r');
        %         else
        %             title(['Block #',int2str(bb)],'Color','k');
        %         end
        %         saveas(full_dens_fig,pfname_de);
        %         close(full_dens_fig);
        
        probe_iso_quality(bb,probe_num,:) = [new_cluster.dprime new_cluster.LL new_cluster.Lratios(1) new_cluster.iso_dists(1)];
        
        Clusters{probe_num} = new_cluster;
        
    end
    
    fprintf('Saving clusters for block %d\n',bb);
    save(cur_dat_name,'Clusters');
end

%% RESCALE SCATTER PLOTS TO HAVE CONSISTENT AXES
% for probe_num = target_probes
%     fprintf('Rescaling plots for probe %d\n',probe_num);
%
%      %look for existing full scatter figure and open if it exists
%     pfname_sc = [full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
%     if exist(pfname_sc,'file') %&& ~ishandle(full_scatter_fig)
%         full_scatter_fig = open(pfname_sc);
%     end
%
%     %look for existing full scatter figure and open if it exists
%     pfname_de = [full_save_dir sprintf('/Probe%d_fullclust_dens.fig',probe_num)];
%     if exist(pfname_de,'file') %&& ~ishandle(full_dens_fig)
%         full_dens_fig = open(pfname_de); set(full_dens_fig,'visible','on');
%     end
%
%     cur_base_block = RefClusters{probe_num}.base_block;
%
%     set(0,'CurrentFigure',full_scatter_fig);
%     subplot(n_cols,n_rows,cur_base_block);
%     xl = xlim(); yl = ylim();
%     for bb = target_blocks
%         subplot(n_cols,n_rows,bb);
%         xlim(xl); ylim(yl);
%     end
%     saveas(full_scatter_fig,pfname_sc);
%     close(full_scatter_fig);
%
% end

%% PRINT OUT PLOTS ACROSS ALL BLOCKS FOR EACH PROBE
% close all
% for probe_num = target_probes
%     fprintf('Probe %d\n',probe_num);
%
%     %look for existing full scatter figure and open if it exists
%     pfname_sc = [full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
%     if exist(pfname_sc,'file') %&& ~ishandle(full_scatter_fig)
%         full_scatter_fig = open(pfname_sc);
%     end
%     fillPage(full_scatter_fig,'papersize',[14 14]);
%     fname = [full_save_dir sprintf('/Probe%d_fullclust_scatter',probe_num)];
%     print(full_scatter_fig,fname,'-dpng');
%     close(full_scatter_fig);
%
% end

%% PRINT FULL DENSITY FIGURES
for probe_num = target_probes
    fprintf('Probe %d of %d\n',probe_num,length(target_probes));
    regenerate_allblock_xydensity(probe_num,target_blocks);
    regenerate_allblock_xyscatters(probe_num,target_blocks);
end

%% PRINT ALL-PROBE PLOTS FOR EACH BLOCK
for bb = target_blocks
    fprintf('Probe %d of %d\n',probe_num,length(target_probes));
    regenerate_allprobe_xydensity(target_probes,bb);
    regenerate_allprobe_xyscatters(target_probes,bb);
end

%% CHECK CLUSTER ASSIGNMENTS ACROSS ALL BLOCKS FOR ALL PROBES HAVING AT LEAST 2 SUs
for probe_num = target_probes
    N_sus = max(RefClusters{probe_num}.cluster_labels) - 1;
    if N_sus > 1
        check_cluster_alignment(probe_num,target_blocks);
        regenerate_allblock_xyscatters(probe_num,target_blocks);
    end
end
