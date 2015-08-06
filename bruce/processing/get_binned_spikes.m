function [all_binned_mua,all_binned_sua,Clust_data,all_su_spk_times,all_su_spk_inds,all_mu_spk_times] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params,rec_block_range)
% [all_binned_mua,all_binned_sua,Clust_data,all_su_spk_times,all_su_spk_inds,all_mu_spk_times] = ...
%     get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
%     all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params,<rec_block_range>)
% grab clustered spike data and bin it
% INPUTS:
%     cluster_dir: name of directory containing final_cluster.mat file
%     all_spk_times: Nx1 cell array (N is number of channels), containing all the spike times from each channel
%     all_clust_ids: Nx1 cell array containing the cluster classification for each spike
%     all_spk_inds: index value, within the original fullV file, for each spike
%     all_t_axis: time bin centers
%     all_t_bin_edges: bin edges of time axis
%     all_bin_edge_pts: indices of bin edges occuring on trial boundaries
%     cur_block_set: set of block numbers being used
%     all_blockvec: time series of block numbers
%     clust_params: struct of clustering parameters, at least needs to have an 'exclude_adjacent' field
%     <rec_block_range>: if using a subset of blocks for this experiment, needed to handle that bookkeeping
% OUTPUTS:
%     all_binned_mua: TxN binned spike counts for MUA on each channel
%     all_binned_sua: binned spike counts for each single unit
%     Clust_data: struct containing cluster info for each SU
%     all_su_spk_times: cell array with spike times for each SU
%     all_su_spk_inds: cell array containing spk indices (relative to FullVs) for each SU spk
%     all_mu_spk_times: cell array with spike times for each MU

if nargin < 11
    rec_block_range = nan;
end

exclude_adjacent = clust_params.exclude_adjacent;
n_probes = clust_params.n_probes;

% LOAD FINAL CLUSTER FILE
fname = [cluster_dir '/final_cluster.mat'];
if exist(fname,'file')
    load(fname);
    SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
    fprintf('%d SUs Clustered\n',length(SU_numbers));
else
    disp('No Cluster data found.');
end

%if the recording spans a specific range of blocks, select these
if ~isnan(rec_block_range)
   SU_ID_mat = SU_ID_mat(rec_block_range,:);
   SU_allBlock_Data = SU_allBlock_Data(:,rec_block_range);
end

%for SU probes
all_binned_sua = nan(length(all_t_axis),length(SU_numbers));
cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
SU_ID_mat = SU_ID_mat(cur_used_blocks,:);
[CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));

%cycle over SUs
su_probes = nan(1,length(SU_numbers));
all_su_spk_times = cell(length(SU_numbers),1);
all_su_spk_inds = cell(length(SU_numbers),1);
for ss = 1:length(SU_numbers)
    used_clust_set = unique(CC(SU_ID_mat==SU_numbers(ss))); %set of clusters used to capture this SU
    SU_block_probes(ss,:) = nan(1,length(cur_block_set));
    cur_su_spk_times = [];
    cur_su_spk_inds = [];
    cur_blocks = [];
    for cc = 1:length(used_clust_set) %cycle over clusters used for this SU
        cur_clust = used_clust_set(cc);
        cur_probe = SU_clust_data(cur_clust).probe_num;
        cur_clust_label = SU_clust_data(cur_clust).cluster_label;
        cur_blocks = [cur_blocks; find(SU_ID_mat(:,cur_clust) == SU_numbers(ss))]; %blocks with this cluster isolated
        SU_block_probes(ss,cur_blocks) = cur_probe; %store which probe was used for this cluster on these blocks
                
        all_su_inds = all_clust_ids{cur_probe} == cur_clust_label; %all spikes in this cluster
        cur_su_spk_times = all_spk_times{cur_probe}(all_su_inds); 
        cur_su_spk_inds = all_spk_inds{cur_probe}(all_su_inds);
        spk_block_inds = round(interp1(all_t_axis,all_blockvec,cur_su_spk_times));
        cur_su_spk_times = cur_su_spk_times(ismember(spk_block_inds,cur_blocks));
        cur_su_spk_inds = cur_su_spk_inds(ismember(spk_block_inds,cur_blocks));
        
        all_su_spk_times{ss} = cat(1,all_su_spk_times{ss},cur_su_spk_times(:));
        all_su_spk_inds{ss} = cat(1,all_su_spk_inds{ss},cur_su_spk_inds(:));
    end
    
    %now get binned SU spike times
    if ~isempty(all_su_spk_times{ss})
        cur_suahist = histc(all_su_spk_times{ss},all_t_bin_edges);
        cur_suahist(all_bin_edge_pts) = []; %account for trial edges in binning
        cur_id_set = ismember(all_blockvec,cur_blocks);
        all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
        su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:)))); %use the mode as a best rep of the probe associated with this SU
    end
end

double_spike_buffer = 3; %number of samples (in either direction) to exclude double spikes from adjacent-probe SUs
all_binned_mua = nan(length(all_t_axis),n_probes);
all_mu_spk_times = cell(n_probes,1);
for cc = 1:n_probes
    %this is the set of blocks where this probe had an SU, and the
    %corresponding SU numbers
    cur_set = find(SU_block_probes == cc);
    if ~isempty(cur_set)
        [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
    else
        cur_SS = [];
    end
    unique_su_nums = unique(cur_SS); %su numbers clustered on this probe
    cur_mua_inds = find(all_clust_ids{cc} >= 1); %anything that's not noise counts as MU (clust_id of 1 or above)
    
    %remove spikes from isolated SUs on the same probe from the MU
    for ss = 1:length(unique_su_nums)
        cur_mua_inds(ismember(all_spk_inds{cc}(cur_mua_inds),all_su_spk_inds{unique_su_nums(ss)})) = [];
    end
    
    %if also excluding spikes from SUs on adjacent probes
    if exclude_adjacent
        nearby_probes = [cc-1 cc+1]; nearby_probes(nearby_probes < 1 | nearby_probes > n_probes) = []; %set of adjacent probes
        cur_set = find(ismember(SU_block_probes,nearby_probes)); 
        if ~isempty(cur_set)
            [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
        else
            cur_SS = [];
        end
        unique_su_nums = unique(cur_SS); %set of SU numbers picked up on adjacent probes
        if ~isempty(unique_su_nums)
            double_spikes = [];
            for ss = 1:length(unique_su_nums)
                cur_blocked_inds = bsxfun(@plus,all_su_spk_inds{unique_su_nums(ss)},-double_spike_buffer:double_spike_buffer); %these indices define a buffer around adjacent SU spikes where we shouldnt count MU
                double_spikes = [double_spikes; find(ismember(all_spk_inds{cc}(cur_mua_inds),cur_blocked_inds))];
            end
            double_spikes = unique(double_spikes);
            fprintf('Eliminating %d of %d double spikes in MUA\n',length(double_spikes),length(cur_mua_inds));
            cur_mua_inds(double_spikes) = [];
        end
    end
    
    %bin our MU
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    all_mu_spk_times{cc} = cat(1,all_mu_spk_times{cc},all_spk_times{cc}(cur_mua_inds));
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end

%% store cluster info

Clust_data.SU_block_probes = SU_block_probes;
Clust_data.SU_numbers = SU_numbers;
Clust_data.SU_ID_mat = SU_ID_mat;
Clust_data.cluster_stats = SU_clust_data;
Clust_data.SU_probes = su_probes;

if exist('SU_allBlock_Data','var')
    SU_allBlock_Data = SU_allBlock_Data(:,cur_used_blocks);
    
    Clust_data.SU_Lratios = reshape([SU_allBlock_Data.Lratios],length(SU_numbers),[]);
    Clust_data.SU_isodists = reshape([SU_allBlock_Data.isoDists],length(SU_numbers),[]);
    Clust_data.SU_isoRel = reshape([SU_allBlock_Data.isoReliable],length(SU_numbers),[]);
    Clust_data.SU_dprimes = reshape([SU_allBlock_Data.dprimes],length(SU_numbers),[]);
    Clust_data.SU_refract = reshape([SU_allBlock_Data.refract],length(SU_numbers),[]);
    
else
    fprintf('NO SU_allBlock_Data, Ignoring!!\n');
    Clust_data.SU_Lratios = nan;
    Clust_data.SU_isodists = nan;
    Clust_data.SU_isoRel = nan;
    Clust_data.SU_dprimes = nan;
    Clust_data.SU_refract = nan;
end

