function [] = check_cluster_alignment(probe_num,target_blocks)

global base_save_dir     

fprintf('Aligning cluster assignments across blocks for probe %d\n',probe_num);
n_blocks = length(target_blocks);

fprintf('Loading RefClusters\n');
rclust_dat_name = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_dat_name);

N_sus = max(RefClusters{probe_num}.cluster_labels) - 1;

cmap = cluster_cmap(N_sus);
N_samps = length(RefClusters{probe_num}.params.spk_pts);
N_chs = length(RefClusters{probe_num}.use_chs);
ref_mean_spike = RefClusters{probe_num}.mean_spike(:,2:end);
ms_size = size(ref_mean_spike);
all_mean_spike = nan(n_blocks,ms_size(1),ms_size(2));
for bb = 1:length(target_blocks)
    cur_data = [base_save_dir sprintf('/Block%d_Clusters.mat',target_blocks(bb))];
    load(cur_data,'Clusters');
    if ~Clusters{probe_num}.failed
        all_mean_spike(bb,:,:) = Clusters{probe_num}.mean_spike(:,2:end);
    end
end

max_n_sus = 8;
if N_sus > max_n_sus
    error('Too many possible SU permutations to consider!');
end
perm_set = perms(1:N_sus);
n_perms = size(perm_set,1);
block_Ctrace = nan(1,n_perms);
best_perms = nan(1,n_blocks);
for bb = 1:length(target_blocks)
    uset = find(~isnan(all_mean_spike(bb,1,:)));
    Cmat = corr(ref_mean_spike,squeeze(all_mean_spike(bb,:,:)));
    for ii = 1:n_perms
        block_Ctrace(ii) = nansum(diag(Cmat(perm_set(ii,:),:)));
    end
    [~,best_perms(bb)] = max(block_Ctrace);
end

for bb = 1:length(target_blocks)
    fprintf('Relabeling clusters for block %d\n',target_blocks(bb));
    cur_data = [base_save_dir sprintf('/Block%d_Clusters.mat',target_blocks(bb))];
    load(cur_data,'Clusters');
    if ~Clusters{probe_num}.failed
        prev_labels = Clusters{probe_num}.cluster_labels;
        prev_clusts = Clusters{probe_num}.spike_clusts;
        cur_perm = perm_set(best_perms(bb),:);
        new_labels = prev_labels;
        label_set = unique(prev_labels); label_set(label_set == 1) = [];
        for ii = 1:length(cur_perm)
            new_labels(prev_labels==ii+1) = cur_perm(ii)+1;
        end
        su_spks = prev_clusts > 1;
        new_clusts = prev_clusts;
        new_clusts(su_spks) = cur_perm(prev_clusts(su_spks)-1) + 1;
        Clusters{probe_num}.spike_clusts = new_clusts;
        Clusters{probe_num}.cluster_labels = new_labels;
        map_to = cur_perm + 1; map_to(~ismember(map_to,label_set)) = [];
        Clusters{probe_num}.mean_spike(:,label_set) = Clusters{probe_num}.mean_spike(:,map_to);
        Clusters{probe_num}.std_spike(:,label_set) = Clusters{probe_num}.std_spike(:,map_to);
        Clusters{probe_num}.n_spks(label_set) = Clusters{probe_num}.n_spks(map_to);
        Clusters{probe_num}.refract = Clusters{probe_num}.refract(map_to-1,:);
        Clusters{probe_num}.Lratios = Clusters{probe_num}.Lratios(map_to-1);
        Clusters{probe_num}.iso_dists = Clusters{probe_num}.iso_dists(map_to-1);
        
        save(cur_data,'Clusters');
        
    end
end

