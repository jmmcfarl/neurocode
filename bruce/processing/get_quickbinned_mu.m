function [all_binned_mua] = get_quickbinned_mu(all_spk_times,all_clust_ids,all_t_axis,all_t_bin_edges,all_bin_edge_pts,clust_params)


n_probes = clust_params.n_probes;

all_binned_mua = nan(length(all_t_axis),n_probes);
all_mu_spk_times = cell(n_probes,1);
for cc = 1:n_probes
    cur_mua_inds = find(all_clust_ids{cc} >= 1);
        
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    all_mu_spk_times{cc} = cat(1,all_mu_spk_times{cc},all_spk_times{cc}(cur_mua_inds));
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end


