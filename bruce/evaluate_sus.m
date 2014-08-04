expt_nums = [81 85 86 87 88 89 91 92 93 95];

clust_dname = 'fin_aclust_data';

dir_prefix = '~';

mah_thresh = 2;
min_used_blocks = 10;

for ee = 1:length(expt_nums)
    data_dir = [dir_prefix '/Data/bruce/G0' num2str(expt_nums(ee))];
    cd(data_dir);
    
    load(sprintf('jbeG0%dExpts.mat',expt_nums(ee))); %load in Expts struct
    n_blocks = length(Expts);
    
    missing_blocks{ee} = [];
    for bb = 1:n_blocks
        fname = sprintf('Expt%dClusterTimes.mat',bb);
        if ~exist(fname,'file')
            missing_blocks{ee} = [missing_blocks{ee} bb];
        end
    end
    if expt_nums(ee) == 93
        missing_blocks{ee} = [missing_blocks{ee} 56 57];
    end
    
    load(clust_dname);
    
    used_blocks = setdiff(1:n_blocks,missing_blocks{ee});
    
    [cur_nblocks,cur_probes] = size(autoclust);
    aclust_mahals = nan(size(autoclust));
    aclust_mancodes = nan(size(autoclust));
    target_probes = zeros(cur_probes,1);
    for pp = 1:cur_probes
        target_probes(pp) = autoclust(used_blocks(1),pp).probenum;
        aclust_mancodes(used_blocks,pp) = [autoclust(used_blocks,pp).man_code];
        aclust_mahals(used_blocks,pp) = [autoclust(used_blocks,pp).mahal_d];
    end
    
    aclust_CellList = zeros(n_blocks,96);
    for pp = 1:length(target_probes)
        good_set = find(aclust_mahals(:,pp) > mah_thresh);
        if length(good_set) < min_used_blocks
            good_set = [];
        end
        aclust_CellList(good_set,target_probes(pp)) = pp;
    end
    cur_sus = any(aclust_CellList);
    all_is_su(ee,:) = false(1,96);
    all_is_su(ee,:) = cur_sus;
%     cur_avg_mahal = nanmean(aclust_mahals);
    use_set = find(cur_sus(target_probes));
    all_avg_mahal(ee,:) = nan(1,96);
    for pp = 1:length(use_set)
        cur_blocks = find(aclust_mahals(:,use_set(pp)) > mah_thresh);
        all_avg_mahal(ee,target_probes(use_set(pp))) = nanmean(aclust_mahals(cur_blocks,use_set(pp)));
    end
    
end

su_data.expt_nums = expt_nums;
su_data.is_su = all_is_su;
su_data.avg_mahal = all_avg_mahal;
su_data.mah_thresh = mah_thresh;
su_data.min_used_blocks = min_used_blocks;
%% SAVE SUMMARY SU DATA
sum_dataname = ['~/Analysis/bruce/summary_analysis/su_data'];
save(sum_dataname,'su_data');


%% PLOT AVERAGE WAVEFORMS FOR BEST AND WORST BLOCKS, ALONG WITH SCATTERPLOTS AND DENSITIES


