clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
n_expts = length(Expts);
for i = 1:n_expts
    load(sprintf('Expt%dClusterTimes.mat',i));
    
    n_trials(i) = length(Expts{i}.Trials);
    cur_trial_starts = [Expts{i}.Trials(:).TrialStart]/1e4;
    cur_trial_durs = [Expts{i}.Trials(:).dur]/1e4;
    cur_trial_ends = cur_trial_starts + cur_trial_durs;
    
    trial_avg_rates = nan(96,length(cur_trial_starts));
    for j = 1:n_trials(i)
        for cc = 1:96
            temp = sum(Clusters{cc}.times >= cur_trial_starts(j) & Clusters{cc}.times <= cur_trial_ends(j));
            trial_avg_rates(cc,j) = temp/cur_trial_durs(j);
        end
    end
    block_avg_rates(i,:) = nanmean(trial_avg_rates,2);
end
%%
for i = 1:n_expts
    load(sprintf('Expt%dClusterTimesDetails.mat',i));
    
    for cc = 1:96
            back_set = find(ClusterDetails{cc}.clst == 1);
            clust_set = find(ClusterDetails{cc}.clst == 2);
        if length(back_set) > 50 & length(clust_set) > 50
            
            
            back_cov = cov(ClusterDetails{cc}.xy(back_set,:));
            spk_X = ClusterDetails{cc}.xy(clust_set,:);
            mah_dist = spk_X*inv(back_cov);
            mah_dist = sum(mah_dist.*spk_X,2);
            avg_mah_dist(i,cc) = mean(mah_dist);
        else
            avg_mah_dist(i,cc) = nan;
        end
        %         plot(ClusterDetails{cc}.xy(:,1),ClusterDetails{cc}.xy(:,2),'.')
        %         hold on
        %
        %         plot(ClusterDetails{cc}.xy(clust_set,1),ClusterDetails{cc}.xy(clust_set,2),'r.')
        %         avg_mah_dist(i,cc)
        %         pause
        %         clf
    end
end
ov_mah_dist = nanmean(avg_mah_dist);
%%
save spk_stability_data block_avg_rates ov_mah_dist 