clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
% loc_glob_blocks = [8 15 20 27 33 39 49] - 6;
forty_blocks = [7 13 17 19 23 26 29 31 36 38 48] - 6;
repeat_inds = [1:9]*1e3;
all_expt_id = [];
for bb = 1:length(forty_blocks)
    n_trials(bb) = length(Expts{forty_blocks(bb)}.Trials);
    trial_start_times{bb} = [Expts{forty_blocks(bb)}.Trials(:).Start]/1e4;
    trial_stop_times{bb} = [Expts{forty_blocks(bb)}.Trials(:).End]/1e4;
    trial_seof{bb} = [Expts{forty_blocks(bb)}.Trials(:).seof];
    trial_completed{bb} = [Expts{forty_blocks(bb)}.Trials(:).Result];
    trial_durs{bb} = trial_stop_times{bb} - trial_start_times{bb};
    all_expt_id = [all_expt_id forty_blocks(bb)*ones(1,n_trials(bb))];
end
all_trial_start_times = cell2mat(trial_start_times);
all_trial_stop_times = cell2mat(trial_stop_times);
all_trial_seof = cell2mat(trial_seof);
all_trial_completed = cell2mat(trial_completed);
all_trial_durs = cell2mat(trial_durs);

%%
stim_fs = 1e4/117.5;
use_clusters = 1:96;
bin_size = 0.005;
T = 3.75;
bin_edges = 0:bin_size:T;
nbins = length(bin_edges);
spk_smwin = 2;
for ss = 1:length(repeat_inds)
    all_binned_spks{ss} = [];
    all_smoothed_spks{ss} = [];
end
for bb = 1:length(forty_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(forty_blocks));
    load(sprintf('Expt%dClusterTimes.mat',forty_blocks(bb)));
    
    cur_trials = find(all_expt_id==forty_blocks(bb));
    cur_trials = cur_trials(all_trial_completed(cur_trials)==1);
    
    for ss = 1:length(repeat_inds)
        cur_trial_set = cur_trials(all_trial_seof(cur_trials)==repeat_inds(ss));
        
        binned_spks = nan(length(cur_trial_set),length(cur_t_axis),length(use_clusters));
        smoothed_spks = nan(length(cur_trial_set),length(cur_t_axis),length(use_clusters));
        for tt = 1:length(cur_trial_set)
            cur_t_axis = all_trial_start_times(cur_trial_set(tt)):bin_size:all_trial_stop_times(cur_trial_set(tt));
            cur_t_axis(nbins+1:end) = [];
            for cc = 1:length(use_clusters)
                binned_spks(tt,:,cc) = histc(Clusters{use_clusters(cc)}.times,cur_t_axis);
                smoothed_spks(tt,:,cc) = jmm_smooth_1d_cor(binned_spks(tt,:,cc),spk_smwin);
            end
%             smoothed_spks(tt,:,:) = zscore(squeeze(smoothed_spks(tt,:,:)));
        end
        all_binned_spks{ss} = cat(1,all_binned_spks{ss},binned_spks);
        all_smoothed_spks{ss} = cat(1,all_smoothed_spks{ss},smoothed_spks);
        
    end
    
end
%%
for ss =1:length(repeat_inds)
    avg_binned(ss,:,:) = squeeze(nanmean(all_binned_spks{ss}));
    sem_binned(ss,:,:) = squeeze(nanstd(all_binned_spks{ss}))/sqrt(sum(~isnan(all_binned_spks{ss}(:,1,1))));
    avg_smoothed(ss,:,:) = squeeze(nanmean(all_smoothed_spks{ss}));
    sem_smoothed(ss,:,:) = squeeze(nanstd(all_smoothed_spks{ss}))/sqrt(sum(~isnan(all_smoothed_spks{ss}(:,1,1))));
end
avg_binned = avg_binned/bin_size;
sem_binned = sem_binned/bin_size;
avg_smoothed = avg_smoothed/bin_size;
sem_smoothed = sem_smoothed/bin_size;
%%
stim_fs = 1e4/117.5;
stim_times = (160:160:320)/stim_fs;
stim1 = 4;
stim2 = 5;
stim3 = 6;
for cur_el = 1:96;
    cur_el
    shadedErrorBar(bin_edges,squeeze(avg_smoothed(stim1,:,cur_el)),squeeze(sem_smoothed(stim1,:,cur_el)),{'color','b'});
    hold on
    shadedErrorBar(bin_edges,squeeze(avg_smoothed(stim2,:,cur_el)),squeeze(sem_smoothed(stim2,:,cur_el)),{'color','r'});
    shadedErrorBar(bin_edges,squeeze(avg_smoothed(stim3,:,cur_el)),squeeze(sem_smoothed(stim3,:,cur_el)),{'color','k'});
%     shadedErrorBar(bin_edges,squeeze(avg_binned(stim1,:,cur_el)),squeeze(sem_binned(stim1,:,cur_el)),{'color','b'});
%     hold on
%     shadedErrorBar(bin_edges,squeeze(avg_binned(stim2,:,cur_el)),squeeze(sem_binned(stim2,:,cur_el)),{'color','r'});
%     shadedErrorBar(bin_edges,squeeze(avg_binned(stim3,:,cur_el)),squeeze(sem_binned(stim3,:,cur_el)),{'color','k'});
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    xlim([0.5 1.5])
    pause
    clf
end