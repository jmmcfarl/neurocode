clear all
close all
cd ~/Data/bruce/7_15_12/

cd G029/
load ./jbeG029.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G029Expts.mat
Expt_nu = [10:13]; %these are the grating expts

% cd G035/
% load ./jbeG035.em.mat
% em_data = Expt; clear Expt
% load ./CellList.mat
% load ./G035Expts.mat
% Expt_nu = [6 7 9]; %these are the grating expts

n_allunits = 96;

%%
all_spk_counts = [];
all_oris = [];
all_sfs = [];
for ee = 1:length(Expt_nu)
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    
    single_units{ee} = find(CellList(Expt_nu(ee),:,1) > 0);
    n_sus = length(single_units);
    multi_units{ee} = setdiff(1:n_allunits,single_units{ee});
    
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    spatial_freqs = [Expts{Expt_nu(ee)}.Trials(:).sf];
    grating_oris = [Expts{Expt_nu(ee)}.Trials(:).or];
    se = [Expts{Expt_nu(ee)}.Trials(:).se]; %changes
    st = [Expts{Expt_nu(ee)}.Trials(:).st]; %changes
    rw = [Expts{Expt_nu(ee)}.Trials(:).rw]; %changes
%     nf = [Expts{Expt_nu(ee)}.Trials(:).nf]; %changes
    completed_trials = find(Trial_durs > 0.4);
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    
    %% weed out trials with blinks in them
    has_a_blink = zeros(size(completed_trials));
    for i = 1:length(completed_trials)
        cur_set = find(eye_ts > Trial_starts(completed_trials(i)) & eye_ts < Trial_ends(completed_trials(i)));
        if any(in_blink(cur_set))
            has_a_blink(i) = 1;
        end
    end
    completed_trials(has_a_blink==1) = [];
    
    %%
    beg_off = 0.025;
    dt = .05;
    
    Trial_starts = Trial_starts + beg_off;
    Trial_ends = Trial_starts + dt;

    n_used_trials(ee) = length(completed_trials);
    spk_counts = nan(n_allunits,n_used_trials(ee));
%     binned_spks = nan(n_allunits,n_used_trials(ee),0.4/dt+1);
    for i = 1:n_used_trials(ee)
%         cur_t_bin_edges = Trial_starts(completed_trials(i)):dt:(Trial_starts(completed_trials(i))+0.4);
        for su = 1:n_allunits
            cur_spk_times = Clusters{su}.times(Clusters{su}.times > Trial_starts(completed_trials(i)) & ...
                Clusters{su}.times < Trial_ends(completed_trials(i)));
%             binned_spks(su,i,:) = histc(cur_spk_times,cur_t_bin_edges);
            spk_counts(su,i) = length(cur_spk_times);
        end
    end
%     spk_rates = bsxfun(@rdivide,spk_counts,Trial_durs(completed_trials));
    spk_rates = spk_counts;
all_spk_counts = [all_spk_counts; spk_counts'];
    all_oris = [all_oris; grating_oris(completed_trials)'];
    all_sfs = [all_sfs; spatial_freqs(completed_trials)'];
    
    %%
    unique_oris = unique(grating_oris);
    unique_oris(unique_oris==11.25) = []; %these orientations are out of place
    unique_sfs = unique(spatial_freqs);
    ori_avg(ee,:,:) = nan(n_allunits,length(unique_oris));
    sf_avg(ee,:,:) = nan(n_allunits,length(unique_sfs));
    n_used_oris(ee,:) = nan(1,length(unique_oris));
    n_used_sfs(ee,:) = nan(1,length(unique_sfs));
    
    for i = 1:length(unique_oris)
        cur_set = find(grating_oris(completed_trials) == unique_oris(i));
        ori_avg(ee,:,i) = mean(spk_rates(:,cur_set),2);
        n_used_oris(ee,i) = length(cur_set);
    end
    for i = 1:length(unique_sfs)
        cur_set = find(spatial_freqs(completed_trials) == unique_sfs(i));
        sf_avg(ee,:,i) = mean(spk_rates(:,cur_set),2);
        n_used_sfs(ee,i) = length(cur_set);
    end
    
end
% ori_set1 = find(unique_oris <= 0);
% ori_set2 = find(unique_oris > 0);
% pref_ori_mod = mod(pref_ori,180);
% mean_dir1 = mean(avg_ori_profile(:,ori_set1),2);
% mean_dir2 = mean(avg_ori_profile(:,ori_set2),2);
% ds_index = (mean_dir1-mean_dir2)./(mean_dir1+mean_dir2);

%%
avg_mu_sf_profile = squeeze(mean(sf_avg));
avg_mu_ori_profile = squeeze(mean(ori_avg));
norm_mu_ori_tuning = bsxfun(@rdivide,avg_mu_ori_profile,max(avg_mu_ori_profile,[],2));
[peak_mu_ori_rate,pref_mu_ori_ind] = max(avg_mu_ori_profile,[],2);
pref_mu_ori = unique_oris(pref_mu_ori_ind);
antipref_ori = pref_mu_ori + 180;
antipref_ori(antipref_ori > 180) = antipref_ori(antipref_ori > 180) - 360;
for i = 1:n_allunits
    cur_loc = find(unique_oris == antipref_ori(i));
    antipref_rate(i) = avg_mu_ori_profile(i,cur_loc);
end
mu_ds_index = (peak_mu_ori_rate - antipref_rate')./(peak_mu_ori_rate + antipref_rate');
pref_mu_ori = mod(pref_mu_ori,180);


%%
for t = 1:96
    X = 0.5*(cos(2*degtorad(all_oris-pref_mu_ori(t)))+1);
    Y = all_spk_counts(:,t);
    unique_spk_cnts = unique(Y);
    spikebins = [];
    for i = 2:length(unique_spk_cnts)
        cur_set = find(Y == unique_spk_cnts(i));
        spikebins = [spikebins; repmat(cur_set,unique_spk_cnts(i),1)];
    end
    spikebins = sort(spikebins);
    
    [fitp,grad] = GLMsolve_jmm(X, spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
    beta(t,:) = fitp.k;
    offset_rate(t) = log(1+exp(beta(t,2)))/dt;
    max_rate(t) = log(1+exp(beta(t,2) + beta(t,1)))/dt;
    [b(t,:),dev(t),stats(t)] = glmfit(X,Y,'poisson');
    
%     
%     pred_X = 0.5*(cos(2*degtorad(unique_oris-pref_mu_ori(t)))+1);
%     pred_r = beta(t,1)*pred_X + beta(t,2);
%     pred_r = log(1+exp(pred_r));
%     
%     figure
%     plot(unique_oris,avg_mu_ori_profile(t,:))
%     hold on
%     plot(unique_oris,pred_r,'r')
%     pause
%     close all
end

%%
cd ~/Data/bruce/7_15_12/G029
save grating_glm_fits beta offset_rate max_rate stats 