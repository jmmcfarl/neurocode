clear all
close all
cd ~/Data/bruce/7_15_12/

% cd G029/
% load ./jbeG029.em.mat
% em_data = Expt; clear Expt
% load ./CellList.mat
% load ./G029Expts.mat
% Expt_nu = [10:13]; %these are the grating expts

% cd G035/
% load ./jbeG035.em.mat
% em_data = Expt; clear Expt
% load ./CellList.mat
% load ./G035Expts.mat
% Expt_nu = [6 7 9]; %these are the grating expts

cd G034/
load ./jbeG034.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G034Expts.mat
Expt_nu = [26 27]; %these are the grating expts

n_allunits = 96;

single_units = find(CellList(1,:,1) > 0);
n_sus = length(single_units);
multi_units = setdiff(1:n_allunits,single_units);

%%
for ee = 1:length(Expt_nu)
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    
    
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
%     dt = 5e-3;
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
    spk_rates = bsxfun(@rdivide,spk_counts,Trial_durs(completed_trials));
    
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
    
%     [~,pref_ori_ind] = max(ori_avg,[],3);
%     [~,pref_sf_ind] = max(sf_avg,[],3);
%     for i = 1:n_allunits
%         cur_set = find(grating_oris(completed_trials) == unique_oris(pref_ori_ind(i)) & ...
%             spatial_freqs(completed_trials) == unique_sfs(pref_sf_ind(i)));
%         avg_binned(i,:) = mean(binned_spks(i,cur_set,:),2);
%     end
        
    %%
%     weights = squeeze(bsxfun(@rdivide,ori_avg(ee,:,:),sum(ori_avg(ee,:,:),3)));
%     unique_ori_rads = deg2rad(unique_oris);
%     for i = 1:n_allunits
%         [~,pref_ori(ee,i)] = max(ori_avg(ee,i,:));
% %         mean_ori(ee,i) = circ_mean(unique_ori_rads,weights(i,:));
% %         kappa_ori(ee,i) = circ_kappa(unique_ori_rads,weights(i,:));
%     end
%     pref_ori(ee,:) = unique_oris(pref_ori(ee,:));
    
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

cd ~/Data/bruce/7_15_12/
% cd G029/
cd G034/
save grating_mu_data *mu* unique* 

%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
ori_rad = degtorad(unique_oris);
full_profiles = [avg_mu_ori_profile(:,end) avg_mu_ori_profile];
full_oris = [-unique_oris(end) unique_oris];
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt)
            plot(full_oris,full_profiles(cur,:));
%             polar(ori_rad,avg_mu_ori_profile(cur,:)));
            axis tight; 
            xl = xlim();
            set(gca,'xtick',[-180 0 180])
            line(xl,[0 0],'color','k')
        end
        cnt = cnt + 1;
    end
end

