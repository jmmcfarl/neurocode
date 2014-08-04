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
    dsf = 40;
    filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),1);
    [~,~,Fs] = load_lfp_data(filename,dsf);
    alpha_lcf = 10;
    alpha_hcf = 20;
    [b_alpha,a_alpha] = butter(2,[alpha_lcf alpha_hcf]/(Fs/2));
    lgamm_lcf = 20;
    lgamm_hcf = 40;
    [b_lg,a_lg] = butter(2,[lgamm_lcf lgamm_hcf]/(Fs/2));
    hgamm_lcf = 40;
    hgamm_hcf = 90;
    [b_hg,a_hg] = butter(2,[hgamm_lcf hgamm_hcf]/(Fs/2));
    clear lfp_data lfp_alpha_amp lfp_lg_amp lfp_hg_amp
    for i = 1:96
        fprintf('Loading LFP %d of %d\n',i,96);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),i);
        [lfp_data(i,:),lfp_timestamps,Fs] = load_lfp_data(filename,dsf);
        lfp_alpha_amp(i,:) = filtfilt(b_alpha,a_alpha,lfp_data(i,:));
        lfp_alpha_amp(i,:) = abs(hilbert(lfp_alpha_amp(i,:)));
        lfp_lg_amp(i,:) = filtfilt(b_lg,a_lg,lfp_data(i,:));
        lfp_lg_amp(i,:) = abs(hilbert(lfp_lg_amp(i,:)));
        lfp_hg_amp(i,:) = filtfilt(b_hg,a_hg,lfp_data(i,:));
        lfp_hg_amp(i,:) = abs(hilbert(lfp_hg_amp(i,:)));
    end
   
    lfp_alpha_amp = zscore(lfp_alpha_amp')';
    lfp_lg_amp = zscore(lfp_lg_amp')';
    lfp_hg_amp = zscore(lfp_hg_amp')';
    
    %%
    buffer_win = 0.05;
    clear avg_gam_pow
    n_used_trials(ee) = length(completed_trials);
    spk_counts = nan(96,n_used_trials(ee));
    for i = 1:n_used_trials(ee)
        cur_set = find(lfp_timestamps > (Trial_starts(completed_trials(i)) + buffer_win) & lfp_timestamps < Trial_ends(completed_trials(i)));
        avg_alpha_pow(i,:) = mean(lfp_alpha_amp(:,cur_set),2);
        avg_lg_pow(i,:) = mean(lfp_lg_amp(:,cur_set),2);
        avg_hg_pow(i,:) = mean(lfp_hg_amp(:,cur_set),2);
        
        for su = 1:n_allunits
            cur_spk_times = Clusters{su}.times(Clusters{su}.times > Trial_starts(completed_trials(i)) & ...
                Clusters{su}.times < Trial_ends(completed_trials(i)));
            spk_counts(su,i) = length(cur_spk_times);
        end
    end
%     avg_alpha_pow = zscore(avg_alpha_pow);
%     avg_lg_pow = zscore(avg_lg_pow);
%     avg_hg_pow = zscore(avg_hg_pow);
    spk_rates = bsxfun(@rdivide,spk_counts,Trial_durs(completed_trials));

    %%
    unique_oris = unique(grating_oris);
    unique_oris(unique_oris==11.25) = []; %these orientations are out of place
    unique_sfs = unique(spatial_freqs);
    alpha_ori_avg(ee,:,:) = nan(96,length(unique_oris));
    alpha_sf_avg(ee,:,:) = nan(96,length(unique_sfs));
    lg_ori_avg(ee,:,:) = nan(96,length(unique_oris));
    lg_sf_avg(ee,:,:) = nan(96,length(unique_sfs));
    hg_ori_avg(ee,:,:) = nan(96,length(unique_oris));
    hg_sf_avg(ee,:,:) = nan(96,length(unique_sfs));
    mu_ori_avg(ee,:,:) = nan(96,length(unique_oris));
    mu_sf_avg(ee,:,:) = nan(96,length(unique_sfs));
    n_used_oris(ee,:) = nan(1,length(unique_oris));
    n_used_sfs(ee,:) = nan(1,length(unique_sfs));
    
    for i = 1:length(unique_oris)
        cur_set = find(grating_oris(completed_trials) == unique_oris(i));
        alpha_ori_avg(ee,:,i) = mean(avg_alpha_pow(cur_set,:));
        lg_ori_avg(ee,:,i) = mean(avg_lg_pow(cur_set,:));
        hg_ori_avg(ee,:,i) = mean(avg_hg_pow(cur_set,:));
        mu_ori_avg(ee,:,i) = mean(spk_rates(:,cur_set),2);
        n_used_oris(ee,i) = length(cur_set);
    end
    for i = 1:length(unique_sfs)
        cur_set = find(spatial_freqs(completed_trials) == unique_sfs(i));
        alpha_sf_avg(ee,:,i) = mean(avg_alpha_pow(cur_set,:));
        lg_sf_avg(ee,:,i) = mean(avg_lg_pow(cur_set,:));
        hg_sf_avg(ee,:,i) = mean(avg_hg_pow(cur_set,:));
        mu_sf_avg(ee,:,i) = mean(spk_rates(:,cur_set),2);
        n_used_sfs(ee,i) = length(cur_set);
    end           
end

%%
avg_alpha_ori_profile = squeeze(mean(alpha_ori_avg));
norm_alpha_ori_tuning = bsxfun(@rdivide,avg_alpha_ori_profile,max(avg_alpha_ori_profile,[],2));
avg_lg_ori_profile = squeeze(mean(lg_ori_avg));
norm_lg_ori_tuning = bsxfun(@rdivide,avg_lg_ori_profile,max(avg_lg_ori_profile,[],2));
avg_hg_ori_profile = squeeze(mean(hg_ori_avg));
norm_hg_ori_tuning = bsxfun(@rdivide,avg_hg_ori_profile,max(avg_hg_ori_profile,[],2));
avg_mu_ori_profile = squeeze(mean(mu_ori_avg));
norm_mu_ori_tuning = bsxfun(@rdivide,avg_mu_ori_profile,max(avg_mu_ori_profile,[],2));

avg_alpha_sf_profile = squeeze(mean(alpha_sf_avg));
norm_alpha_sf_tuning = bsxfun(@rdivide,avg_alpha_sf_profile,max(avg_alpha_sf_profile,[],2));
avg_lg_sf_profile = squeeze(mean(lg_sf_avg));
norm_lg_sf_tuning = bsxfun(@rdivide,avg_lg_sf_profile,max(avg_lg_sf_profile,[],2));
avg_hg_sf_profile = squeeze(mean(hg_sf_avg));
norm_hg_sf_tuning = bsxfun(@rdivide,avg_hg_sf_profile,max(avg_hg_sf_profile,[],2));
avg_mu_sf_profile = squeeze(mean(mu_sf_avg));
norm_mu_sf_tuning = bsxfun(@rdivide,avg_mu_sf_profile,max(avg_mu_sf_profile,[],2));

% cd ~/Data/bruce/7_15_12/G034
% save grating_avg_profiles avg_* norm_* unique_oris unique_sfs
%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
%%
ori_rad = degtorad(unique_oris);
full_alpha_profiles = [avg_alpha_ori_profile(:,end) avg_alpha_ori_profile];
full_lg_profiles = [avg_lg_ori_profile(:,end) avg_lg_ori_profile];
full_hg_profiles = [avg_hg_ori_profile(:,end) avg_hg_ori_profile];
full_oris = [-unique_oris(end) unique_oris];
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt)
            hold on
%             plot(full_oris,full_alpha_profiles(cur,:));
%             plot(full_oris,full_lg_profiles(cur,:),'r');
            plot(full_oris,full_hg_profiles(cur,:),'k');
            axis tight;
            xl = xlim();
            line(xl,[0 0],'color','k')
        end
        cnt = cnt + 1;
    end
end
%%
full_alpha_profiles = [avg_alpha_sf_profile];
full_lg_profiles = [avg_lg_sf_profile];
full_hg_profiles = [avg_hg_sf_profile];
full_sfs = [unique_sfs];
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt)
            hold on
            plot(full_sfs,full_alpha_profiles(cur,:));
            plot(full_sfs,full_lg_profiles(cur,:),'r');
            plot(full_sfs,full_hg_profiles(cur,:),'k');
            axis tight;
            xl = xlim();
            line(xl,[0 0],'color','k')
        end
        cnt = cnt + 1;
    end
end

%%
ori_rad = degtorad(unique_oris);
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt)
            plot(avg_mu_ori_profile(cur,:));
%             polar(ori_rad,avg_mu_ori_profile(cur,:)));
            axis tight; 
            xl = xlim();
            line(xl,[0 0],'color','k')
        end
        cnt = cnt + 1;
    end
end

%%
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt)
            plot(avg_gam_sf_profile(cur,:));
            axis tight; 
            xl = xlim();
            line(xl,[0 0],'color','k')
        end
        cnt = cnt + 1;
    end
end
%%
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt)
            plot(avg_mu_sf_profile(cur,:));
            axis tight; 
            xl = xlim();
            line(xl,[0 0],'color','k')
        end
        cnt = cnt + 1;
    end
end
