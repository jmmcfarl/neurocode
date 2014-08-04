clear all
close all
cd
addpath('/Users/James/Data/bruce/7_15_12/')

cd /Users/James/Data/bruce/7_15_12/G029/

load ./jbeG029.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G029Expts.mat
load ./eye_calibration_data

dsf = 10;
n_pos_bins = 15;
x_range_edges = linspace(-0.185,0.865,n_pos_bins+1);
y_range_edges = linspace(-.96,0.1,n_pos_bins+1);
s_range_edges = linspace(-0.55,0.55,n_pos_bins+1);
%%
Expt_nu = [30 31]; %these are the grating expts
n_allunits = 96;
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
    line_oris = [Expts{Expt_nu(ee)}.Trials(:).or];
    xo = [Expts{Expt_nu(ee)}.Trials(:).xo];
    yo = [Expts{Expt_nu(ee)}.Trials(:).yo];
    sO = [Expts{Expt_nu(ee)}.Trials(:).sO]; %changes
    ar = [Expts{Expt_nu(ee)}.Trials(:).ar]; %changes
    %     completed_trials = find(Trial_durs > 0.4);
    completed_trials = 1:length(Trial_durs); %use all trials
    
    %% GET NEEDED EYE TRACKING DATA
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,1:2),eye_vals(:,3:4),eye_ts);
    
    [sac_data,in_sac,eye_speed] = get_saccades(eye_vals(:,1:2),eye_vals(:,3:4),eye_ts);
    
    [fixation_data,in_fixation] = parse_fixations(in_blink,in_sac,eye_ts);
    
    in_window = (abs(eye_vals(:,1)) < 1 & abs(eye_vals(:,2)) < 1);
    
    used_points = in_window' == 1 & in_fixation == 1;
    
    eye_dt = mode(diff(eye_ts));
    eye_fs = 1/eye_dt;
    eye_fsd = eye_fs/dsf;
    eye_vals_d = downsample(eye_vals,dsf);
    eye_ts_d = downsample(eye_ts,dsf);
    in_blink_d = downsample(in_blink,dsf);
    used_points_d = downsample(used_points,dsf);
    
    %% GET LFP GAMMA DATA
%     dsf = 40;
%     gam_lcf = 80;
%     gam_hcf = 20;
    
    dsf = 30;
    gam_lcf = 50;
    gam_hcf = 100;
    pow_smooth = 1/eye_fsd;
    clear lfp_data lfp_gam_amp
    for i = 1:96
        fprintf('Loading LFP %d of %d\n',i,96);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),i);
        [lfp_data(i,:),lfp_timestamps,Fs] = load_lfp_data(filename,dsf);        
        [b,a] = butter(2,[gam_lcf gam_hcf]/(Fs/2));
        lfp_gam_amp(i,:) = filtfilt(b,a,lfp_data(i,:));
        lfp_gam_amp(i,:) = abs(hilbert(lfp_gam_amp(i,:)));
        lfp_gam_amp(i,:) = smooth(lfp_gam_amp(i,:),round(pow_smooth*Fs));
    end
    
%     lfp_timestamps_interp = lfp_timestamps(1):1/eye_fsd:lfp_timestamps(end);
    lfp_gam_amp_interp = interp1(lfp_timestamps,lfp_gam_amp',eye_ts_d)';
    
    %% BIN SPIKES
    binned_spks = nan(n_allunits,size(eye_vals_d,1));
    for su = 1:n_allunits
        cur_spk_times = Clusters{su}.times;
        binned_spks(su,:) = hist(cur_spk_times,eye_ts_d);
    end
    avg_rates(ee,:) = mean(binned_spks(:,used_points_d),2);
    
    %%
    lfp_gam_amp_interp = bsxfun(@minus,lfp_gam_amp_interp,nanmean(lfp_gam_amp_interp,2));
    lfp_gam_amp_interp = bsxfun(@rdivide,lfp_gam_amp_interp,nanstd(lfp_gam_amp_interp,[],2));
    binned_spks = bsxfun(@minus,binned_spks,nanmean(binned_spks,2));
    binned_spks = bsxfun(@rdivide,binned_spks,nanstd(binned_spks,[],2));
    
    %% COMPUTE STIMULUS PARAMETER VECTORS
    interp_trial_starts = round(interp1(eye_ts_d,1:length(eye_ts_d),Trial_starts));
    interp_trial_stops = round(interp1(eye_ts_d,1:length(eye_ts_d),Trial_ends));
    unique_oris = unique(line_oris);
    unique_ars = unique(ar);
    
    bar_ori_vec = nan(size(eye_vals_d,1),1);
    bar_ar_vec = nan(size(eye_vals_d,1),1);
    bar_so_vec = nan(size(eye_vals_d,1),1);
    bar_xo_vec = nan(size(eye_vals_d,1),1);
    bar_yo_vec = nan(size(eye_vals_d,1),1);
    for i = 1:length(completed_trials)
        cur_inds = interp_trial_starts(completed_trials(i)):interp_trial_stops(completed_trials(i));
        bar_ori_vec(cur_inds) = line_oris(completed_trials(i));
        bar_ar_vec(cur_inds) = ar(completed_trials(i));
        bar_so_vec(cur_inds) = sO(completed_trials(i));
        bar_xo_vec(cur_inds) = xo(completed_trials(i));
        bar_yo_vec(cur_inds) = yo(completed_trials(i));
    end
    
    bar_ori_vec(used_points_d==0) = nan;
    bar_ar_vec(used_points_d==0) = nan;
    
    %%  
    min_samps = 5;
    delay = round(0.05*eye_fsd);
    para_mu_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    orth_mu_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    para_gam_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    orth_gam_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    for i = 1:length(unique_oris)
        cur_set = find(bar_ori_vec == unique_oris(i) & bar_ar_vec == unique_ars(1));
        [bins,indx] = histc(bar_so_vec(cur_set),s_range_edges);
        
        n_ov_para(i) = length(cur_set);
        for j = 1:n_pos_bins
            cur_inds = cur_set(indx==j);
            cur_inds = cur_inds + delay;
            cur_inds(cur_inds > length(bar_so_vec)) = [];
            if length(cur_inds) > min_samps
                para_mu_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
                para_gam_profile(ee,:,i,j) = nanmean(lfp_gam_amp_interp(:,cur_inds),2);
            end
            n_para_bins(i,j) = length(cur_inds);
            avg_para_xbin(i,j) = nanmean(bar_xo_vec(cur_inds));
            avg_para_ybin(i,j) = nanmean(bar_yo_vec(cur_inds));
            avg_para_sbin(i,j) = nanmean(bar_so_vec(cur_inds));
        end
        
        cur_set = find(bar_ori_vec == unique_oris(i) & bar_ar_vec == unique_ars(2));
        n_ov_orth(i) = length(cur_set);
        [bins,indx] = histc(bar_so_vec(cur_set),s_range_edges);
        for j = 1:n_pos_bins
            cur_inds = cur_set(indx==j);
            cur_inds = cur_inds + delay;
            cur_inds(cur_inds > length(bar_so_vec)) = [];
            if length(cur_inds) > min_samps
                orth_mu_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
                orth_gam_profile(ee,:,i,j) = nanmean(lfp_gam_amp_interp(:,cur_inds),2);
            end
            n_orth_bins(i,j) = length(cur_inds);
            avg_orth_xbin(i,j) = nanmean(bar_xo_vec(cur_inds));
            avg_orth_ybin(i,j) = nanmean(bar_yo_vec(cur_inds));
            avg_orth_sbin(i,j) = nanmean(bar_so_vec(cur_inds)); 
        end        
    end
    
end

%%
ov_para_mu_profile = squeeze(nanmean(para_mu_profile));
ov_orth_mu_profile = squeeze(nanmean(orth_mu_profile));
ov_para_gam_profile = squeeze(nanmean(para_gam_profile));
ov_orth_gam_profile = squeeze(nanmean(orth_gam_profile));
ov_avg_rates = mean(avg_rates);

%%
cd ~/Data/bruce/7_15_12/G029
save hgam_mu_1dnoise_profile unique* para_* orth_* avg_* eye_fsd
%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

%%
cur_ori = 4;
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt);hold on  
            plot(avg_para_sbin(cur_ori,:),squeeze(ov_para_gam_profile(cur,cur_ori,:)));
hold on
            plot(avg_orth_sbin(cur_ori,:),squeeze(ov_orth_gam_profile(cur,cur_ori,:)),'r');
            xlim([-0.6 0.6]);%ylim([-0.7 -0.1])
%             line([0.1 0.7],[-0.4 -0.4],'color','k')
%             line([0.4 0.4],[-0.7 -0.1],'color','k')
        end
        cnt = cnt + 1;
    end
end

%%
mu_ov_avg_profile = cat(2,ov_para_mu_profile,ov_orth_mu_profile);
ov_x_profile = [avg_para_xbin; avg_orth_xbin];
ov_y_profile = [avg_para_ybin; avg_orth_ybin];
ov_s_profile = [avg_para_sbin; avg_orth_sbin];
mu_max_profiles = max(mu_ov_avg_profile,[],3);
mu_min_profiles = min(mu_ov_avg_profile,[],3);
avg_bar_oris = [unique_oris(:); unique_oris(:)];
avg_bar_ars = [zeros(4,1); 10*ones(4,1)];
for i = 1:length(unique_oris)
    used_set = find(avg_bar_oris == unique_oris(i));
    mod_index(i,:) = (max(mu_max_profiles(:,used_set),[],2) - min(mu_min_profiles(:,used_set),[],2));
end
[bar_ori_maxrates,mu_pref_bar_oris] = max(mod_index);
mu_pref_bar_oris = unique_oris(mu_pref_bar_oris);
for i = 1:96
    
    cur_set = find(avg_bar_ars < 1 & avg_bar_oris == mu_pref_bar_oris(i));
    cur_par_x = ov_x_profile(cur_set,:);
    cur_par_y = ov_y_profile(cur_set,:);
    cur_par_so = ov_s_profile(cur_set,:);
    mu_parallel_profile(i,:) =  squeeze(mu_ov_avg_profile(i,cur_set,:));
    
    cur_set = find(avg_bar_ars > 1 & avg_bar_oris == mu_pref_bar_oris(i));
    cur_orth_x = ov_x_profile(cur_set,:);
    cur_orth_y = ov_y_profile(cur_set,:);
    cur_orth_so = ov_s_profile(cur_set,:);
    mu_orth_profile(i,:) = squeeze(mu_ov_avg_profile(i,cur_set,:));
            
    par_prof_norm = mu_parallel_profile(i,:) - min(mu_parallel_profile(i,:));
    par_prof_norm = par_prof_norm/nansum(par_prof_norm);
    orth_prof_norm = mu_orth_profile(i,:) - min(mu_orth_profile(i,:));
    orth_prof_norm = orth_prof_norm/nansum(orth_prof_norm);
    
    init_mean_x(i) = 0.5*nansum(cur_par_x.*par_prof_norm) + 0.5*nansum(cur_orth_x.*orth_prof_norm);
    init_mean_y(i) = 0.5*nansum(cur_par_y.*par_prof_norm) + 0.5*nansum(cur_orth_y.*orth_prof_norm);
    
    init_mean_par_s(i) = nansum(cur_par_so.*par_prof_norm);
    init_mean_orth_s(i) = nansum(cur_orth_so.*orth_prof_norm);
    
    init_std_par(i) = sqrt(nansum((cur_par_so - init_mean_par_s(i)).^2.*par_prof_norm));
    init_std_orth(i) = sqrt(nansum((cur_orth_so - init_mean_orth_s(i)).^2.*orth_prof_norm));

    
    beta0(1) = min(mu_parallel_profile(i,:)); beta0(2) = range(mu_parallel_profile(i,:));
    beta0(3) = init_mean_par_s(i); beta0(4) = init_std_par(i)^2;
    [mu_gfit_par(i,:),~,J,~,mu_par_mse(i)] = nlinfit(cur_par_so,mu_parallel_profile(i,:),@(X,cur_para_so) my_gaussfun(X,cur_par_so),beta0);
    mu_par_Jcond(i) = cond(J);
    mu_gpred_par(i,:) = my_gaussfun(mu_gfit_par(i,:),cur_par_so);
    mu_dom_par(i) = (max(mu_gpred_par(i,:)) - min(mu_gpred_par(i,:)));
    
    beta0(1) = min(mu_orth_profile(i,:)); beta0(2) = range(mu_orth_profile(i,:));
    beta0(3) = init_mean_orth_s(i); beta0(4) = init_std_orth(i)^2;
    [mu_gfit_orth(i,:),~,J,~,mu_orth_mse(i)] = nlinfit(cur_orth_so,mu_orth_profile(i,:),@(X,cur_orth_so) my_gaussfun(X,cur_orth_so),beta0);
    mu_orth_Jcond(i) = cond(J);
    mu_gpred_orth(i,:) = my_gaussfun(mu_gfit_orth(i,:),cur_orth_so);
    mu_dom_orth(i) = (max(mu_gpred_orth(i,:)) - min(mu_gpred_orth(i,:)));
    
    [~,parloc] = min(abs(cur_par_so-mu_gfit_par(i,3)));
    [~,orthloc] = min(abs(cur_orth_so-mu_gfit_orth(i,3)));    
    mu_mean_x(i) = 0.5*cur_par_x(parloc) + 0.5*cur_orth_x(orthloc);
    mu_mean_y(i) = 0.5*cur_par_y(parloc) + 0.5*cur_orth_y(orthloc);
    
%     beta0(1) = 0; beta0(2) = max(orth_profile(i,:));
%     beta0(3) = mean_orth_s(i); beta0(4) = std_orth(i)^2;
%     gfit_orth(i,:) = nlinfit(cur_orth_so,orth_profile(i,:),@(X,cur_orth_so) my_gaussfun(X,cur_orth_so),beta0);
end

%%
gam_ov_avg_profile = cat(2,ov_para_gam_profile,ov_orth_gam_profile);
ov_x_profile = [avg_para_xbin; avg_orth_xbin];
ov_y_profile = [avg_para_ybin; avg_orth_ybin];
ov_s_profile = [avg_para_sbin; avg_orth_sbin];
gam_max_profiles = max(gam_ov_avg_profile,[],3);
gam_min_profiles = min(gam_ov_avg_profile,[],3);
avg_bar_oris = [unique_oris(:); unique_oris(:)];
avg_bar_ars = [zeros(4,1); 10*ones(4,1)];
for i = 1:length(unique_oris)
    used_set = find(avg_bar_oris == unique_oris(i));
    mod_index(i,:) = (max(gam_max_profiles(:,used_set),[],2) - min(gam_min_profiles(:,used_set),[],2));
end
[bar_ori_maxrates,gam_pref_bar_oris] = max(mod_index);
gam_pref_bar_oris = unique_oris(gam_pref_bar_oris);
for i = 1:96
    
    cur_set = find(avg_bar_ars < 1 & avg_bar_oris == gam_pref_bar_oris(i));
    cur_par_x = ov_x_profile(cur_set,:);
    cur_par_y = ov_y_profile(cur_set,:);
    cur_par_so = ov_s_profile(cur_set,:);
    gam_parallel_profile(i,:) =  squeeze(gam_ov_avg_profile(i,cur_set,:));
    
    cur_set = find(avg_bar_ars > 1 & avg_bar_oris == gam_pref_bar_oris(i));
    cur_orth_x = ov_x_profile(cur_set,:);
    cur_orth_y = ov_y_profile(cur_set,:);
    cur_orth_so = ov_s_profile(cur_set,:);
    gam_orth_profile(i,:) = squeeze(gam_ov_avg_profile(i,cur_set,:));
            
    par_prof_norm = gam_parallel_profile(i,:) - min(gam_parallel_profile(i,:));
    par_prof_norm = par_prof_norm/nansum(par_prof_norm);
    orth_prof_norm = gam_orth_profile(i,:) - min(gam_orth_profile(i,:));
    orth_prof_norm = orth_prof_norm/nansum(orth_prof_norm);
    
    init_mean_x(i) = 0.5*nansum(cur_par_x.*par_prof_norm) + 0.5*nansum(cur_orth_x.*orth_prof_norm);
    init_mean_y(i) = 0.5*nansum(cur_par_y.*par_prof_norm) + 0.5*nansum(cur_orth_y.*orth_prof_norm);
    
    init_mean_par_s(i) = nansum(cur_par_so.*par_prof_norm);
    init_mean_orth_s(i) = nansum(cur_orth_so.*orth_prof_norm);
    
    init_std_par(i) = sqrt(nansum((cur_par_so - init_mean_par_s(i)).^2.*par_prof_norm));
    init_std_orth(i) = sqrt(nansum((cur_orth_so - init_mean_orth_s(i)).^2.*orth_prof_norm));

    
    beta0(1) = min(gam_parallel_profile(i,:)); beta0(2) = range(gam_parallel_profile(i,:));
    beta0(3) = init_mean_par_s(i); beta0(4) = init_std_par(i)^2;
    [gam_gfit_par(i,:),~,J,~,gam_par_mse(i)] = nlinfit(cur_par_so,gam_parallel_profile(i,:),@(X,cur_para_so) my_gaussfun(X,cur_par_so),beta0);
    gam_par_Jcond(i) = cond(J);
    gam_gpred_par(i,:) = my_gaussfun(gam_gfit_par(i,:),cur_par_so);
    gam_dom_par(i) = (max(gam_gpred_par(i,:)) - min(gam_gpred_par(i,:)));
    
    beta0(1) = min(gam_orth_profile(i,:)); beta0(2) = range(gam_orth_profile(i,:));
    beta0(3) = init_mean_orth_s(i); beta0(4) = init_std_orth(i)^2;
    [gam_gfit_orth(i,:),~,J,~,gam_orth_mse(i)] = nlinfit(cur_orth_so,gam_orth_profile(i,:),@(X,cur_orth_so) my_gaussfun(X,cur_orth_so),beta0);
    gam_orth_Jcond(i) = cond(J);
    gam_gpred_orth(i,:) = my_gaussfun(gam_gfit_orth(i,:),cur_orth_so);
    gam_dom_orth(i) = (max(gam_gpred_orth(i,:)) - min(gam_gpred_orth(i,:)));
    
    [~,parloc] = min(abs(cur_par_so-gam_gfit_par(i,3)));
    [~,orthloc] = min(abs(cur_orth_so-gam_gfit_orth(i,3)));    
    gam_mean_x(i) = 0.5*cur_par_x(parloc) + 0.5*cur_orth_x(orthloc);
    gam_mean_y(i) = 0.5*cur_par_y(parloc) + 0.5*cur_orth_y(orthloc);
    
%     beta0(1) = 0; beta0(2) = max(orth_profile(i,:));
%     beta0(3) = mean_orth_s(i); beta0(4) = std_orth(i)^2;
%     gfit_orth(i,:) = nlinfit(cur_orth_so,orth_profile(i,:),@(X,cur_orth_so) my_gaussfun(X,cur_orth_so),beta0);
end

%%
mu_std_par = sqrt(mu_gfit_par(:,4));
mu_std_orth = sqrt(mu_gfit_orth(:,4));
gam_std_par = sqrt(gam_gfit_par(:,4));
gam_std_orth = sqrt(gam_gfit_orth(:,4));

%% FOR MU
bad_fits = find(mu_par_Jcond > 1e3 | mu_orth_Jcond > 1e3);
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt);hold on  
            if ~ismember(cur,bad_fits)
            ellipse(mu_std_orth(cur),mu_std_par(cur),degtorad(mu_pref_bar_oris(cur)),...
                mu_mean_x(cur),mu_mean_y(cur));
            plot(mu_mean_x(cur),mu_mean_y(cur),'.','markersize',20);
            xlim([0.1 0.7]);ylim([-0.7 -0.1])
            line([0.1 0.7],[-0.4 -0.4],'color','k')
            line([0.4 0.4],[-0.7 -0.1],'color','k')
            end
        end
        cnt = cnt + 1;
    end
end

%% FOR GAM
bad_fits = find(gam_par_Jcond > 1e3 | gam_orth_Jcond > 1e3);
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt);hold on  
            if ~ismember(cur,bad_fits)
            ellipse(gam_std_orth(cur),gam_std_par(cur),degtorad(gam_pref_bar_oris(cur)),...
                gam_mean_x(cur),gam_mean_y(cur));
            plot(gam_mean_x(cur),gam_mean_y(cur),'.','markersize',20);
            xlim([0.1 0.7]);ylim([-0.7 -0.1])
            line([0.1 0.7],[-0.4 -0.4],'color','k')
            line([0.4 0.4],[-0.7 -0.1],'color','k')
            end
        end
        cnt = cnt + 1;
    end
end