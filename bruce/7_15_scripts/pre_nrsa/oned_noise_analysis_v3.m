clear all
close all
addpath('~/Data/bruce/7_15_12/')
cd ~/Data/bruce/7_15_12/G029/

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
Expt_nu = [30 31]; %these are the ond noise expts
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

    %% BIN SPIKES
    binned_spks = nan(n_allunits,size(eye_vals_d,1));
    for su = 1:n_allunits
        cur_spk_times = Clusters{su}.times;
        binned_spks(su,:) = hist(cur_spk_times,eye_ts_d);
    end
    avg_rates(ee,:) = mean(binned_spks(:,used_points_d),2);
    
    %% COMPUTE STIMULUS PARAMETER VECTORS
    interp_trial_starts = round(interp1(eye_ts_d,1:length(eye_ts_d),Trial_starts));
    interp_trial_stops = round(interp1(eye_ts_d,1:length(eye_ts_d),Trial_ends));
    
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
    
    unique_oris = unique(line_oris);
    unique_ars = unique(ar);

    %%  
    min_samps = 5;
    delay = round(0.05*eye_fsd);
    para_avg_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    orth_avg_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    for i = 1:length(unique_oris)
        cur_set = find(bar_ori_vec == unique_oris(i) & bar_ar_vec == unique_ars(1));
        [bins,indx] = histc(bar_so_vec(cur_set),s_range_edges);
        
        n_ov_para(i) = length(cur_set);
        for j = 1:n_pos_bins
            cur_inds = cur_set(indx==j);
            cur_inds = cur_inds + delay;
            cur_inds(cur_inds > length(bar_so_vec)) = [];
            if length(cur_inds) > min_samps
                para_avg_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
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
                orth_avg_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
            end
            n_orth_bins(i,j) = length(cur_inds);
            avg_orth_xbin(i,j) = nanmean(bar_xo_vec(cur_inds));
            avg_orth_ybin(i,j) = nanmean(bar_yo_vec(cur_inds));
            avg_orth_sbin(i,j) = nanmean(bar_so_vec(cur_inds));
        end        
    end    
end

%%
% cd ~/Data/bruce/7_15_12/G029
% save oned_noise_mapping_v3 *profile unique* *bin* 

%%
ov_para_profile = squeeze(nanmean(para_avg_profile));
ov_orth_profile = squeeze(nanmean(orth_avg_profile));
ov_avg_rates = mean(avg_rates);
%%
ov_avg_profile = cat(2,ov_para_profile,ov_orth_profile);
ov_x_profile = [avg_para_xbin; avg_orth_xbin];
ov_y_profile = [avg_para_ybin; avg_orth_ybin];
ov_s_profile = [avg_para_sbin; avg_orth_sbin];
max_profiles = max(ov_avg_profile,[],3);
min_profiles = min(ov_avg_profile,[],3);
avg_bar_oris = [unique_oris(:); unique_oris(:)];
avg_bar_ars = [zeros(4,1); 10*ones(4,1)];
for i = 1:length(unique_oris)
    used_set = find(avg_bar_oris == unique_oris(i));
    mod_index(i,:) = (max(max_profiles(:,used_set),[],2) - min(min_profiles(:,used_set),[],2));
end
[bar_ori_maxrates,pref_bar_oris] = max(mod_index);
pref_bar_oris = unique_oris(pref_bar_oris);

%%
% AR ORI
% 0 -45
% 0 0
% 0 45
% 0 90
% 10 -45
% 10 0
% 10 45
% 10 90

vec_angles = [
    135
    180
    -135
    -90
    45
    90
    135
    180];

xcent = 0.34; ycent = -0.43;

for i = 1:96
    
    cur_set_par = find(avg_bar_ars < 1 & avg_bar_oris == pref_bar_oris(i));
    cur_par_x = ov_x_profile(cur_set_par,:);
    cur_par_y = ov_y_profile(cur_set_par,:);
    cur_par_so = ov_s_profile(cur_set_par,:);
    parallel_profile(i,:) =  squeeze(ov_avg_profile(i,cur_set_par,:));
    parallel_x(i,:) = cur_par_x;
    parallel_y(i,:) = cur_par_y;
    
    cur_set_orth = find(avg_bar_ars > 1 & avg_bar_oris == pref_bar_oris(i));
    cur_orth_x = ov_x_profile(cur_set_orth,:);
    cur_orth_y = ov_y_profile(cur_set_orth,:);
    cur_orth_so = ov_s_profile(cur_set_orth,:);
    orth_profile(i,:) = squeeze(ov_avg_profile(i,cur_set_orth,:));
    orth_x(i,:) = cur_orth_x;
    orth_y(i,:) = cur_orth_y;
    
    par_prof_norm = parallel_profile(i,:) - min(parallel_profile(i,:));
    par_prof_norm = par_prof_norm/nansum(par_prof_norm);
    orth_prof_norm = orth_profile(i,:) - min(orth_profile(i,:));
    orth_prof_norm = orth_prof_norm/nansum(orth_prof_norm);
        
    init_mean_par_s(i) = nansum(cur_par_so.*par_prof_norm);
    init_mean_orth_s(i) = nansum(cur_orth_so.*orth_prof_norm);
    init_std_par(i) = sqrt(nansum((cur_par_so - init_mean_par_s(i)).^2.*par_prof_norm));
    init_std_orth(i) = sqrt(nansum((cur_orth_so - init_mean_orth_s(i)).^2.*orth_prof_norm));
    
    beta0(1) = min(parallel_profile(i,:)); beta0(2) = range(parallel_profile(i,:));
    beta0(3) = init_mean_par_s(i); beta0(4) = init_std_par(i)^2;
    [gfit_par(i,:),~,J,~,par_mse(i)] = nlinfit(cur_par_so,parallel_profile(i,:),@(X,cur_para_so) my_gaussfun(X,cur_par_so),beta0);
    par_Jcond(i) = cond(J);
    gpred_par(i,:) = my_gaussfun(gfit_par(i,:),cur_par_so);
    dom_par(i) = (max(gpred_par(i,:)) - min(gpred_par(i,:)))/ov_avg_rates(i);
    
    beta0(1) = min(orth_profile(i,:)); beta0(2) = range(orth_profile(i,:));
    beta0(3) = init_mean_orth_s(i); beta0(4) = init_std_orth(i)^2;
    [gfit_orth(i,:),~,J,~,orth_mse(i)] = nlinfit(cur_orth_so,orth_profile(i,:),@(X,cur_orth_so) my_gaussfun(X,cur_orth_so),beta0);
    orth_Jcond(i) = cond(J);
    gpred_orth(i,:) = my_gaussfun(gfit_orth(i,:),cur_orth_so);
    dom_orth(i) = (max(gpred_orth(i,:)) - min(gpred_orth(i,:)))/ov_avg_rates(i);
    
    [~,parloc] = min(abs(cur_par_so-gfit_par(i,3)));
    [~,orthloc] = min(abs(cur_orth_so-gfit_orth(i,3)));    
    mean_x_old(i) = 0.5*cur_par_x(parloc) + 0.5*cur_orth_x(orthloc);
    mean_y_old(i) = 0.5*cur_par_y(parloc) + 0.5*cur_orth_y(orthloc);
     mean_x(i) = cur_par_so(parloc)*cos(degtorad(vec_angles(cur_set_par))) + cur_orth_so(orthloc)*cos(degtorad(vec_angles(cur_set_orth)));
     mean_y(i) = cur_par_so(parloc)*sin(degtorad(vec_angles(cur_set_par))) + cur_orth_so(orthloc)*sin(degtorad(vec_angles(cur_set_orth)));
   
%     beta0(1) = 0; beta0(2) = max(orth_profile(i,:));
%     beta0(3) = mean_orth_s(i); beta0(4) = std_orth(i)^2;
%     gfit_orth(i,:) = nlinfit(cur_orth_so,orth_profile(i,:),@(X,cur_orth_so) my_gaussfun(X,cur_orth_so),beta0);
end

mean_x = mean_x + xcent;
mean_y = mean_y + ycent;

ov_peak_rates = max(max(parallel_profile,[],2),max(orth_profile,[],2));
norm_parallel_profile = bsxfun(@rdivide,parallel_profile,ov_peak_rates);
norm_orth_profile = bsxfun(@rdivide,orth_profile,ov_peak_rates);

std_par = sqrt(gfit_par(:,4));
std_orth = sqrt(gfit_orth(:,4));

good_fits = find(par_Jcond < 1000 & orth_Jcond < 1000);

save oned_fixation_fits_v3 std_* *_profile mean_* gpred_* gfit_* dom_* parallel_* orth_* good_fits

%%
for i = 1:96
 
    cur_set_par = find(avg_bar_ars < 1 & avg_bar_oris == pref_bar_oris(i));
    cur_par_x = ov_x_profile(cur_set_par,:);
    cur_par_y = ov_y_profile(cur_set_par,:);
    cur_par_so = ov_s_profile(cur_set_par,:);
    
    cur_set_orth = find(avg_bar_ars > 1 & avg_bar_oris == pref_bar_oris(i));
    cur_orth_x = ov_x_profile(cur_set_orth,:);
    cur_orth_y = ov_y_profile(cur_set_orth,:);
    cur_orth_so = ov_s_profile(cur_set_orth,:);

    subplot(2,2,1)
    plot(cur_par_x,parallel_profile(i,:))
    subplot(2,2,2)
    plot(cur_par_y,parallel_profile(i,:))
    subplot(2,2,3)
    plot(cur_orth_x,orth_profile(i,:))
    subplot(2,2,4)
    plot(cur_orth_y,orth_profile(i,:))
    
    fprintf('Mean X: %.3f  Mean Y: %.3f\n',mean_x(i),mean_y(i));
    
    
    pause
    clf
end

    
    
    
