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
    
    [b,a] = butter(2,0.01/(eye_fsd/2),'low'); %remove low-frequency changes in eye signal
    
    lEye_d = bsxfun(@plus,eye_vals_d(:,1:2)*left_gain,left_offset);
    rEye_d = bsxfun(@plus,eye_vals_d(:,3:4)*right_gain,right_offset);
    lEye_di = interp1(find(used_points_d==1),lEye_d(used_points_d==1,:),1:length(used_points_d));
    lEye_di(isnan(lEye_di(:,1)),1) = 0;lEye_di(isnan(lEye_di(:,2)),2) = 0;
    rEye_di = interp1(find(used_points_d==1),rEye_d(used_points_d==1,:),1:length(used_points_d));
    rEye_di(isnan(rEye_di(:,1)),1) = 0;rEye_di(isnan(rEye_di(:,2)),2) = 0;
    
    low_Fh = filtfilt(b,a,lEye_di(:,1));
    low_Fv = filtfilt(b,a,lEye_di(:,2));
    lEye_di(:,1) = lEye_di(:,1) - low_Fh;
    lEye_di(:,2) = lEye_di(:,2) - low_Fv;
    low_Fh = filtfilt(b,a,rEye_di(:,1));
    low_Fv = filtfilt(b,a,rEye_di(:,2));
    rEye_di(:,1) = rEye_di(:,1) - low_Fh;
    rEye_di(:,2) = rEye_di(:,2) - low_Fv;
    avg_di = 0.5*lEye_di + 0.5*rEye_di;
    
    %     low_Fh = filtfilt(b,a,avg_di(:,1));
    %     low_Fv = filtfilt(b,a,avg_di(:,2));
    %     avg_di(:,1) = avg_di(:,1) - low_Fh;
    %     avg_di(:,2) = avg_di(:,2) - low_Fv;
    
    win_size = 0.5;
    out_of_bounds = find(abs(avg_di(:,1)) > win_size | abs(avg_di(:,2)) > win_size);
    used_points_d(ismember(used_points_d,out_of_bounds)) = 0;
    
    lEye_d(used_points_d==0,:) = nan;
    rEye_d(used_points_d==0,:) = nan;
    avg_di(used_points_d==0,:) = nan;
    %
    
    %     lEye_d = eye_vals_d(:,1:2); lEye_d(:,1) = -lEye_d(:,1);
    %     rEye_d = eye_vals_d(:,3:4); rEye_d(:,1) = -rEye_d(:,1);
    %      lEye_d(used_points_d==0,:) = nan;
    %     rEye_d(used_points_d==0,:) = nan;
    
    %     lEye_d = bsxfun(@minus,lEye_d,nanmedian(lEye_d));
    %     rEye_d = bsxfun(@minus,rEye_d,nanmedian(rEye_d));
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
    
    %% ADD IN Eye Tracking CORRECTION
    unique_oris = unique(line_oris);
    unique_ars = unique(ar);
    
    %     bar_xo_vec_cor = bar_xo_vec - rEye_d(:,1);
    %     bar_yo_vec_cor = bar_yo_vec - rEye_d(:,2);
    bar_xo_vec_cor = bar_xo_vec - avg_di(:,1);
    bar_yo_vec_cor = bar_yo_vec - avg_di(:,2);
    x_correction(ee) = -nanmean(avg_di(used_points_d==1,1));
    y_correction(ee) = -nanmean(avg_di(used_points_d==1,2));
    
    bar_so_vec_cor = bar_so_vec;
    for i = 1:length(unique_oris)
        cur_set = find(bar_ori_vec == unique_oris(i) & bar_ar_vec == unique_ars(1));
        cur_angle = degtorad(unique_oris(i) + 180);
        %         cur_angle = degtorad(unique_oris(i));
        bar_so_vec_cor(cur_set) = bar_so_vec_cor(cur_set) - avg_di(cur_set,1)*cos(cur_angle);
        bar_so_vec_cor(cur_set) = bar_so_vec_cor(cur_set) - avg_di(cur_set,2)*sin(cur_angle);
        
        cur_set = find(bar_ori_vec == unique_oris(i) & bar_ar_vec == unique_ars(2));
        cur_angle = degtorad(unique_oris(i) + 90);
        %         cur_angle = degtorad(unique_oris(i) + 90);
        bar_so_vec_cor(cur_set) = bar_so_vec_cor(cur_set) - avg_di(cur_set,1)*cos(cur_angle);
        bar_so_vec_cor(cur_set) = bar_so_vec_cor(cur_set) - avg_di(cur_set,2)*sin(cur_angle);
    end
    
    %     bar_so_vec(isnan(lEye_d(:,1))) = nan;
    %%  
    min_samps = 5;
    delay = round(0.05*eye_fsd);
    para_avg_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    orth_avg_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    para_cavg_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    orth_cavg_profile(ee,:,:,:) = nan(96,length(unique_oris),n_pos_bins);
    for i = 1:length(unique_oris)
        cur_set = find(bar_ori_vec == unique_oris(i) & bar_ar_vec == unique_ars(1));
        [bins,indx_cor] = histc(bar_so_vec_cor(cur_set),s_range_edges);
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
            
            cur_inds = cur_set(indx_cor==j);
            cur_inds = cur_inds + delay;
            cur_inds(cur_inds > length(bar_so_vec)) = [];
            if length(cur_inds) > min_samps
                para_cavg_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
            end
            avg_cpara_xbin(i,j) = nanmean(bar_xo_vec_cor(cur_inds));
            avg_cpara_ybin(i,j) = nanmean(bar_yo_vec_cor(cur_inds));
            avg_cpara_sbin(i,j) = nanmean(bar_so_vec_cor(cur_inds));
            n_para_cbins(i,j) = length(cur_inds);
        end
        %         [bins,indx_cor] = histc(bar_yo_vec_cor(cur_set),y_range_edges);
        %         [bins,indx] = histc(bar_yo_vec(cur_set),y_range_edges);
        %
        %         n_ov_para(i) = length(cur_set);
        %         for j = 1:n_pos_bins
        %             cur_inds = cur_set(indx==j);
        %             cur_inds = cur_inds + delay;
        %             cur_inds(cur_inds > length(bar_so_vec)) = [];
        %             if length(cur_inds) > 20
        %                 para_avg_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
        %             end
        %             n_para_bins(i,j) = length(cur_inds);
        %
        %             cur_inds = cur_set(indx_cor==j);
        %             cur_inds = cur_inds + delay;
        %             cur_inds(cur_inds > length(bar_so_vec)) = [];
        %             if length(cur_inds) > 20
        %                 para_cavg_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
        %             end
        %             n_para_cbins(i,j) = length(cur_inds);
        %         end
        
        cur_set = find(bar_ori_vec == unique_oris(i) & bar_ar_vec == unique_ars(2));
        n_ov_orth(i) = length(cur_set);
        [bins,indx_cor] = histc(bar_so_vec_cor(cur_set),s_range_edges);
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
            
            cur_inds = cur_set(indx_cor==j);
            cur_inds = cur_inds + delay;
            cur_inds(cur_inds > length(bar_so_vec)) = [];
            if length(cur_inds) > min_samps
                orth_cavg_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
            end
            n_orth_cbins(i,j) = length(cur_inds);
            avg_corth_xbin(i,j) = nanmean(bar_xo_vec_cor(cur_inds));
            avg_corth_ybin(i,j) = nanmean(bar_yo_vec_cor(cur_inds));
            avg_corth_sbin(i,j) = nanmean(bar_so_vec_cor(cur_inds));
        end
        %         cur_set = find(bar_ori_vec == unique_oris(i) & bar_ar_vec == unique_ars(2));
        %         n_ov_orth(i) = length(cur_set);
        %         [bins,indx_cor] = histc(bar_yo_vec_cor(cur_set),y_range_edges);
        %         [bins,indx] = histc(bar_yo_vec(cur_set),y_range_edges);
        %         for j = 1:n_pos_bins
        %             cur_inds = cur_set(indx==j);
        %             cur_inds = cur_inds + delay;
        %             cur_inds(cur_inds > length(bar_so_vec)) = [];
        %             if length(cur_inds) > 20
        %                 orth_avg_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
        %             end
        %             n_orth_bins(i,j) = length(cur_inds);
        %
        %             cur_inds = cur_set(indx_cor==j);
        %             cur_inds = cur_inds + delay;
        %             cur_inds(cur_inds > length(bar_so_vec)) = [];
        %             if length(cur_inds) > 20
        %                 orth_cavg_profile(ee,:,i,j) = mean(binned_spks(:,cur_inds),2);
        %             end
        %             n_orth_cbins(i,j) = length(cur_inds);
        %         end
        
    end
    
end

%%
cd ~/Data/bruce/7_15_12/G029
save oned_noise_mapping *profile unique* *bin* 
%%
ov_para_profile = squeeze(nanmean(para_avg_profile));
ov_orth_profile = squeeze(nanmean(orth_avg_profile));
ov_para_cprofile = squeeze(nanmean(para_cavg_profile));
ov_orth_cprofile = squeeze(nanmean(orth_cavg_profile));
ov_avg_rates = mean(avg_rates);
%%
% for c = single_units{1};
% for c = 1:96;
%     close all
%     figure
%     subplot(2,2,1)
%     plot(squeeze(ov_para_profile(c,:,:))')
%     yl = ylim();
%     subplot(2,2,3)
%     plot(squeeze(ov_para_cprofile(c,:,:))')
%     ylim(yl)
%     subplot(2,2,2)
%     plot(squeeze(ov_orth_profile(c,:,:))')
%     yl = ylim();
%     subplot(2,2,4)
%     plot(squeeze(ov_orth_cprofile(c,:,:))')
%     ylim(yl)
%     pause
% end
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
for i = 1:96
    
    cur_set = find(avg_bar_ars < 1 & avg_bar_oris == pref_bar_oris(i));
    cur_par_x = ov_x_profile(cur_set,:);
    cur_par_y = ov_y_profile(cur_set,:);
    cur_par_so = ov_s_profile(cur_set,:);
    parallel_profile(i,:) =  squeeze(ov_avg_profile(i,cur_set,:));
    parallel_x(i,:) = cur_par_x;
    parallel_y(i,:) = cur_par_y;
    
    cur_set = find(avg_bar_ars > 1 & avg_bar_oris == pref_bar_oris(i));
    cur_orth_x = ov_x_profile(cur_set,:);
    cur_orth_y = ov_y_profile(cur_set,:);
    cur_orth_so = ov_s_profile(cur_set,:);
    orth_profile(i,:) = squeeze(ov_avg_profile(i,cur_set,:));
    orth_x(i,:) = cur_orth_x;
    orth_y(i,:) = cur_orth_y;
    
    par_prof_norm = parallel_profile(i,:) - min(parallel_profile(i,:));
    par_prof_norm = par_prof_norm/nansum(par_prof_norm);
    orth_prof_norm = orth_profile(i,:) - min(orth_profile(i,:));
    orth_prof_norm = orth_prof_norm/nansum(orth_prof_norm);
    
    init_mean_x(i) = 0.5*nansum(cur_par_x.*par_prof_norm) + 0.5*nansum(cur_orth_x.*orth_prof_norm);
    init_mean_y(i) = 0.5*nansum(cur_par_y.*par_prof_norm) + 0.5*nansum(cur_orth_y.*orth_prof_norm);
    
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
    mean_x(i) = 0.5*cur_par_x(parloc) + 0.5*cur_orth_x(orthloc);
    mean_y(i) = 0.5*cur_par_y(parloc) + 0.5*cur_orth_y(orthloc);
    
%     beta0(1) = 0; beta0(2) = max(orth_profile(i,:));
%     beta0(3) = mean_orth_s(i); beta0(4) = std_orth(i)^2;
%     gfit_orth(i,:) = nlinfit(cur_orth_so,orth_profile(i,:),@(X,cur_orth_so) my_gaussfun(X,cur_orth_so),beta0);
end

ov_peak_rates = max(max(parallel_profile,[],2),max(orth_profile,[],2));
norm_parallel_profile = bsxfun(@rdivide,parallel_profile,ov_peak_rates);
norm_orth_profile = bsxfun(@rdivide,orth_profile,ov_peak_rates);

std_par = sqrt(gfit_par(:,4));
std_orth = sqrt(gfit_orth(:,4));

save oned_fixation_fits std_* *_profile mean_* gpred_* gfit_* dom_* parallel_* orth_*
%%
% ov_avg_cprofile = cat(2,ov_para_cprofile,ov_orth_cprofile);
% ov_x_cprofile = [avg_cpara_xbin; avg_corth_xbin];
% ov_y_cprofile = [avg_cpara_ybin; avg_corth_ybin];
% ov_s_cprofile = [avg_cpara_sbin; avg_corth_sbin];
% max_cprofiles = max(ov_avg_cprofile,[],3);
% min_cprofiles = min(ov_avg_cprofile,[],3);
% avg_bar_oris = [unique_oris(:); unique_oris(:)];
% avg_bar_ars = [zeros(4,1); 10*ones(4,1)];
% for i = 1:length(unique_oris)
%     used_set = find(avg_bar_oris == unique_oris(i));
%     mod_cindex(i,:) = (max(max_cprofiles(:,used_set),[],2) - min(min_cprofiles(:,used_set),[],2));
% end
% [bar_ori_maxrates,cpref_bar_oris] = max(mod_cindex);
% cpref_bar_oris = unique_oris(cpref_bar_oris);
% for i = 1:96
%     cur_set = find(avg_bar_ars < 1 & avg_bar_oris == cpref_bar_oris(i));
%     cur_par_x = ov_x_cprofile(cur_set,:);
%     cur_par_y = ov_y_cprofile(cur_set,:);
%     cur_par_so = ov_s_cprofile(cur_set,:);
%     parallel_cprofile(i,:) =  squeeze(ov_avg_cprofile(i,cur_set,:));
% %     cur_offset_par = min(parallel_cprofile(i,:));
% %     parallel_cprofile(i,:) = parallel_cprofile(i,:) - cur_offset_par;
%     
%     cur_set = find(avg_bar_ars > 1 & avg_bar_oris == cpref_bar_oris(i));
%     cur_orth_x = ov_x_cprofile(cur_set,:);
%     cur_orth_y = ov_y_cprofile(cur_set,:);
%     cur_orth_so = ov_s_cprofile(cur_set,:);
%     orth_cprofile(i,:) = squeeze(ov_avg_cprofile(i,cur_set,:));
% %     cur_offset_orth = min(orth_cprofile(i,:));
% %     orth_cprofile(i,:) = orth_cprofile(i,:) - cur_offset_orth;
%     
%     par_prof_norm = parallel_cprofile(i,:) - min(parallel_cprofile(i,:));
%     par_prof_norm = par_prof_norm/nansum(par_prof_norm);
%     orth_prof_norm = orth_cprofile(i,:) - min(orth_cprofile(i,:));
%     orth_prof_norm = orth_prof_norm/nansum(orth_prof_norm);
%     
%     init_cmean_x(i) = 0.5*nansum(cur_par_x.*par_prof_norm) + 0.5*nansum(cur_orth_x.*orth_prof_norm);
%     init_cmean_y(i) = 0.5*nansum(cur_par_y.*par_prof_norm) + 0.5*nansum(cur_orth_y.*orth_prof_norm);
%     
%     init_cmean_par_s(i) = nansum(cur_par_so.*par_prof_norm);
%     init_cmean_orth_s(i) = nansum(cur_orth_so.*orth_prof_norm);
%     
%     init_cstd_par(i) = sqrt(nansum((cur_par_so - init_cmean_par_s(i)).^2.*par_prof_norm));
%     init_cstd_orth(i) = sqrt(nansum((cur_orth_so - init_cmean_orth_s(i)).^2.*orth_prof_norm));
%     
%     beta0(1) = 0; beta0(2) = max(parallel_cprofile(i,:));
%     beta0(3) = init_cmean_par_s(i); beta0(4) = init_cstd_par(i)^2;
%     gfit_cpar(i,:) = nlinfit(cur_par_so,parallel_cprofile(i,:),@(X,cur_para_so) my_gaussfun(X,cur_par_so),beta0);
%     gpred_cpar(i,:) = my_gaussfun(gfit_cpar(i,:),cur_par_so);
%     dom_cpar(i) = (max(gpred_cpar(i,:)) - min(gpred_cpar(i,:)))/ov_avg_rates(i);
%     
%     options.Display = 'Off';
%     beta0(1) = 0; beta0(2) = max(orth_cprofile(i,:));
%     beta0(3) = init_cmean_orth_s(i); beta0(4) = init_cstd_orth(i)^2;
%     gfit_corth(i,:) = nlinfit(cur_orth_so,orth_cprofile(i,:),@(X,cur_orth_so) my_gaussfun(X,cur_orth_so),beta0,options);
%     gpred_corth(i,:) = my_gaussfun(gfit_corth(i,:),cur_orth_so);
%     dom_corth(i) = (max(gpred_corth(i,:)) - min(gpred_corth(i,:)))/ov_avg_rates(i);
% 
%     [~,parloc] = min(abs(cur_par_so-gfit_cpar(i,3)));
%     [~,orthloc] = min(abs(cur_orth_so-gfit_corth(i,3)));    
%     cmean_x(i) = 0.5*cur_par_x(parloc) + 0.5*cur_orth_x(orthloc);
%     cmean_y(i) = 0.5*cur_par_y(parloc) + 0.5*cur_orth_y(orthloc);
% 
% end
% 
% ov_cpeak_rates = max(max(parallel_cprofile,[],2),max(orth_cprofile,[],2));
% norm_cparallel_profile = bsxfun(@rdivide,parallel_cprofile,ov_cpeak_rates);
% norm_corth_profile = bsxfun(@rdivide,orth_cprofile,ov_cpeak_rates);
% 
%%
std_par = sqrt(gfit_par(:,4));
std_orth = sqrt(gfit_orth(:,4));
cstd_par = sqrt(gfit_cpar(:,4));
cstd_orth = sqrt(gfit_corth(:,4));

% dom_par = gfit_par(:,2);
% dom_orth = gfit_orth(:,2);
% cdom_par = gfit_cpar(:,2);
% cdom_orth = gfit_corth(:,2);

%%
ov_used_set = single_units{1};
max_std = 0.5;
used_set = ov_used_set(std_par(ov_used_set) < max_std & std_orth(ov_used_set) < max_std);
% used_set = 1:96;
figure; hold on
cmap = jet(length(used_set));
for i = 1:length(used_set)
    ellipse(std_orth(used_set(i)),std_par(used_set(i)),degtorad(pref_bar_oris(used_set(i))),...
        mean_x(used_set(i)),mean_y(used_set(i)),cmap(i,:));
    plot(mean_x(used_set(i)),mean_y(used_set(i)),'o','color',cmap(i,:));
end
xlim([0 0.7]);ylim([-0.8 -0.1])

used_set = ov_used_set(cstd_par(ov_used_set) < max_std & cstd_orth(ov_used_set) < max_std);
figure; hold on
cmap = jet(length(used_set));
for i = 1:length(used_set)
    
    ellipse(cstd_orth(used_set(i)),cstd_par(used_set(i)),degtorad(cpref_bar_oris(used_set(i))),...
        cmean_x(used_set(i)),cmean_y(used_set(i)),cmap(i,:));
    plot(cmean_x(used_set(i)),cmean_y(used_set(i)),'o','color',cmap(i,:));
end
xlim([0 0.7]);ylim([-0.8 -0.1])

%%
figure
subplot(2,1,1)
plot(std_par,cstd_par,'o')
hold on
plot(std_par(single_units{1}),cstd_par(single_units{1}),'r.')
line([0.0 0.35],[0. 0.35],'color','k')
xlim([0. 0.35]);ylim([0 0.35])
subplot(2,1,2)
plot(std_orth,cstd_orth,'o')
hold on
plot(std_orth(single_units{1}),cstd_orth(single_units{1}),'r.')
line([0.0 0.35],[0. 0.35],'color','k')
xlim([0. 0.35]);ylim([0 0.35])

figure
subplot(2,1,1)
plot(dom_par,dom_cpar,'o')
hold on
plot(dom_par(single_units{1}),dom_cpar(single_units{1}),'r.')
line([0.0 15],[0. 15],'color','k')
xlim([0. 15]);ylim([0 15])
subplot(2,1,2)
plot(dom_orth,dom_corth,'o')
hold on
plot(dom_orth(single_units{1}),dom_corth(single_units{1}),'r.')
line([0.0 15],[0. 15],'color','k')
xlim([0. 15]);ylim([0 15])
fprintf('Parallel std diff: %.3f\n',signrank(std_par,cstd_par));
fprintf('Orth std diff: %.3f\n',signrank(std_orth,cstd_orth));
fprintf('Parallel unit std diff: %.3f\n',signrank(std_par(single_units{1}),cstd_par(single_units{1})));
fprintf('Orth unit std diff: %.3f\n',signrank(std_orth(single_units{1}),cstd_orth(single_units{1})));

fprintf('Parallel dom diff: %.3f\n',signrank(dom_par,dom_cpar));
fprintf('Orth dom diff: %.3f\n',signrank(dom_orth,dom_corth));
fprintf('Parallel unit dom diff: %.3f\n',signrank(dom_par(single_units{1}),dom_cpar(single_units{1})));
fprintf('Orth unit dom diff: %.3f\n',signrank(dom_orth(single_units{1}),dom_corth(single_units{1})));

%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

%%
bad_fits = find(par_Jcond > 1e3 | orth_Jcond > 1e3);
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt);hold on  
            if ~ismember(cur,bad_fits)
            ellipse(std_orth(cur),std_par(cur),degtorad(pref_bar_oris(cur)),...
                mean_x(cur),mean_y(cur));
            plot(mean_x(cur),mean_y(cur),'.','markersize',20);
            xlim([0.1 0.7]);ylim([-0.7 -0.1])
            line([0.1 0.7],[-0.4 -0.4],'color','k')
            line([0.4 0.4],[-0.7 -0.1],'color','k')
            end
        end
        cnt = cnt + 1;
    end
end
%%
% for c = 1:96;
%     
%     cur_set = find(avg_bar_ars < 1 & avg_bar_oris == cpref_bar_oris(i));
%     cur_par_x = ov_x_cprofile(cur_set,:);
%     cur_par_y = ov_y_cprofile(cur_set,:);
%     cur_par_so = ov_s_cprofile(cur_set,:);
% %     gaus_fun = normpdf(cur_par_so,cmean_par_s(c),cstd_par(c));
% %     gaus_fun = gaus_fun/max(gaus_fun)*max(parallel_cprofile(c,:));
%     temp = gfit_cpar(c,:); 
%     pred = my_gaussfun(temp,cur_par_so);
%     
%     figure
%     set(gcf,'Position',[500 600 800 800]);subplot(2,2,1)
%     plot(cur_par_so,parallel_cprofile(c,:))
%     hold on
% %     plot(cur_par_so,gaus_fun,'r')
%     plot(cur_par_so,pred,'k')
%     
%     cur_set = find(avg_bar_ars > 1 & avg_bar_oris == cpref_bar_oris(i));
%     cur_orth_x = ov_x_cprofile(cur_set,:);
%     cur_orth_y = ov_y_cprofile(cur_set,:);
%     cur_orth_so = ov_s_cprofile(cur_set,:);
% %     gaus_fun = normpdf(cur_par_so,cmean_orth_s(c),cstd_orth(c));
% %     gaus_fun = gaus_fun/max(gaus_fun)*max(orth_cprofile(c,:));
%     temp = gfit_corth(c,:); 
%     pred = my_gaussfun(temp,cur_par_so);
%     
%     subplot(2,2,2)
%     plot(cur_par_so,orth_cprofile(c,:))
%     hold on
% %     plot(cur_par_so,gaus_fun,'r')
%     plot(cur_par_so,pred,'k')
%     
%     cur_set = find(avg_bar_ars < 1 & avg_bar_oris == pref_bar_oris(i));
%     cur_para_x = ov_x_profile(cur_set,:);
%     cur_para_y = ov_y_profile(cur_set,:);
%     cur_para_so = ov_s_profile(cur_set,:);
% %     gaus_fun = normpdf(cur_par_so,mean_par_s(c),std_par(c));
% %     gaus_fun = gaus_fun/max(gaus_fun)*max(parallel_profile(c,:));
%     temp = gfit_par(c,:); 
%     pred = my_gaussfun(temp,cur_par_so);
%     
%     subplot(2,2,3)
%     plot(cur_par_so,parallel_profile(c,:))
%     hold on
% %     plot(cur_par_so,gaus_fun,'r')
%     plot(cur_par_so,pred,'k')
% 
%      cur_set = find(avg_bar_ars > 1 & avg_bar_oris == pref_bar_oris(i));
%     cur_orth_x = ov_x_profile(cur_set,:);
%     cur_orth_y = ov_y_profile(cur_set,:);
%     cur_orth_so = ov_s_profile(cur_set,:);
% %     gaus_fun = normpdf(cur_par_so,mean_orth_s(c),std_orth(c));
% %     gaus_fun = gaus_fun/max(gaus_fun)*max(orth_profile(c,:));
%     temp = gfit_orth(c,:); 
%     pred = my_gaussfun(temp,cur_par_so);
%     
%     subplot(2,2,4)
%     plot(cur_par_so,orth_profile(c,:))
%     hold on
% %     plot(cur_par_so,gaus_fun,'r')
%     plot(cur_par_so,pred,'k')
% 
%    pause
%     close all
% end