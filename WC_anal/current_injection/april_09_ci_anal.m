clear all

% depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_1s_POS_CI';
depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_1s_POS_CI';
depol_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_1s_POS';
% depol_array{4} = 'C:\WC_Germany\april_09_data\2009-04-21_CWC_LFP\2009-4-21_12s_1s_POS_CI';
% depol_array{5} = 'C:\WC_Germany\april_09_data\2009-05-16_CWC_LFP\2009-5-16_1s_POS_CI';

% hypol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_1s_NEG_CI';
hypol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_1s_NEG_CI';
hypol_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_1s_NEG';
% hypol_array{4} = 'C:\WC_Germany\april_09_data\2009-04-21_CWC_LFP\2009-4-21_12s_1s_NEG_CI';
% hypol_array{5} = 'C:\WC_Germany\april_09_data\2009-05-16_CWC_LFP\2009-5-16_1s_NEG_CI';

% prior_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_spontaneous';
prior_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_spontaneous';
prior_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_spontaneous';
% prior_array{4} = 'C:\WC_Germany\april_09_data\2009-04-21_CWC_LFP\2009-4-21_12s_spontaneous';
% prior_array{5} = 'C:\WC_Germany\april_09_data\2009-05-16_CWC_LFP\2009-5-16_spontaneous';

% save_name{1} = '2009-4-5';
save_name{1} = '2009-4-7';
save_name{2} = '2009-04-13';
% save_name{4} = '2009-04-21';
% save_name{5} = '2009-05-16';

ci_time = 1;
wait_time = 2;
first_pulse = 1;

Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;

mp_min = -100;
mp_max = 10;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

load new_1s_dep_data_minus_spike

for d = 1:2

    %% get prior distribution
%     load(prior_array{d})
% 
%     prior_dist(d,:) = gpkde(data,2,[mp_min;mp_max;nbins]);
%     prior_dist_beg(d,:) = gpkde(data(1:6e6),2,[mp_min;mp_max;nbins]);
%     prior_dist_end(d,:) = gpkde(data(end-6e6:end),2,[mp_min;mp_max;nbins]);
%     
%     mean_prior(d) = mean(data);
%     std_prior(d) = std(data);

    %% get depolarizing CI
    load(depol_array{d})
    
    data = data_minus_spike{d};
    
    time = (1:length(data))/Fs;
    data = downsample(data,dsf);
    time = downsample(time,dsf);

    num_cis(d) = floor((max(time)-2*first_pulse)/(ci_time+wait_time))

    sweep_time = (ci_time+wait_time)*Fsd;

    depol_mat = zeros(num_cis(d),sweep_time);

    for i = 1:num_cis(d)

        begpt = (i-1)*3*Fsd+1;
        endpt = begpt+3*Fsd;
        depol_mat(i,:) = data(begpt:endpt-1);

    end

    tvec = (1:sweep_time)/Fsd;

%     tdep_dist = zeros(length(tvec),length(mp_range));
% 
%     for t = 1:length(tvec)
%         tdep_dist(t,:) = gpkde(depol_mat(:,t),2,[mp_min;mp_max;nbins]);
%     end

%     norm_depol_mat = (depol_mat - mean_prior(d))/std_prior(d);

    dep_end = find(tvec > 2.0,1,'first');
    dep_start = find(tvec > 1.0,1,'first');
%     norm_mean_afterdep_traj(d,:) = mean(norm_depol_mat);
    mean_afterdep_traj(d,:) = mean(depol_mat);
%     mean_afterdep(d,:) = mean(tdep_dist(dep_end+100:dep_end+round(0.5*Fsd),:));
%     mean_beforedep(d,:) = mean(tdep_dist(round(0.5*Fsd):dep_start-100,:));
%     mean_duringdep(d,:) = mean(tdep_dist(dep_start:dep_end,:));
    
    
%         load(hypol_array{d})
% 
%     time = (1:length(data))/Fs;
%     data = downsample(data,dsf);
%     time = downsample(time,dsf);
% 
%     num_cis = floor((max(time)-2*first_pulse)/(ci_time+wait_time))
% 
%     sweep_time = (ci_time+wait_time)*Fsd;
% 
%     hypol_mat = zeros(num_cis,sweep_time);
% 
%     for i = 1:num_cis(d)
% 
%         begpt = (i-1)*3*Fsd+1;
%         endpt = begpt+3*Fsd;
%         hypol_mat(i,:) = data(begpt:endpt-1);
% 
%     end
% 
%     tvec = (1:sweep_time)/Fsd;
% 
% %     thyp_dist = zeros(length(tvec),length(mp_range));
% % 
% %     for t = 1:length(tvec)
% %         thyp_dist(t,:) = gpkde(hypol_mat(:,t),2,[mp_min;mp_max;nbins]);
% %     end
% 
% %     norm_hypol_mat = (hypol_mat - mean_prior(d))/std_prior(d);
% 
%     hyp_end = find(tvec > 2.0,1,'first');
%     hyp_start = find(tvec > 1.0,1,'first');
%     mean_afterhyp_mp(d,:) = mean(mean(hypol_mat(:,hyp_end+100:hyp_end+round(0.5*Fsd))));
% %     norm_mean_afterhyp_traj(d,:) = mean(norm_hypol_mat);
%     mean_afterhyp_traj(d,:) = mean(hypol_mat);
% %     mean_afterhyp(d,:) = mean(thyp_dist(hyp_end+100:hyp_end+round(0.5*Fsd),:));
% %     mean_beforehyp(d,:) = mean(thyp_dist(round(0.5*Fsd):hyp_start-100,:));
% %     mean_duringhyp(d,:) = mean(thyp_dist(hyp_start:hyp_end,:));
% %   
%     
%     avg_amp_before_dep(d) = mean(mean(depol_mat(:,round(0.5*Fsd):dep_start-100)));
%     avg_amp_during_dep(d) = mean(mean(depol_mat(:,dep_start:dep_end)));
%     avg_amp_before_hyp(d) = mean(mean(hypol_mat(:,round(0.5*Fsd):hyp_start-100)));
%     avg_amp_during_hyp(d) = mean(mean(hypol_mat(:,hyp_start:hyp_end)));
%     
    
%     plot(mp_range,prior_dist(d,:),'linewidth',2)
%     hold on
%     plot(mp_range,prior_dist_beg(d,:))
%     plot(mp_range,prior_dist_end(d,:),'--')
%     plot(mp_range,mean_afterdep(d,:),'r')
%     plot(mp_range,mean_beforedep(d,:),'r--')
%     plot(mp_range,mean_afterhyp(d,:),'g')
%     plot(mp_range,mean_beforehyp(d,:),'g--')
%     xlim([-90 -20])
%     t_name = ['C:\WC_Germany\current_injection\persistent_ci\ci_mp_dist\' save_name{d}];
%     print('-dpng',t_name);close
    
%     plot(tvec,mean_afterdep_traj(d,:),'r')
%     hold on
%     plot(tvec,mean_afterhyp_traj(d,:),'g')
%     t_name = ['C:\WC_Germany\current_injection\persistent_ci\ci_mp_traj\' save_name{d}];
%     print('-dpng',t_name);close

    
    
    
end


cd C:\WC_Germany\current_injection\persistent_ci
% save april_09_1s_ci_data tvec mean_* avg_amp* norm_mean* prior*
% save april_09_1s_ci_data tvec mean_* 

%% comparison
% for d = 1:2
% 
%     [ad_peak_amps,ad_peak_locs] = findpeaks(mean_afterhyp(d,:),'minpeakheight',0.01);
%     [pr_peak_amps,pr_peak_locs] = findpeaks(prior_dist(d,:),'minpeakheight',0.01);
% 
%     if length(ad_peak_locs) == 2
%         n_ad_down_mp(d) = mp_range(ad_peak_locs(1));
%         n_ad_up_mp(d) = mp_range(ad_peak_locs(2));
%         [dummy,minloc] = min(mean_afterhyp(d,ad_peak_locs(1):ad_peak_locs(2)));
%         n_ad_midpt(d) = minloc+ad_peak_locs(1);
%     elseif length(ad_peak_locs) > 2
%         [dummy, peakorder] = sort(ad_peak_amps,'descend');
%         ad_peak_locs(peakorder(3:end)) = [];
%         n_ad_down_mp(d) = mp_range(ad_peak_locs(1));
%         n_ad_up_mp(d) = mp_range(ad_peak_locs(2));
%         [dummy,minloc] = min(mean_afterhyp(d,ad_peak_locs(1):ad_peak_locs(2)));
%         n_ad_midpt(d) = minloc+ad_peak_locs(1);
%     else
%         n_ad_down_mp(d) = mp_range(ad_peak_locs(1));
%         n_ad_up_mp(d) = nan;
%         n_ad_midpt(d) = length(mp_range);
%     end
% 
%     if length(pr_peak_locs) == 2
%         n_pr_down_mp(d) = mp_range(pr_peak_locs(1));
%         n_pr_up_mp(d) = mp_range(pr_peak_locs(2));
%         [dummy,minloc] = min(prior_dist(d,pr_peak_locs(1):pr_peak_locs(2)));
%         n_pr_midpt(d) = minloc+pr_peak_locs(1);
%     elseif length(pr_peak_locs) > 2
%         [dummy, peakorder] = sort(pr_peak_amps,'descend');
%         pr_peak_locs(peakorder(3:end)) = [];
%         n_pr_down_mp(d) = mp_range(pr_peak_locs(1));
%         n_pr_up_mp(d) = mp_range(pr_peak_locs(2));
%         [dummy,minloc] = min(prior_dist(d,pr_peak_locs(1):pr_peak_locs(2)));
%         n_pr_midpt(d) = minloc+pr_peak_locs(1);
%     else
%         n_pr_down_mp(d) = mp_range(pr_peak_locs(1));
%         n_pr_up_mp(d) = nan;
%         n_pr_midpt(d) = length(mp_range);
%     end
% 
%     n_ad_down_fract(d) = trapz(mean_afterhyp(d,1:n_ad_midpt(d)))/trapz(mean_afterhyp(d,:));
%     n_pr_down_fract(d) = trapz(prior_dist(d,1:n_pr_midpt(d)))/trapz(prior_dist(d,:));
% 
%     d
% end
% 
% cd C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata
% save april_09_1s_hyp_data *fract *down_mp *up_mp
