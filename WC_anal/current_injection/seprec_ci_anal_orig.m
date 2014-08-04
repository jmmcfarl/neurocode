clear all

depol_dir{1} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-05-24_CWC_LFP_B';
depol_dir{2} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-05-31_CWC_LFP';
depol_dir{3} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-06-04_CWC_LFP_B';
depol_dir{4} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-03_CWC_LFP_A';
depol_dir{5} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-03_CWC_LFP_B';
depol_dir{6} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-04_CWC_LFP_B';
depol_dir{7} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-04_CWC_LFP_D';
depol_dir{8} = 'C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-07-05_CWC_LFP';

prior_dir{1} = 'C:\WC_Germany\EC_MPasci\A2007_05_24_CWC_LFP_B';
prior_dir{2} = 'C:\WC_Germany\EC_MPasci\A2007_05_31_CWC_LFP';
prior_dir{3} = 'C:\WC_Germany\EC_MPasci\A2007_06_04_CWC_LFP_B';
prior_dir{4} = 'C:\WC_Germany\EC_MPasci\A2007_07_03_CWC_LFP_A';
prior_dir{5} = 'C:\WC_Germany\EC_MPasci\A2007_07_03_CWC_LFP_B';
prior_dir{6} = 'C:\WC_Germany\EC_MPasci\A2007_07_04_CWC_LFP_B';
prior_dir{7} = 'C:\WC_Germany\EC_MPasci\A2007_07_04_CWC_LFP_D';
prior_dir{8} = 'C:\WC_Germany\EC_MPasci\A2007_07_05_CWC_LFP';

save_name{1} = '2007-5-24';
save_name{2} = '2007-5-31';
save_name{3} = '2007-6-4';
save_name{4} = '2007-7-3_A';
save_name{5} = '2007-7-3_B';
save_name{6} = '2007-7-4_B';
save_name{7} = '2007-7-4_D';
save_name{8} = '2007-7-5';





mp_min = -1.0;
mp_max = 0.1;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;

% prior_dist = zeros(8,length(mp_range));
% for d = 1:8
%     
%     load(prior_dir{d})
%     data_string = prior_dir{d}(25:end);
%     data_string = strcat(data_string,'_MP');
%     eval(['cur_data = ' data_string ';'])
%     prior_dist(d,:) = gpkde(cur_data,.02,[mp_min;mp_max;nbins]);
%     prior_dist_beg(d,:) = gpkde(cur_data(1:3e6),.02,[mp_min;mp_max;nbins]);
%     prior_dist_end(d,:) = gpkde(cur_data(end-3e6:end),.02,[mp_min;mp_max;nbins]);
% 
%     mean_prior(d) = mean(cur_data);
%     std_prior(d) = std(cur_data);
%     
%    d
% end
% 

load old_1s_dep_data_minus_spike

for d = 1:8
    
    load(depol_dir{d})
    data = data_minus_spike{d};

    send_times = find(time==1.9999);
    sbeg_times = find(time==0);

    num_sweeps = length(sbeg_times);
    tvec = 0:1/Fsd:(2-1/Fsd);
    sweep_mat = zeros(num_sweeps,length(tvec));
    for i = 1:num_sweeps-1
        curdata = data(sbeg_times(i):send_times(i));
        sweep_mat(i,:) = downsample(curdata,dsf);
    end

%     figure
%     imagesc(sweep_mat); ylim([55 65])
%     t_name = ['C:\WC_Germany\current_injection\persistent_ci\ci_mat_seprec\' save_name{d}];
%     print('-dpng',t_name); close
    
%     tdep_dist = zeros(length(tvec),length(mp_range));
%     thyp_dist = zeros(length(tvec),length(mp_range));
%     for t = 1:length(tvec)
%         tdep_dist(t,:) = gpkde(sweep_mat(61:120,t),.02,[mp_min;mp_max;nbins]);
%         thyp_dist(t,:) = gpkde(sweep_mat(1:60,t),.02,[mp_min;mp_max;nbins]);
%     end
    
%     norm_sweep_mat = (sweep_mat - mean_prior(d))/std_prior(d);
    
    dep_end = find(tvec > 1.6,1,'first');
    dep_start = find(tvec > 0.5,1,'first');
%     norm_mean_afterdep_traj(d,:) = mean(norm_sweep_mat(61:120,:));
    mean_afterdep_traj(d,:) = mean(sweep_mat(61:120,:));
%     mean_afterdep(d,:) = mean(tdep_dist(dep_end:end,:));
%     mean_beforedep(d,:) = mean(tdep_dist(1:dep_start,:));
%     mean_duringdep(d,:) = mean(tdep_dist(dep_start+100:dep_end-100,:));
    
    mean_afterhyp_traj(d,:) = mean(sweep_mat(1:60,:));
%     norm_mean_afterhyp_traj(d,:) = mean(norm_sweep_mat(1:60,:));
%     mean_afterhyp(d,:) = mean(thyp_dist(dep_end:end,:));
%     mean_beforehyp(d,:) = mean(thyp_dist(1:dep_start,:));
%     mean_duringhyp(d,:) = mean(thyp_dist(dep_start+100:dep_end-100,:));

    avg_amp_before_dep(d) = mean(mean(sweep_mat(61:120,1:dep_start-100)));
    avg_amp_during_dep(d) = mean(mean(sweep_mat(61:120,dep_start:dep_end)));
    avg_amp_before_hyp(d) = mean(mean(sweep_mat(1:60,1:dep_start-100)));
    avg_amp_during_hyp(d) = mean(mean(sweep_mat(1:60,dep_start:dep_end)));

    
%         plot(mp_range*100,prior_dist(d,:),'linewidth',2)
%     hold on
%     plot(mp_range*100,prior_dist_beg(d,:))
%     plot(mp_range*100,prior_dist_end(d,:),'--')
%     plot(mp_range*100,mean_afterdep(d,:),'r')
%     plot(mp_range*100,mean_beforedep(d,:),'r--')
%     plot(mp_range*100,mean_afterhyp(d,:),'g')
%     plot(mp_range*100,mean_beforehyp(d,:),'g--')
%     xlim([-90 -20])
%     t_name = ['C:\WC_Germany\current_injection\persistent_ci\ci_mp_dist_seprec\' save_name{d}];
%     print('-dpng',t_name);close
    
    plot(tvec,mean_afterdep_traj(d,:),'r')
    hold on
    plot(tvec,mean_afterhyp_traj(d,:),'g')
    t_name = ['C:\WC_Germany\current_injection\persistent_ci\ci_mp_traj_seprec\' save_name{d}];
    print('-dpng',t_name);close

    
d
end

cd C:\WC_Germany\current_injection\persistent_ci
% save seprec_1s_ci_data tvec mean_* avg_amp* norm_mean* prior*
save seprec_1s_ci_data tvec mean_* 


