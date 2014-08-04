clear all
close all
cd C:\WC_Germany\current_injection\persistent_ci

load april_09_1s_ci_data 

% used_norm_mean_dep_traj = norm_mean_afterdep_traj(1:3,501:2500);
% used_norm_mean_hyp_traj = norm_mean_afterhyp_traj(1:3,501:2500);
used_mean_dep_traj = mean_afterdep_traj(:,501:2500);
used_mean_hyp_traj = mean_afterhyp_traj(:,501:2500);

% used_mean_beforedep = mean_beforedep(1:3,:);
% used_mean_afterdep = mean_afterdep(1:3,:);

load seprec_1s_ci_data

% norm_mean_afterdep_traj = [norm_mean_afterdep_traj; used_norm_mean_dep_traj];
% norm_mean_afterhyp_traj = [norm_mean_afterhyp_traj; used_norm_mean_hyp_traj];
mean_afterdep_traj = [mean_afterdep_traj*100; used_mean_dep_traj];
mean_afterhyp_traj = [mean_afterhyp_traj*100; used_mean_hyp_traj];

% mean_beforedep = [mean_beforedep/100; used_mean_beforedep];
% mean_afterdep = [mean_afterdep/100; used_mean_afterdep];

% mp_min = -100;
% mp_max = 10;
% nbins = 500;
% mp_range = linspace(mp_min,mp_max,nbins);
% 
% avg_dist_beforedep = mean(mean_beforedep);
% lci_dist_beforedep = avg_dist_beforedep+std(mean_beforedep)/sqrt(13);
% uci_dist_beforedep = avg_dist_beforedep-std(mean_beforedep)/sqrt(13);
% 
% avg_dist_afterdep = mean(mean_afterdep);
% lci_dist_afterdep = avg_dist_afterdep+std(mean_afterdep)/sqrt(13);
% uci_dist_afterdep = avg_dist_afterdep-std(mean_afterdep)/sqrt(13);
% 
% figure
% plot(mp_range,avg_dist_beforedep,'linewidth',2)
% hold on
% plot(mp_range,lci_dist_beforedep,'--')
% plot(mp_range,uci_dist_beforedep,'--')
% % X = [mp_range fliplr(mp_range)];
% % Y = [lci_dist_beforedep fliplr(uci_dist_beforedep)];
% % fill(X,Y,'b')
% plot(mp_range,avg_dist_afterdep,'r','linewidth',2)
% % Y = [lci_dist_afterdep fliplr(uci_dist_afterdep)];
% % fill(X,Y,'r')
% plot(mp_range,lci_dist_afterdep,'r--')
% plot(mp_range,uci_dist_afterdep,'r--')


%%
% split_points_b = [-65 -60 -69 -54 -60 nan -64 -72 -62 -63 -61];
% split_points_a = [-65 -63 -69 -55 -62 nan -63 -69 -62 -65 -64];
% 
% for i = 1:length(split_points_a)
%     cur = find(mp_range > split_points_b(i),1,'first');
%     if ~isempty(cur) & ~isnan(cur)
%         split_id_b(i) = cur;
%         split_id_a(i) = find(mp_range > split_points_a(i),1,'first');
%     else
%         split_id_a(i) = nan;
%         split_id_b(i) = nan;
%     end
% end
% 
% prob_up_b = trapz(mean_beforedep(:,split_id_b:end),2)./trapz(mean_beforedep,2);
% prob_up_a = trapz(mean_afterdep(:,split_id_b:end),2)./trapz(mean_afterdep,2);
% prob_up_b(6) = [];
% prob_up_a(6) = [];

%% averages in mV
avg_afterdep_traj = mean(mean_afterdep_traj);
uci_afterdep_traj = avg_afterdep_traj+std(mean_afterdep_traj)/sqrt(11);
lci_afterdep_traj = avg_afterdep_traj-std(mean_afterdep_traj)/sqrt(11);

avg_afterhyp_traj = mean(mean_afterhyp_traj);
uci_afterhyp_traj = avg_afterhyp_traj+std(mean_afterhyp_traj)/sqrt(11);
lci_afterhyp_traj = avg_afterhyp_traj-std(mean_afterhyp_traj)/sqrt(11);

plot(tvec,avg_afterdep_traj)
hold on
X = [tvec fliplr(tvec)];
Y = [uci_afterdep_traj fliplr(lci_afterdep_traj)];
fill(X,Y,'b')

figure
plot(tvec,avg_afterhyp_traj,'r')
hold on
X = [tvec fliplr(tvec)];
Y = [uci_afterhyp_traj fliplr(lci_afterhyp_traj)];
fill(X,Y,'r')



%% averages 12s
load all_12s_ci_data
starttime = find(time > 8, 1,'first')+1
endtime = find(time > 24, 1, 'first')-1
new_mean_norm_sweep = new_mean_norm_sweep(:,starttime:endtime);
new_mean_sweep = new_mean_sweep(:,starttime:endtime);

ov_mean_norm_sweep = [mean_norm_sweep; new_mean_norm_sweep];
ov_mean_sweep = [mean_sweep*100; new_mean_sweep];

av_mean_norm_sweep = mean(ov_mean_norm_sweep);
se_mean_norm_sweep = std(ov_mean_norm_sweep)/sqrt(6);
av_mean_sweep = mean(ov_mean_sweep);
se_mean_sweep = std(ov_mean_sweep)/sqrt(6);

t_axis = t_axis';

figure
plot(t_axis,av_mean_sweep)
hold on
X = [t_axis fliplr(t_axis)];
Y = [(av_mean_sweep - se_mean_sweep) fliplr(av_mean_sweep + se_mean_sweep)];
fill(X,Y,'b')
% ylim([-2.5 4.5])
xlim([0 16])
% line([1 1],[-2 -1],'Color','k')
% line([13 13],[-2 -1],'Color','k')
