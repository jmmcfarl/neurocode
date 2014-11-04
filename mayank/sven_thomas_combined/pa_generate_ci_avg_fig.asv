clear all
close all
cd C:\WC_Germany\persistent_9_27_2010\

load ./april09_1s_ci_data.mat
ap_stim_rate = stim_rate;
ap_non_stim_rate = non_stim_rate;

keep_inds = 501:2500;
used_mean_dep_traj = mean_afterdep_traj(:,keep_inds);
% used_mean_hyp_traj = mean_afterhyp_traj(:,keep_inds);

load ./seprec_1s_ci_data.mat

all_stim_rate = [stim_rate ap_stim_rate];

mean_afterdep_traj = [mean_afterdep_traj*100; used_mean_dep_traj];
% mean_afterhyp_traj = [mean_afterhyp_traj*100; used_mean_hyp_traj];

onesec_stim_data = [4 6 8 11 12 13 14 16 22 18 19];
cd C:\WC_Germany\sven_thomas_combined
load ./spike_rate_data_fin
backgnd_rate_onesec = backgnd_rate(onesec_stim_data);

percent_increase = (all_stim_rate - backgnd_rate_onesec)./backgnd_rate_onesec*100;

%% averages in mV
n_cells = 11;

avg_afterdep_traj = mean(mean_afterdep_traj);
uci_afterdep_traj = avg_afterdep_traj+std(mean_afterdep_traj)/sqrt(n_cells);
lci_afterdep_traj = avg_afterdep_traj-std(mean_afterdep_traj)/sqrt(n_cells);

junct_pot = 7; %junction potential in mV

tvec = (1:size(mean_afterdep_traj,2))/Fsd;
figure
plot(tvec,avg_afterdep_traj-junct_pot)
hold on
X = [tvec fliplr(tvec)];
Y = [uci_afterdep_traj-junct_pot fliplr(lci_afterdep_traj)-junct_pot];
fill(X,Y,'b')
s = zeros(size(tvec));
s(tvec >= 0.5 & tvec <= 1.5) = 1;
plot(tvec,s*10-80,'k')

%% averages 12s

twelve_data = [6 8 18 19 20];
cd C:\WC_Germany\persistent_9_27_2010
load ./all_12s_CI_data_minus_spk.mat

backgnd_rate_twelvesec = backgnd_rate(twelve_data);

percent_increase_twelve = (stim_rate - backgnd_rate_twelvesec)./backgnd_rate_twelvesec*100;

starttime = find(t_axis > 8, 1,'first')+1;
endtime = find(t_axis > 24, 1, 'first')-1;
new_mean_sweep = new_mean_sweep(:,starttime:endtime);

ov_mean_sweep = [mean_sweep*100; new_mean_sweep];

av_mean_sweep = mean(ov_mean_sweep);
se_mean_sweep = std(ov_mean_sweep)/sqrt(5);

t_axis = (1:size(ov_mean_sweep,2))/Fsd;
figure
plot(t_axis,av_mean_sweep-junct_pot)
hold on
X = [t_axis fliplr(t_axis)];
Y = [(av_mean_sweep - se_mean_sweep)-junct_pot fliplr(av_mean_sweep + se_mean_sweep)-junct_pot];
fill(X,Y,'b')
% ylim([-2.5 4.5])
xlim([0 16])
% line([1 1],[-2 -1],'Color','k')
% line([13 13],[-2 -1],'Color','k')
s = zeros(size(t_axis));
s(t_axis >= 1 & t_axis <= 13) = 1;
plot(t_axis,s*10-60,'k')


%% averages 12s with new MEC recs
cd C:\WC_Germany\sven_thomas_combined
load new_mec_ci_recs 
starttime = find(t_axis > 8, 1,'first')+1;
endtime = find(t_axis > 24, 1, 'first')-1;
new_mean_sweep2 = new_mean_sweep(:,starttime:endtime);

twelve_data = [6 8 18 19 20];
cd C:\WC_Germany\persistent_9_27_2010
load ./all_12s_CI_data_minus_spk.mat

new_mean_sweep = new_mean_sweep(:,starttime:endtime);

ov_mean_sweep = [mean_sweep*100; new_mean_sweep; new_mean_sweep2];

av_mean_sweep = mean(ov_mean_sweep);
se_mean_sweep = std(ov_mean_sweep)/sqrt(8);

junct_pot = 7;

t_axis = (1:size(ov_mean_sweep,2))/Fsd;
figure
plot(t_axis,av_mean_sweep-junct_pot)
hold on
X = [t_axis fliplr(t_axis)];
Y = [(av_mean_sweep - se_mean_sweep)-junct_pot fliplr(av_mean_sweep + se_mean_sweep)-junct_pot];
fill(X,Y,'b')
% ylim([-2.5 4.5])
xlim([0 16])
% line([1 1],[-2 -1],'Color','k')
% line([13 13],[-2 -1],'Color','k')
s = zeros(size(t_axis));
s(t_axis >= 1 & t_axis <= 13) = 1;
plot(t_axis,s*10-60,'k')


