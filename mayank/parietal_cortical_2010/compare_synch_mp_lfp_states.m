function [up_hist,down_hist,lag_hist,rob_stats] = compare_synch_mp_lfp_states(state_seq1,state_seq2,desynch_times,Fs,f_names)
addpath('G:\Code\LIBRA_27aug09')

up_trans1 = find(state_seq1(1:end-1) == 1 & state_seq1(2:end)==2);
down_trans1 = find(state_seq1(1:end-1)==2 & state_seq1(2:end)==1);

up_trans2 = find(state_seq2(1:end-1) == 1 & state_seq2(2:end)==2);
down_trans2 = find(state_seq2(1:end-1)==2 & state_seq2(2:end)==1);

up_trans1(up_trans1 > down_trans1(end)) = [];
down_trans1(down_trans1 < up_trans1(1)) = [];

up_trans2(up_trans2 > down_trans2(end)) = [];
down_trans2(down_trans2 < up_trans2(1)) = [];

up_durs1 = (down_trans1 - up_trans1)/Fs;
down_durs1 = (up_trans1(2:end) - down_trans1(1:end-1))/Fs;

up_durs2 = (down_trans2 - up_trans2)/Fs;
down_durs2 = (up_trans2(2:end) - down_trans2(1:end-1))/Fs;

desynch_ups1 = [];
desynch_ups2 = [];
for i = 1:size(desynch_times,1)
    desynch_ups1 = [desynch_ups1; ...
        find(up_trans1/Fs > desynch_times(i,1) & up_trans1/Fs < desynch_times(i,2))];
    desynch_ups2 = [desynch_ups2; ...
        find(up_trans2/Fs > desynch_times(i,1) & up_trans2/Fs < desynch_times(i,2))];
end

desynch_downs1 = [];
desynch_downs2 = [];
for i = 1:size(desynch_times,1)
    desynch_downs1 = [desynch_downs1; ...
        find(down_trans1/Fs > desynch_times(i,1) & down_trans1/Fs < desynch_times(i,2))];
    desynch_downs2 = [desynch_downs2; ...
        find(down_trans2/Fs > desynch_times(i,1) & down_trans2/Fs < desynch_times(i,2))];
end

up_trans1(desynch_ups1) = [];
down_trans1(desynch_downs1) = [];
up_trans2(desynch_ups2) = [];
down_trans2(desynch_downs2) = [];

up_durs1(desynch_ups1) = [];
up_durs2(desynch_ups2) = [];
down_durs1(desynch_downs1) = [];
down_durs2(desynch_downs2) = [];

%% compare state duration distributions
dist_range = 0:10/Fs:3;

up_hist(1,:) = histc(up_durs1,dist_range);
up_hist(1,:) = up_hist(1,:)/sum(up_hist(1,:));
up_hist(2,:) = histc(up_durs2,dist_range);
up_hist(2,:) = up_hist(2,:)/sum(up_hist(2,:));

% stairs(dist_range,up_hist(1,:)), hold on
% stairs(dist_range,up_hist(2,:),'r')
% t_names = ['G:\WC_Germany\parietal_cortical_2010\compare_states\up_dist_' f_names];
% print('-dpng',t_names), close
    
down_hist(1,:) = histc(down_durs1,dist_range);
down_hist(1,:) = down_hist(1,:)/sum(down_hist(1,:));
down_hist(2,:) = histc(down_durs2,dist_range);
down_hist(2,:) = down_hist(2,:)/sum(down_hist(2,:));

% stairs(dist_range,down_hist(1,:)), hold on
% stairs(dist_range,down_hist(2,:),'r')
% t_names = ['G:\WC_Germany\parietal_cortical_2010\compare_states\down_dist_' f_names];
% print('-dpng',t_names), close

%% compute relative timing
up_lags = zeros(size(up_trans1));
for i = 1:length(up_trans1)
    [dummy,near_up] = min(abs(up_trans1(i)-up_trans2));
    up_lags(i) = up_trans1(i)-up_trans2(near_up);
end
up_lags = up_lags/Fs;

down_lags = zeros(size(down_trans1));
for i = 1:length(down_trans1)
    [dummy,near_down] = min(abs(down_trans1(i)-down_trans2));
    down_lags(i) = down_trans1(i)-down_trans2(near_down);
end
down_lags = down_lags/Fs;

% lag_range = linspace(-0.2,0.2,50);
% lag_hist(1,:) = histc(up_lags,lag_range);
% lag_hist(1,:) = lag_hist(1,:)/sum(lag_hist(1,:));
% lag_hist(2,:) = histc(down_lags,lag_range);
% lag_hist(2,:) = lag_hist(2,:)/sum(lag_hist(2,:));

lag_range = linspace(-0.5,0.5,400);
lag_hist(1,:) = ksdensity(up_lags,lag_range);
lag_hist(2,:) = ksdensity(down_lags,lag_range);

% stairs(lag_range,lag_hist(1,:)), hold on
% stairs(lag_range,lag_hist(2,:),'r')
% yl = ylim;
% line([0 0],yl,'Color','k')
% t_names = ['G:\WC_Germany\parietal_cortical_2010\compare_states\lag_dist_' f_names];
% print('-dpng',t_names), close

%% now estimate robust mean and variance of each distribution
n_ups = length(up_lags);
n_downs = length(down_lags);
[rew,raw]=mcdcov(up_lags,'plots',0,'h',round(n_ups*0.75));
rob_stats.mean_up_lag = rew.center;
rob_stats.var_up_lag = rew.cov;

[rew,raw]=mcdcov(down_lags,'plots',0,'h',round(n_downs*0.75));
rob_stats.mean_down_lag = rew.center;
rob_stats.var_down_lag = rew.cov;




