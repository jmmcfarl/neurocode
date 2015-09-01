function [up_hist1,up_hist2,down_hist1,down_hist2,up_hist,down_hist] = compare_mp_lfp_states(state_seq1,state_seq2,Fs,f_names)

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



%% compare state duration distributions
dist_range = linspace(0,3,80);

up_hist1 = histc(up_durs1,dist_range);
up_hist1 = up_hist1/sum(up_hist1);
up_hist2 = histc(up_durs2,dist_range);
up_hist2 = up_hist2/sum(up_hist2);

stairs(dist_range,up_hist1), hold on
stairs(dist_range,up_hist2,'r')
t_names = ['G:\WC_Germany\parietal_cortical_2010\compare_states\up_dist_' f_names];
print('-dpng',t_names), close
    
down_hist1 = histc(down_durs1,dist_range);
down_hist1 = down_hist1/sum(down_hist1);
down_hist2 = histc(down_durs2,dist_range);
down_hist2 = down_hist2/sum(down_hist2);

stairs(dist_range,down_hist1), hold on
stairs(dist_range,down_hist2,'r')
t_names = ['G:\WC_Germany\parietal_cortical_2010\compare_states\down_dist_' f_names];
print('-dpng',t_names), close

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

lag_range = linspace(-0.5,0.5,50);
up_hist = histc(up_lags,lag_range);
up_hist = up_hist/sum(up_hist);
down_hist = histc(down_lags,lag_range);
down_hist = down_hist/sum(down_hist);

stairs(lag_range,up_hist), hold on
stairs(lag_range,down_hist,'r')
yl = ylim;
line([0 0],yl,'Color','k')
t_names = ['G:\WC_Germany\parietal_cortical_2010\compare_states\lag_dist_' f_names];
print('-dpng',t_names), close




