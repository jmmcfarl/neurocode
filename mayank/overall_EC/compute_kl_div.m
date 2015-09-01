function [sep_info] = compute_kl_div(data,state_seq)

min_cond1 = min(data(state_seq==1));
min_cond2 = min(data(state_seq==2));
max_cond1 = max(data(state_seq==1));
max_cond2 = max(data(state_seq==2));

amp_grid = linspace(max([min_cond1 min_cond2]),min([max_cond1 max_cond2]),200);
data_dist_1 = ksdensity(data(state_seq==1),amp_grid);
data_dist_2 = ksdensity(data(state_seq==2),amp_grid);

kl1 = sum(data_dist_1.*log(data_dist_1./data_dist_2));
kl2 = sum(data_dist_2.*log(data_dist_2./data_dist_1));

sep_info = (kl1+kl2)/2;