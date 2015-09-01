function [mua_times] = mua_extract(lfp_signal,threshold)

niqf = 2016/2;
[b,a] = butter(2,150/niqf,'high');

lfp_f = filter(b,a,lfp_signal);

lfp_f = -zscore(lfp_f);

surr_dat = lfp_f;

surr_dat(lfp_f > threshold) = 1;
surr_dat(lfp_f < threshold) = 0;
dsurr_dat = [0;diff(surr_dat)];

mua_times = find(dsurr_dat == 1);




