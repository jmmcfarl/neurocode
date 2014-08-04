function [] = print_stats(data)

n_pts = sum(~isnan(data));
data_mean = nanmean(data);
data_sem = nanstd(data)/sqrt(n_pts);
rdata_sem = robust_std_dev(data);
data_median = nanmedian(data);
fprintf('N:%d, Mean: %.4f  SEM: %.4f\n',n_pts,data_mean,data_sem);
fprintf('Median: %.4f  rSEM: %.4f\n',data_median,rdata_sem);