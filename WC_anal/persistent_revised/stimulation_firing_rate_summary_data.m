%stimulation rates, and prior rates, excluding 2009-4-5

stim_rate_avg = [8.96 7.65 9.28 14.1 8.8 11.7 12.03 4.78 4.48 5.1 4.33];
ov_rate_avg = [1.63 2.92 4.96 2.9 1.83 4.47 5.3 1.13 3.95 3.1 4.73];
up_rate_avg = [3.10 4.17 7.98 6.60 4.70 7.71 9.70 5.53 8.61 5.36 9.21];
    
fold = stim_rate_avg./ov_rate_avg;
fold_up = stim_rate_avg./up_rate_avg;
perc_increase = (stim_rate_avg-ov_rate_avg)./ov_rate_avg*100;
pers_increase_up = (stim_rate_avg-up_rate_avg)./up_rate_avg*100;

mean(perc_increase)
std(perc_increase)/sqrt(length(perc_increase))
