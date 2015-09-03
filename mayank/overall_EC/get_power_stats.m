function [pow_stats,peak_pow,tot_pow,peak_freq] = get_power_stats(P,f,f_range)

used_freqs = find(f > f_range(1) & f < f_range(2));
df = f(2)-f(1);
[peak_pow,peak_freq] = max(P(:,used_freqs),[],2);
peak_freq = f(used_freqs(peak_freq));
tot_pow = df*trapz(P(:,used_freqs),2);

pow_stats.avg_lpow = mean(log(peak_pow));
pow_stats.std_lpow = std(log(peak_pow));
pow_stats.avg_pfreq = mean(peak_freq);
pow_stats.std_pfreq = std(peak_freq);
pow_stats.avg_tpow = mean(tot_pow);
pow_stats.std_tpow = std(tot_pow);