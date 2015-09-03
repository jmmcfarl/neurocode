function [n,x] = log_hist(data,range,numBins)

x = logspace(log10(range(1)),log10(range(2)),numBins);

n = histc(data,x);

n = n/sum(n);