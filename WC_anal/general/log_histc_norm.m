function [n,x] = log_histc_norm(data,range,numBins)

if range(1) <= 0
    error('lower range is not positive')
end

start_pt = log10(range(1));
end_pt = log10(range(2));

x = logspace(start_pt,end_pt,numBins);
n = histc(data,x);
n = n/sum(n);