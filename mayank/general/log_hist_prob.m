function [n,x] = log_hist_prob(data,range,numBins)

if range(1) <= 0
    range(1) = 0.01;
end

start_pt = log(range(1));
end_pt = log(range(2));


grid = linspace(start_pt,end_pt,numBins);

x = exp(grid);

n = hist(data,x);

n = n/sum(n);