function [avg_fun,std_fun,N_fun,bin_cents,avg_xvals] = ...
    get_nonparametric_relationship(Y_fac,X_fac,n_bins)
%bins X_fac into equi-populated bins, and computes avg and std of Y in
%those bins

if length(Y_fac) ~= length(X_fac)
    error('alignment mismatch');
end

%get edges of equi-spaced percentile bins in x
x_bin_edges = prctile(X_fac,linspace(0,100,n_bins+1));
bin_cents = 0.5*x_bin_edges(1:end-1) + 0.5*x_bin_edges(2:end); %bin centers

[avg_fun,std_fun,N_fun,avg_xvals] = deal(nan(n_bins,1));
for ii = 1:n_bins
   cur_set = find(X_fac >= x_bin_edges(ii) & X_fac < x_bin_edges(ii+1)); %set of data within this x-bin
   N_fun(ii) = sum(~isnan(Y_fac(cur_set))); %number of non-nan Yvals
   avg_fun(ii) = nanmean(Y_fac(cur_set));
   std_fun(ii) = nanstd(Y_fac(cur_set));
   avg_xvals(ii) = nanmean(X_fac(cur_set));
end
