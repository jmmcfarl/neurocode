function [smoothed_gram] = smooth_gram(gram, time_sigma, freq_sigma)

x_axis = -3*time_sigma:3*time_sigma;
y_axis = -3*freq_sigma:3*freq_sigma;

[xgrid,ygrid] = meshgrid(x_axis,y_axis);

g_kern = exp(-(xgrid.^2/(2*time_sigma^2) + ygrid.^2/(2*freq_sigma^2)));
g_kern = g_kern/sum(sum(g_kern));

smoothed_gram = conv2(gram,g_kern');

x_chop = (length(x_axis)-1)/2;
y_chop = (length(y_axis)-1)/2;

smoothed_gram(1:x_chop,:) = [];
smoothed_gram(end-x_chop+1:end,:) = [];
smoothed_gram(:,1:y_chop) = [];
smoothed_gram(:,end-y_chop+1:end) = [];