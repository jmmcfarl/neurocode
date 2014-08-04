function [sm_data] = hsmm_uds_gsmooth(data,sigma,kernlen)

%1d gaussian smoothing function

%inputs: 
%   data: vector to be smoothed
%   sigma: sigma of gaussian kernel
%   kernlen: number of sigma to extend kernel in either direction
%outputs:
%   sm_data: smoothed data

if nargin < 3
    kernlen = 3; %default
end

if size(data,1) > size(data,2) 
    data = data';
end

kern_grid = [-kernlen*sigma:kernlen*sigma]; %grid over which kernel is defined
gauss_kern = exp(-(kern_grid).^2/(2*sigma^2)); %gaussian smoothing kernel
gauss_kern = gauss_kern/sum(gauss_kern); %normalize to unit area
kern_size = length(gauss_kern);

if kern_size >= length(data)
    error('Data not long enough for smoothing kernel');
end

sm_data = conv(data,gauss_kern);

%clip convolution artifact
sm_data(end-(kern_size-1)/2+1:end) = [];
sm_data(1:(kern_size-1)/2) = [];

%find all points within a half kernel length of an edge
boundary_points = 1:kernlen*sigma;
half_kern_length = kernlen*sigma+1;

%use assymetric smoothing at Boundaries
for i = 1:length(boundary_points)
   
    rescaled_kern = gauss_kern(half_kern_length-i+1:end); %chop and rescale smoothing kernel to unit sum
    rescaled_kern = rescaled_kern/sum(rescaled_kern);
    
    sm_data(i) = data(1:half_kern_length-1+i)*rescaled_kern';
    sm_data(end-i+1) = data(end-half_kern_length-i+2:end)*flipud(rescaled_kern');

end


