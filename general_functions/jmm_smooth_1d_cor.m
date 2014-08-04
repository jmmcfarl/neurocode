function [sm_data] = jmm_smooth_1d_cor(data,sig,windur)

%1d gaussian smoothing function
%inputs: 
%   data: vector to be smoothed
%   sig: sigma of gaussian kernel
%   windur: number of sigma to extend kernel in either direction
%outputs:
%   sm_data: smoothed data

if nargin > 2
    kern_length = windur;
else
    kern_length = 3; %default
end

if size(data,1) > size(data,2) 
    data = data';
end

if sig == 0
    sm_data = data;
    return
end

kern_grid = [-kern_length*ceil(sig):kern_length*ceil(sig)];
gauss_kern = exp(-(kern_grid).^2/(2*sig^2));
gauss_kern = gauss_kern/sum(gauss_kern); %normalize to unit area
kern_size = length(gauss_kern);

if kern_size >= length(data)
    disp('Error data not long enough for smoothing kernel')
    return
end

sm_data = conv(data,gauss_kern);

%clip convolution artifact
sm_data(end-(kern_size-1)/2+1:end) = [];
sm_data(1:(kern_size-1)/2) = [];

%find all points within a half kernel length of an edge
boundary_points = 1:kern_length*ceil(sig);

half_kern_length = kern_length*ceil(sig)+1;

for i = 1:length(boundary_points)
   
    rescaled_kern = gauss_kern(half_kern_length-i+1:end);
    rescaled_kern = rescaled_kern/sum(rescaled_kern);
    
    sm_data(i) = data(1:half_kern_length-1+i)*rescaled_kern';
    sm_data(end-i+1) = data(end-half_kern_length-i+2:end)*flipud(rescaled_kern');

end


