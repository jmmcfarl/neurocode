function [smoothed_data] = boxcar_smooth(data,width)

kernel = ones(width,1)/width;
smoothed_data = convn(data,kernel,'same');