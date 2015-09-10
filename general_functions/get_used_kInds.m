function used_kInds = get_used_kInds(stim_dims,use_nPix)
% used_kInds = get_used_kInds(stim_dims,use_nPix)
% get the set of indices for filter coefficients that fall in a central
% range of used pixels within the time-embedded stim matrix
% INPUTS:
%   stim_dims: vector of stimulus dimensions with the usual format
%       [time_lags, spatial pixx, spatial_pixy]
%   use_nPix: number of pixels we want to use (assumed to be centered)
% OUTPUTS:
%   used_kInds: vector of used indicies relative to the full time-embedded Xmat

[Xinds,~] = meshgrid(1:stim_dims(2),1:stim_dims(1)); %Xindex associated with each filter element
buffer_pix = floor((stim_dims(2) - use_nPix)/2); %number of pixels to ignore at beginning of stim (half the excess pixels)
cur_use_pix = (1:use_nPix) + buffer_pix; %set of pixels we want to use
used_kInds = find(ismember(Xinds(:),cur_use_pix)); %these are there index values within the time-embedded xmat