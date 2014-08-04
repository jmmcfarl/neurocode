function [slicemat] = padfiltermat(k,spixids,sdim)
% USAGE: [slicemat] = padfiltermat(k,spixids,sdim)
%   takes subset of pixel ids and constructs square RF
slice          = zeros(sdim,sdim); 
slice(spixids) = k; 
slicemat       = reshape(slice,sdim,sdim); 
end

