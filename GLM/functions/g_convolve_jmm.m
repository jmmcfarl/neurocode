function g = g_convolve_jmm( stim, kern)
%
% Usage: g = g_convolve( stim, kern, <frac> )
%
% Efficient form of g_convolution
% Assumes first term in kernel is 0-latency stim 
% and that k is reversed in time
%
% Created by DAB a long time ago

% Make stim into a row
[L1 L] = size(stim);
if L1 ~= 1
  stim = stim';  L = L1;
end

g = zeros(1,L);
for i = 1:length(kern)
  g = g + kern(i)*stim;
  stim = [0 stim(1:end-1)];
end

% Make into a column
g = g(:);
