function Fn = dnf_arkal_coef2freq(An, Fs)
% Fn = dnf_arkal_coef2freq(An, Fs)
% 
% Converts the AR coefficients to oscillation frequency.
% the pole phase is returned, to convert phase to freq use
%
% An is PxN matrix.  Make sure AR coefficeints are by column.
%
% 
% Copyright (c) 2009, David P. Nguyen <dpnguyen@mit.edu>
  
if ~exist('Fs')
  Fs = 1;
end

% CONVERT THE AR COEFFICIENTS TO FREQUENCY
S = size(An);
P = S(1);
N = S(2);

tmproots = zeros(P,N);

for nn = 1:N
  tmp = roots(fliplr([1, -An(:,nn)']));
  tmproots(:,nn) = [tmp; zeros(P-length(tmp),1)];
end

Fn = angle(tmproots);
Fn = sort(Fn, 1);
Fn = Fn.*(1/pi).*(Fs/2);

  
