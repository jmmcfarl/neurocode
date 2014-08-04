function poleinfo = dnf_arkal_coef2poleinfo(An)
% poleinfo = dnf_arkal_coef2poleinfo(An)
% 
% Converts the AR coefficients to pole magnitude and pole phase
%
% An is PxN matrix.  Make sure AR coefficeints are by column.
%
% Copyright (c) 2009, David Nguyen <dpnguyen@mit.edu>
  
  
% CONVERT THE AR COEFFICIENTS TO FREQUENCY
  S = size(An);
  P = S(1);
  N = S(2);
  
  tmproots = zeros(P,N);
  
  for nn = 1:N
    tmp = roots(fliplr([1, -An(:,nn)']));
    tmproots(:,nn) = [tmp; zeros(P-length(tmp),1)];
  end
  Fn = -angle(tmproots);
  Mn = abs(tmproots);
  
  %% SORT THE FREQUENCIES
  [Fn,I] = sort(Fn, 1);

  %% SORT THE MAG
  for j = 1:N
    Mn(:,j) = Mn(I(:,j),j);
  end
  
  poleinfo.polemag = 1./Mn;
  poleinfo.polephase = Fn;
  poleinfo.polefreq = 'freq = polephase.*(1/pi).*(Fs/2)';
