function sgram = dnf_arkal_specgram(An,varargin)
% specgram = dnf_arkal_specgram(An,varargin)
%
% Input:
%   An is output of dnf_arkal.  A matrix of AR coefficients
%   in column order, where each column is a time instance.
%   The size of the column is P (the order of the AR model) 
%   
% Options, defaults in brackets:
%  'fs' = sampling freq [1]
%  'index' = vector of discrete time indices to do decomp
%  'psdres' = number of points to eval on freq axis
%
% 
% Copyright (c) 2009, David P. Nguyen <dpnguyen@mit.edu>
      
  An = squeeze(An);
  sa = size(An);
  P = sa(1);
  N = sa(2);

  opt.fs = 1;
  opt.index = 1:N;
  opt.psdres = 100;
  opt.dospecgram = 0;
  opt.timestamp = [];
  
  opt = parsevarargin(varargin, opt);  

  if length(opt.timestamp) == 0
    opt.timestamp = 1:N;
  end
  
  N = length(opt.index);
  N2 = sa(2);
  
  omega = linspace(0.00001, pi, opt.psdres);
  omega = omega(:);
  
  sgram.specgram = zeros(length(omega),N);
  sgram.w = omega.*(opt.fs./2./pi);
  sgram.timestamp = opt.timestamp;
  sgram.fs = opt.fs;
  sgram.plot = 'imagesc(sg.timestamp, sg.w, sg.specgram);';
  
  for n=1:N   
    nn = opt.index(n);
    cp = [1, -An(:,nn)'];
    cp = cp(:)';
    
    for wi = 1:length(omega)
      w = omega(wi);
      z = exp(j*w);
   
      ad1 = z.^(-(0:P));
      ad1 = ad1(:)';
      ad2 = z.^(0:P);
      ad2 = ad2(:)';
        
      sgram.specgram(wi,nn) = 1./(sum(ad1.*cp)*sum(ad2.*cp));    
    end
  end
  
  
  
  
