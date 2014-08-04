function iFreq = dnf_hilbert_instfreq(y, Fs, smoothwin)
%
% iFreq = dnf_hilbert_instfreq(Y, Fs, smoothwinlength)
%
%  Y is the band-limited data signal
%  Fs is the sampling frequency of Y
%  smoothwinlength is the length of gaussian smoothing window
%     Default length: floor(Fs/80)
%
% 
% Copyright (c) 2009, David P. Nguyen <dpnguyen@mit.edu>
%
  
  if ~exist('smoothwin')
    smoothwin = floor(Fs/80);
  end
  
  dt = 1/Fs;
  
  h = hilbert(y(:));
  a = angle(h);
  a = unwrap(a);
  a = abs(a);
  
  iFreq.f = diff(a).*Fs./(2*pi);
  iFreq.f = [iFreq.f(1); iFreq.f(:)];
  
  %% smooth the angle
  winhalf = max(smoothwin, 5);
  smwin = gausswin(winhalf*2+1);
  smwin = smwin./sum(smwin);
  
  %% pad a before smoothing
  a = [a(1).*ones(winhalf,1); a; a(end).*ones(winhalf,1)];
  
  a = conv(a, smwin);
  a = a(winhalf*2+1:end-winhalf*2);
 
  f = diff(a).*Fs./(2*pi);
  iFreq.f_smooth = [f(1); f(:)];

  if 1 == 0
    N = numel(y);
    [pp, ii] = findpeaks(f);
    ii = [1; ii(:); N];
    pp = [f(1); pp(:); f(N)];
    f = interp1(ii, pp, 1:N);
  end
  
