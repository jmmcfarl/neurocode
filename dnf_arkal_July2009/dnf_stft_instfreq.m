function iFreq = dnf_stft_instfreq(y, Fs, window, nFFT)
%
% iFreq = dnf_stft_instfreq(Y, Fs, windowlength, nFFT)
%
%  Y is the band-limited data signal
%  Fs is the sampling frequency of Y
%  windowlength is the length of the data window to consider for each
%     iteration
%  nFFT is the length of the FFT
%
% Copyright (c) 2009, David P. Nguyen <dpnguyen@mit.edu>
%
   
  y = y(:,1);
  
  if ~exist('Fs')
    Fs = 100;
  end
  
  if ~exist('window', 'var'); window = ceil(Fs/5); end
  if ~exist('nFFT', 'var'); nFFT = 512; end
  
  [S, F, T, P] = spectrogram(y, window, window-1, nFFT, Fs, 'yaxis');
  
  if 1 == 0
    figure;
    imagesc(T, F, 10.*log10(abs(P)));
    title(sprintf('%d %d', window, nFFT));
    axis xy;
  end

  % take the maximum peak of the simulation
  [dummy, indmax] = max(P, [], 1);
  iFreq = F(indmax);

  % take the instantaneous frequency of the simulation
  % pad estimates so they equal the right number of points
  L1 = floor(T(1).*Fs) - 1;
  iFreq = [iFreq(1).*ones(L1,1); iFreq(:)];
  L2 = numel(y) - numel(iFreq);
  iFreq = [iFreq(:); iFreq(end).*ones(L2,1)];
  
  % compute time stamps
  dT = T(2)-T(1);  
  T1 = T(1)-dT*L1:dT:T(1)-dT;
  T2 = T(end)+dT:dT:T(end)+dT*L2;
  T = [T1(:); T(:); T2(:)];

  % compute the weighted / central freq
  Pw = P./repmat(sum(P), [size(P,1), 1]);    
  mm = sum(Pw.*repmat(F, [1 size(P,2)]));

  % pad estimates so they equal the right number of points
  mm = [mm(1).*ones(L1,1); mm(:); mm(end).*ones(L2,1)];

  % smooth the estimate
  winhalf = floor(Fs/20);
  win = winhalf*2;
  smwin = gausswin(winhalf*2+1);
  smwin = smwin./sum(smwin);

  iFreq = iFreq(:);
  
  % pad iFreq before smoothing
  iFreq2 = [repmat(iFreq(1), winhalf, 1); iFreq; repmat(iFreq(end), winhalf, 1)];
  iFreq_sm = conv(iFreq2, smwin);
  iFreq_sm = iFreq_sm(win+1:end-win);
  
  iFreq3 = iFreq;
  clear iFreq;
  
  % output the results
  iFreq.f = iFreq3;
  iFreq.f_smooth = iFreq_sm;
  iFreq.f_mean = mm(:);
  iFreq.timestamp = T(:);
  iFreq.srate = Fs;
  iFreq.stft_parms = [window nFFT];
  iFreq.smoothlength = numel(smwin);
  
