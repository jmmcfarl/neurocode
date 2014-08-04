function y = sine_amfm(simopt)
% y = sine_amfm(simopt)
% Create oscillation where amplitude, frequency can change with 
% each sample.  Also, noise can be injected into the signal.
%
% simopt.T = 200;
% simopt.fs = 100;
%%% compute trajectory for frequency
% simopt.ft.f = [40 60];  % 2 steps
% simopt.ft.f = [30 70 40];  % 2 ramps
% simopt.ft.t = [0 100 200];
%%% compute trajectory for amplitude
% simopt.at.a = [1]; % flat envelope
% simopt.at.a = [1 10]; % ramping envelope
% simopt.at.t = [0 200];
%%% compute trajectory for noise
% simopt.nt.n = [0.1];
% simopt.nt.t = [0 200]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2009, David P. Nguyen <dpnguyen@neurostat.mit.edu>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  dt = 1./simopt.fs;
  at = maketimeseries(simopt.at.a, simopt.at.t, dt);
  ft = maketimeseries(simopt.ft.f, simopt.ft.t, dt);
  nt = maketimeseries(simopt.nt.n, simopt.at.t, dt);
  nt = randn(1,length(nt)).*nt;
  
  p = zeros(1,length(ft));
  p(1) = ft(1)*dt;
  for k = 2:length(p)
    p(k) = p(k-1) + ft(k)*dt;
  end
  
  y = at.*sin(2.*pi.*p) + nt;
  
 
  
function y = maketimeseries(a, t, dt)
% a = amplitude points (Mx1)
% t = time points (Mx1)
% dt = delta time (1x1)
%
%
  
  a = a(:);
  t = t(:);
  
  if length(a) == length(t)

    %% in this case, create ramps
    time = 0:dt:t(end);
    clear ti;
    for k = 1:length(t)
      ind = find(time <= t(k));
      ti(k) = ind(end);
    end
    
    L = diff(ti);
    y = [];
    for k = 1:length(L)
      y = [y, linspace(a(k), a(k+1), L(k))];
    end
    
  elseif length(a) == (length(t)-1)
    
    %% in this scenario, create steps
    time = 0:dt:t(end);
    clear ti;
    for k = 1:length(t)
      ind = find(time <= t(k));
      ti(k) = ind(end);
    end
    
    L = diff(ti);
    y = [];
    for k = 1:length(L)
      y = [y, linspace(a(k), a(k), L(k))];
    end
        
  end
  
