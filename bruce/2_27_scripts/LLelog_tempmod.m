function [ll,grad] = LLelog_tempmod(K,X,S,spatial_mod,S_spatial_mod,lamrange,lamrange2,llist, hold_const, nlxrange)
%
% Usage: [ll,grad] = LLelog2(K,X,S,<lamrange>,<lamrange2>,<llist>, <hold_const>)
%
% calculates log-likelihood and gradient for non-linearity: f(x) = log(1+e(x))
%
% lamrange2 is a slope-slope penalty, llist is |k| penalty (but need to impose constr)
%
%  lamrange = [lambda1 range1_beginning range1_end]
%             [lambda2 range2_beginning range2_end
%             [...
%
% if range(1) is less then zero, assume circular with slope
%
% if range(2) is less than zeros, then not range but endpoints (for 2-D etc)

HCmult = 1e6;
if (nargin < 4) || isempty(lamrange)
  els = [];
else
  els = lamrange(:,1);
  if sum(els) > 0
    r1s = lamrange(:,2);
    r2s = lamrange(:,3);
  else
    els = [];  r1s = [];
  end
%   circ = r1s < 0;
  r1s = abs(r1s);
end
if (nargin < 5) || isempty(lamrange2)
  els2 = [];
else
  els2 = lamrange2(:,1);
  if sum(els2) > 0
    r1s2 = lamrange2(:,2);
    r2s2 = lamrange2(:,3);
  else
    els2 = [];  r1s2 = [];
  end
%   circ2 = r1s2 < 0;
  r1s2 = abs(r1s2);
end
if nargin < 6
  llist = [];
end
if (nargin < 7) || isempty(hold_const)
  hold_const = [];
  hcpun = 0;
else
  %K(hold_const(1,:)) = hold_const(2,:);
  hcgrad = make_row(K(hold_const(1,:)))-hold_const(2,:);
  hcpun = HCmult * sum(hcgrad.^2);
  %hold_const = hold_const(1,:);
end

k = K(1:end-1);
b = K(end);

kx = X*k+b+spatial_mod;
ks = S*k+b+S_spatial_mod;

%% PUT CEILING ON kx and ks so it doesnt blow up
kx(kx > 50) = 50;
ks(ks > 50) = 50;

ekx = exp(kx);
r = log(1+ekx);

eks = exp(ks);
rs = log(1+eks);

rs(rs < 1e-10) = 1e-10;

% Standard point-process likelihood
ll = sum(log(rs)) - sum(r);

wg = ekx./(1+ekx);
ers = (eks./(1+eks))./rs;
grad = ers'*S-wg'*X;             % grad w.r.t. k
grad(end+1) = sum(ers)-sum(wg);  % grad w.r.t. b

% Add penalty terms from slope
for i = 1:length(els)
    %     if r2s(i) > 0
    range = (r1s(i):r2s(i))';
    chunk = k(range)';  %'
    xchunk = nlxrange(range);
    
    dy_dx = (chunk(2:end)-chunk(1:end-1))./(xchunk(2:end)-xchunk(1:end-1));
%     ll = ll - els(i) * sum(dy_dx.^2) - hcpun;
            ll = ll - els(i) * sum((chunk(2:end) - chunk(1:end-1)).^2) - hcpun;
    
            grad(range(2:end)) = grad(range(2:end)) - 2*els(i)*(chunk(2:end) - chunk(1:end-1));
            grad(range(1:end-1)) = grad(range(1:end-1)) - 2*els(i)*(chunk(1:end-1) - chunk(2:end));
    
%     grad(range(2:end)) = grad(range(2:end)) - 2*els(i)*(chunk(2:end) - chunk(1:end-1))./(xchunk(2:end)-xchunk(1:end-1));
%     grad(range(1:end-1)) = grad(range(1:end-1)) - 2*els(i)*(chunk(1:end-1) - chunk(2:end))./(xchunk(1:end-1) - xchunk(2:end));
    
    
    %     else
    %         a = r1s(i);  b = -r2s(i);
    %         ll = ll - els(i) * (k(a)-k(b))^2 - hcpun;
    %         gradum = 2*els(i) * (k(a)-k(b));
    %         grad(a) = grad(a) - gradum;  % corrected from +gradnum 2/22/11
    %         grad(b) = grad(b) + gradum;  % corrected from -gradnum 2/22/11
    %         if circ(i) > 0
    %             disp('circ no work, in this instance')
    %         end
    %     end
end

% Add penalty terms from slope-slope
for i = 1:length(els2)
    range = r1s2(i):r2s2(i);
    chunk = k(range)'; %'
    xchunk = nlxrange(range);
    
    xchunk2 = [xchunk(1)-1 xchunk xchunk(end)+1]; %constant slope off front and back
    	  chunk2 = [2*chunk(1)-chunk(2) chunk 2*chunk(end)-chunk(end-1)];%Makes constant slope off the front, and 0 slope off the back
%     chunk2 = [chunk(1)-1 chunk chunk(end)+1];
    
%     dy_dx2 = (chunk2(3:end)-chunk2(2:end-1))./(xchunk2(3:end)-xchunk2(2:end-1));
%     dy_dx1 = (chunk2(2:end-1)-chunk2(1:end-2))./(xchunk2(2:end-1)-xchunk2(1:end-2));
%     xchunk_mids = (xchunk2(2:end)+xchunk2(1:end-1))/2;
%     d2y_d2x = (dy_dx2-dy_dx1)./diff(xchunk_mids);
%     
%     delta = (xchunk2(2:end-1)-xchunk2(3:end)).*xchunk2(1:end-2).^2 + (xchunk2(3:end)-xchunk2(1:end-2)).*xchunk2(2:end-1).^2 + (xchunk2(1:end-2)-xchunk2(2:end-1)).*xchunk2(3:end).^2;
%     fds = 2*((xchunk2(2:end-1)-xchunk2(3:end)).*chunk2(1:end-2) + (xchunk2(3:end)-xchunk2(1:end-2)).*chunk2(2:end-1) + (xchunk2(1:end-2)-xchunk2(2:end-1)).*chunk2(3:end));
    
%     ll = ll - els2(i) * sum(d2y_d2x.^2) - hcpun;
    
      ll = ll - els2(i) * sum((2*chunk2(2:end-1) - chunk2(1:end-2) - chunk2(3:end)).^2) - hcpun;
        
    grad(range(3:end-2)) = grad(range(3:end-2)) - 2*els2(i) * (6*chunk2(4:end-3) - 4*chunk2(3:end-4) - 4*chunk2(5:end-2) + chunk2(2:end-5)+chunk2(6:end-1));
    grad(range(1)) = grad(range(1)) - 2*els2(i) * (chunk2(2) + chunk2(4) - 2*chunk2(3));
    grad(range(end)) = grad(range(end)) - 2*els2(i) * (chunk2(end-1) + chunk2(end-3) - 2*chunk2(end-2));
    grad(range(2)) = grad(range(2)) - 2*els2(i) * (5*chunk2(3) - 2*chunk2(2) - 4*chunk2(4) + chunk2(5) );
    grad(range(end-1)) = grad(range(end-1)) - 2*els2(i) * (5*chunk2(end-2) - 2*chunk2(end-1) - 4*chunk2(end-3) + chunk2(end-4) );
    
%     fds2 = [0 fds 0];
%     delta2 = [-2 delta -2];
%     grad(range(3:end-2)) = grad(range(3:end-2)) -2*els2(i)*(2*fds(2:end-3).*(xchunk(1:end-4)-xchunk(2:end-3))./...
%         delta(2:end-3).^2 + 2*fds(3:end-2).*(xchunk(4:end-1)-xchunk(2:end-3))./delta(3:end-2).^2 + ...
%         2*fds(4:end-1).*(xchunk(4:end-1)-xchunk(5:end))./delta(4:end-1).^2);
%      
%     grad(range(2)) = grad(range(2)) - 2*els2(i) * (2*fds(2)*(xchunk(3)-xchunk(1))/delta(2)^2 + 2*fds(3)*(xchunk(3)-xchunk(4))/delta(3)^2);
%     grad(range(end-1)) = grad(range(end-1)) - 2*els2(i) * (2*fds(end-2)*(xchunk(end-3)-xchunk(end-2))/delta(end-2)^2 + 2*fds(end-1)*(xchunk(end)-xchunk(end-2))/delta(end-1)^2);
%     grad(range(1)) = grad(range(1)) - 2*els2(i) * (2*fds(2)*(xchunk(2)-xchunk(3))/delta(2)^2);
%     grad(range(end)) = grad(range(end)) - 2*els2(i) * (2*fds(end-1)*(xchunk(end-3)-xchunk(end-2))/delta(end-1)^2);
%     grad2(range(3:end-2)) = grad2(range(3:end-2)) - 2*els2(i) * ()


end

if ~isempty(llist)
    
  a = llist(2:end); %list of targets constraints
  %store weather it's a pos or neg constraint
  mults = ones(size(a)); 
  mults(a < 0) = -1;
  a = abs(a);
  
  ll = ll - llist(1) * sum(abs(k(a)));
  grad(a) = grad(a) - mults*llist(1); 
  
  % assume all terms in llist are either greater or less than zero to their current values
%   grad(a) = grad(a) - (2*(k(a) >= 0)-1) * llist(1);
  % grad(a) = grad(a) - constr*llist(1);
end

if ~isempty(hold_const)
  grad(hold_const(1,:)) = -2*HCmult*hcgrad;
end

% Flip signs of ll and grad for minimization
ls = length(ks); %normalize by number of spikes
if ls > 0
  ll=-ll/ls;
  grad=-grad'/ls; %'
end
