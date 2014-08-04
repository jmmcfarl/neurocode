function d2 = dist_shift2d( d1, sh, dim, wrap )
%
% Usage: d2 = dist_shift( d1, sh, dim, <wrap> )
%
% shifts 2-dimensional distribution in dimension specified by dim
%
% Default (wrap = 1) wraps around, otherwiseputs zeros at the end
% negative shifts distribution left

sh = floor(sh+0.5);
if nargin < 4
  wrap = 1;
end

[L(1) L(2)] = size(d1);

if sh > L(dim)
  sh = L(dim);
end
if sh < -L(dim)
  sh = -L(dim);
end

mapper = mod((0:(L(dim)-1))-sh,L(dim))+1;
if dim == 1
  d2 = d1(mapper,:);
else
  d2 = d1(:,mapper);
end

if wrap == 0
  if dim == 1
    if sh > 0
      d2(1:sh,:) = 0;
    else
      d2(L(1)-(0:(-sh-1)),:) = 0;
    end
  else
    if sh > 0
      d2(:,1:sh) = 0;
    else
      d2(:,L(2)-(0:(-sh-1))) = 0;
    end
  end
end
