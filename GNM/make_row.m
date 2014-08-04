function [nv,L] = make_row( v )
% 
% Usage: [nv,L] = make_row( v )
%
% Makes row vector out of ambiguously formatted vector

[a b] = size(v);
if a > b
  nv = v';  %'
else
  nv = v;
end

L = length(nv);
