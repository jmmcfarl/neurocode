function [C,p] = nancorr(data,data2)

if nargin == 1
udata = find(~any(isnan(data),2));
[C,p] = corr(data(udata,:));
elseif nargin == 2
 udata = find(~any(isnan(data),2) & ~any(isnan(data2),2));
[C,p] = corr(data(udata,:),data2(udata,:));
   
else
    error('Too many inputs');
end