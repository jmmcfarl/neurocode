function [C,p] = nancorr(data,data2,ctype)

if nargin < 2
    data2 = [];
end
if nargin < 3
    ctype = 'pearson';
end

if isempty(data2)
    udata = find(~any(isnan(data),2));
    [C,p] = corr(data(udata,:),'type',ctype);
else
    udata = find(~any(isnan(data),2) & ~any(isnan(data2),2));
    [C,p] = corr(data(udata,:),data2(udata,:),'type',ctype);
end