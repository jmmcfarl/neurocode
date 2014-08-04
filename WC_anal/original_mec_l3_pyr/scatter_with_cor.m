function [corr,p_val] = scatter_with_cor(data1,data2)

if size(data1,1) ~= size(data2,1)
    data1 = data1';
end

good_pts = ~isnan(data1) & ~isnan(data2);

plot(data1(good_pts),data2(good_pts),'o')

[a,b] = corrcoef(data1(good_pts),data2(good_pts));
title(sprintf('C = %0.2g  p = %0.2g',a(2,1),b(2,1)))

corr = a(2,1);
p_val = b(2,1);