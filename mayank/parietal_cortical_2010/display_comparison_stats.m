function [] = display_comparison_stats(data1,data2)

m1 = median(data1);
m2 = median(data2);
lqr1 = prctile(data1,25);
lqr2 = prctile(data2,25);
uqr1 = prctile(data1,75);
uqr2 = prctile(data2,75);

p = signrank(data1,data2);
p1 = signrank(data1);
p2 = signrank(data2);

fprintf('data 1 median %.3e  iqr  %.3e - %.3e  p = %.3e\n',m1,lqr1,uqr1,p1)
fprintf('data 2 median %.3e  iqr  %.3e - %.3e  p = %.3e\n',m2,lqr2,uqr2,p2)
fprintf('comparison  %.3e\n',p)