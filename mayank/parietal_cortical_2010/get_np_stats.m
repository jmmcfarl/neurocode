function [] = get_np_stats(vec1,vec2)

med1 = median(vec1);
med2 = median(vec2);
lq1 = prctile(vec1,25);
lq2 = prctile(vec2,25);
uq1 = prctile(vec1,75);
uq2 = prctile(vec2,75);

p = ranksum(vec1,vec2);

fprintf('vector 1.  median: %.2e  lqr: %.2e  uqr: %.2e\n',med1,lq1,uq1);
fprintf('vector 2.  median: %.2e  lqr: %.2e  uqr: %.2e\n',med2,lq2,uq2);
fprintf('p value  %.2e\n',p);
if length(vec1) == length(vec2)
    p2 = signrank(vec1,vec2);
    fprintf('paired p value  %.2e\n',p2)
end
