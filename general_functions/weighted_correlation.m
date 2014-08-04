function [weighted_r,raw_r] = weighted_correlation(X,Y,weights)

X = X(:);
Y = Y(:);
weights = weights(:);
weights = weights/sum(weights);

mx = sum(weights.*X);
my = sum(weights.*Y);
cov_xy = sum(weights.*(X - mx).*(Y - my));
cov_xx = sum(weights.*(X - mx).*(X - mx));
cov_yy = sum(weights.*(Y - my).*(Y - my));

weighted_r = cov_xy./sqrt(cov_xx*cov_yy);
raw_r = corr(X,Y);

