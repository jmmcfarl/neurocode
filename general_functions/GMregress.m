function [slope,intercept] = GMregress(X,Y)
%[slope,intercept] = gmregress(X,Y)
%compute the geometric-mean regression of Y on X

%convert to column vector format
if size(X,2) > size(X,1)
    X = X';
end
if size(Y,2) > size(Y,1)
    Y = Y';
end

slope_sign = sign(corr(X,Y));
SSX = sum((X - mean(X)).^2);
SSY = sum((Y - mean(Y)).^2);
slope = slope_sign*sqrt(SSY/SSX);
intercept = mean(Y) - mean(X)*slope;
