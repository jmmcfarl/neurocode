function g = my_gaussfun(beta,X)

min_var = 0.02^2;
if beta(4) < min_var
    beta(4) = min_var;
    pen = 100;
else
    pen = 0;
end
g = beta(1) + beta(2)*exp(-(X-beta(3)).^2/(2*beta(4))) + pen;