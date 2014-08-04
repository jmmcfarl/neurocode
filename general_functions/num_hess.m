function [H, NFV] = num_hess(func, X , h)

H = zeros(length(X));

% for each dimension of objective function
for i=1:length(X)
    % derivative at first point (left)
    x1 = X;
    x1(i) = X(i) - h;
    [df1, NFV] = num_grad(func, x1, NFV, h);
    
    % derivative at second point (right)
    x2 = X;
    x2(i) = X(i) + h;
    [df2, NFV] = num_grad(func, x2, NFV, h);
    
    % differentiate between the two derivatives
    d2f = (df2-df1) / (2*h);
    
    % assign as row i of Hessian
    H(i,:) = d2f';
end