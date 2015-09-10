function L = logsumexp(x,dim)
% L = logsumexp(x,dim)
% compute the log(sum(exp(x))), summing across dimension dim

if nargin < 2
    dim = 1;
end

sz = size(x);

if sum(sz ~= 1) > 1
    m = max(x,[],dim);
    resh_dims = ones(1,length(sz));
    resh_dims(dim) = sz(dim);
    rm = repmat(m,resh_dims);
else
    x = x(:);
    m = max(x);
    rm = m;
end
L = log(sum(exp(x-rm),dim)) + m;