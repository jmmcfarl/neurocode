function d = compDklgaussian(mu1,C1,mu2,C2)
%  d = compDklgaussian(mu1, C1, mu2, C2)
%
%  Computes the KL divergence between two multi-dimensional Gaussians
%
%  Inputs:
%    mu1 = mean of 1st Gaussian
%    C1 = covariance of 1st Gaussian
%    mu2,C2 = mean,cov of 2nd Gaussian
%
%  Notes:
%     D_KL = Integral(p1 log (p1/p2))
%  Analytically:  (where |*| = Determinant, Tr = trace, n = dimensionality
%     =  1/2 log(|C2|/|C1|) + 1/2 Tr(C1^.5*C2^(-1)*C1^.5)
%        + 1/2 (mu2-mu1)^T*C2^(-1)*(mu2-mu1) - 1/2 n

if nargin == 1;
    DD = mu1;
    mu1 = DD.mu1;
    mu2 = DD.mu2;
    C1 = DD.v1;
    C2 = DD.v2;
end

n = length(mu1);
b = mu2-mu1;
C2inv = inv(C2);
C1sqrt = mysqrtm(C1);

Term1 = C1sqrt*C2inv*C1sqrt;
Term2 = b'*C2inv*b;


det1 = det(C1);
det2 = det(C2);
tol = 1e8;
if (cond(C1)>tol) | (cond(C2)>tol)
    fprintf('Determinants not non-zero: %.3f %.3f\n', det1, det2);
    fprintf('Ignoring determinants...\n');
    Term3 = 0;
else
    Term3 = .5*log(det2/det1);
end

d = .5*trace(Term1) + .5*Term2 - .5*n + Term3;


