function [vecs, vals, GaussParams] = compiSTAC(mu1,A1,mu0,A0,ndims);
% [vecs, vals, GaussParams] = compiSTAC(mu1,A1,mu0,A0,ndims);
%
% Computes a basis for the subspace which captures the maximal KL
% divergence between the two Gaussians N(mu1,A1) and N(mu0,A0).
%
% Implicitly, N(mu1,A1) is a Gaussian description of the spike-triggered
% ensemble, while N(mu0,A0) is a description of the raw stim distribution
%
% Whitens with respect to N(mu0,A0), so computation is simplified.
%
% Inputs:  
%   mu1 = mean of 1st Gaussian (column vector)
%   A1 = covariance of 1st Gaussian 
%   mu2, A2 = mean, cov of 2nd Gaussian
%   ndims =  number of dimensions to peel off
%
% Ouptuts: 
%   vecs - (n x ndims) matrix, whose columns form an (ordered) basis for the 
%          maximal information-preserving subspace of degree ndims
%   vals - value of the KL divergence as subspace dimensionality increases
%          from 1 to ndims 
%   GaussParams - structure containing the means and covariances of the two
%          Gaussians projected into the subspace of interest.  
%          (Useful if we wish to use ratio-of-Gaussians to describe the
%          nonlinearity).
%
% Last updated: J Pillow 06/2006

% Initialize some optimization params
vecs = [];
n = length(mu1);
UB = ones(n,1);
LB = -ones(n,1);
opts = optimset('display', 'off', 'gradobj', 'off', 'largescale', 'off', ...
    'maxfunevals', 200000, 'maxiter', 50);

% Whiten, Convert to zero-mean coordinates
A0whiten = mysqrtm(inv(A0));
mu = A0whiten*(mu1-mu0);
A = A0whiten*A1*A0whiten;

% Compute SVD of covariance, for initial guesses
[u,s,v] = svd(A);
a = min(ndims,floor(n/2));
k0s = [u(:,[1:a end-a+1:end]) mu./norm(mu)];

bv = [];
vA = [];
vAv = [];

j = 1;
while j <= min(ndims,n-1)
    BackingUP = 0;
    %  Start by finding best starting point for optimization
    kstrt = orthogonalsubset(vecs, k0s);
    v0s = 0;
    for ii = 1:size(kstrt,2);
        v0s(ii) = negKLsubspace(kstrt(:,ii),mu, A, bv, vA, vAv,vecs);
    end
    imin = find(v0s == min(v0s));  imin = imin(1);
    k0 = kstrt(:,imin);
    
    % Perform optimization -- restart if optimization doesn't terminate
    Beq = zeros(j-1,1);
    [k,v0,exitflag] = fmincon(@negKLsubspace, k0,[],[],vecs',Beq,LB,UB,...
        @NormEq1,opts,mu,A,bv,vA,vAv,vecs);
    if exitflag<1  % Check convergence
        %fprintf(1, 'iSTAC-- possible error: optimization not terminated; j=%d\n',j);
        % Note: Up the optimization parameter 'niter' if worried about
        % convergence
    end
    if j > 1  % normalize k with respect to previous vecs
        k = k-vecs*(vecs'*k);
        k = k./norm(k);
    end

    % Compte KL divergence along this dimension
    vecs(:,j) = k;
    
    vals(j,1) = compDklgaussian(vecs'*mu, vecs'*A*vecs, zeros(j,1), eye(j));
    valdiffs = [vals(1); diff(vals)];
    valmarginals(j,1) = compDklgaussian(k'*mu, k'*A*k, 0, 1);
    
    % Check that vals is smaller than all previous values
    if BackingUP >= 3
        BackingUP = 0;
    elseif (valdiffs(j) > min(valdiffs(1:j-1))) & (j < n/2)
        jj = find(valdiffs(1:j-1) < valdiffs(j));
        k0s = [k k0s];
        vecsorig = vecs;
        vecs = vecs(:,1:jj(1)-1);
        vals = vals(1:jj(1)-1);
        j = jj(1);
        fprintf(1, 'Going back to iter #%d (valdiff=%.4f)\n', j,valdiffs(end));
        BackingUP = 1;
        %
    elseif j>1 
        vv = vecs(:,[1:j-2 j]);
        valtst = compDklgaussian(vv'*mu, vv'*A*vv, zeros(j-1,1), eye(j-1));
        if valtst > vals(j-1)
            fprintf(1, 'Wrong dim possibly stripped off [%.4f %.4f]; going back to prev dim\n', ...
                vals(j-1), valtst);
            k0s = [k k0s];
            vecs = vecs(:,1:j-2);
            vals = vals(1:j-2);
            j = j-1;
            BackingUP = BackingUP+1;
        end
    end
    if ~BackingUP
        fprintf(1,' Stripped off dimension %d, KL div=[%2.4f %2.4f]\n', ...
            j, valdiffs(j), valmarginals(j));
        j = j+1;
    end

    % compute projection of A and mu onto vecs
    bv = vecs'*mu;
    vA = vecs'*A;
    vAv = vecs'*A*vecs;
end

vecs = A0whiten*vecs;
vecs = gsorth(vecs);

GaussParams.mu1 = vecs'*mu1;
GaussParams.mu0 = vecs'*mu0;
GaussParams.v1 = vecs'*A1*vecs;
GaussParams.v2 = vecs'*A0*vecs;

%  -------------------------------
function vorth = orthogonalsubset(B, vecs)
%  orthogonalize set of vectors with respect to columns of B ;
%  remove those in the column space of B

etol = 1e-10;

if isempty(B)
    vorth = vecs;
    return;
end

[n,m] = size(B);
Binv = inv(B'*B);
vorth = [];
nv = 0;
for j = 1:size(vecs,2);
    k = vecs(:,j) - B*(Binv*B'*vecs(:,j));
    if norm(k) > etol;
        nv = nv+1;
        vorth(:,nv) = k./norm(k);
    end
end

%  -------------------------------
function [c,ceq,dc,dceq] = NormEq1(x, varargin);
% [c,ceq] = NormEq1(x, varargin);
%
% nonlinear function for implementing the constraint norm(x) = 1;

c = [];
ceq = x'*x-1;

if nargout>2
    dc = [];
    dceq = 2*x;
end
