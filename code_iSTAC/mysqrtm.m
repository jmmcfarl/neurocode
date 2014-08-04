function Msqrt = mysqrtm(M);
% Msqrt = mysqrtm(M);
% 
% Function for computing the sqrt of symmetric matrices
% Handles ill-conditioned matrices gracefully
% Note: M must be symmetric!

[u,s,v] = svd(M);
Msqrt = u*sqrt(s)*u';
