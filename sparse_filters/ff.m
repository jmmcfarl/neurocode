
function f = ff(T)

global A;

[p,k]=size(A);
L=A*inv(T');
L2=L.^2;
N=ones(k,k)-eye(k);
f=sum(sum(L2.*(L2*N)));