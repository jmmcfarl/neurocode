
function G=Gf(T)

k=size(T,1);
ep=.0001;
Z=zeros(k,k);
G=Z;
for r=1:k
    for s=1:k
        dT=Z;
        dT(r,s)=ep;
        G(r,s)=(ff(T+dT)-ff(T-dT))/(2*ep);
    end
end