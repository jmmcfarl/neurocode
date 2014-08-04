function G = localization_2drf_entropy_deriv(T,loading_mat)

k=size(T,1);
ep=.0001;
Z=zeros(k,k);
G=Z;
for r=1:k
    for s=1:k
        dT=Z;
        dT(r,s)=ep;
        G(r,s)=(localization_2drf_entropy(T+dT,loading_mat)-localization_2drf_entropy(T-dT,loading_mat))/(2*ep);
    end
end