function [T,table]=find_best_oblique_rotation(T,loading_mat)

global A;
A = loading_mat;

al=1;
table=[];
for iter=0:30
    f=localization_2drf_entropy(T,loading_mat);
    G=localization_2drf_entropy_deriv(T,loading_mat);
%     f=ff(T);
%     G=Gf(T);
    Gp=G-T*diag(sum(T.*G));
    s=norm(Gp,'fro');
    table=[table;iter f log10(s) al];
    if s<10^(-5),break,end
    al=2*al;
    for i=0:10
        X=T-al*Gp;
        v=1./sqrt(sum(X.^2));
        Tt=X*diag(v);
        ft=localization_2drf_entropy(Tt,loading_mat);
%         ft=ff(Tt);
        if ft<f-.5*s^2*al,break,end
        al=al/2;
    end
    T=Tt;
    iter
end