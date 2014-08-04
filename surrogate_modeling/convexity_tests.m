
%%
syms x R k1 k2

%logexp spk NL and quadratic internal
% y = -(R*log(log(1+exp(k1*x+k2^2*x^2))) - log(1+exp(k1*x+k2^2*x^2)));

% exp spk NL and quadratic internal
y = -(R*(k1*x+k2^2*x^2) - exp(k1*x+k2^2*x^2));

%exp spk NL and exp + lin internal
% y = -(R*(k1*x+exp(k2*x)) - exp(k1*x+exp(k2*x)));

%logexp spk NL and exp+lin internal
% y = -(R*(log(log(1+exp(k1*x+exp(k2*x))))) - log(1+exp(k1*x+exp(k2*x))));

%exp spk NL and 2 exp internals
% y = -(R*(exp(k1*x)+exp(k2*x)) - exp(exp(k1*x)+exp(k2*x)));

%logexp spk NL and 2 exp internals
% y = -(R*(log(log(1+exp(exp(k1*x)+exp(k2*x))))) - log(1+exp(exp(k1*x)+exp(k2*x))));

%exp spk NL and 2 logexp internals
% y = -(R*(log(1+exp(k1*x))+log(1+exp(k2*x))) - exp(log(1+exp(k1*x))+log(1+exp(k2*x))));

%logexp spk NL and 2 logexp internals
y = -(R*(log(log(1+exp(log(1+exp(k1*x))+log(1+exp(k2*x)))))) - log(1+exp(log(1+exp(k1*x))+log(1+exp(k2*x)))));

dydk1 = simplify(diff(y,k1));
dy2dk12 = simplify(diff(dydk1,k1));
dydk2 = simplify(diff(y,k2));
dy2dk22 = simplify(diff(dydk2,k2));
dy2dk1k2 = simplify(diff(dydk1,k2));
dy2dk2k1 = simplify(diff(dydk2,k1));

%%
f = matlabFunction(dy2dk12)
% f = matlabFunction(dy2dk22)
% f = matlabFunction(dy2dk1k2)

fdet = matlabFunction(dy2dk12*dy2dk22 - dy2dk1k2*dy2dk2k1)

%%
close all
% xvals = linspace(-5,5,500);
xv = 3;
Rv = 0.5;
kv = 1;
k2v = 1;
kvals = linspace(-5,5,500);

yv = f(Rv,kv,kvals,xv);
% yv = f(kv,k2v,xvals);
yvd = fdet(Rv,kv,kvals,xv);

figure
plot(kvals,yvd)
set(gca,'yscale','log')
% xlim([-3 3])

figure
plot(kvals,yv)
set(gca,'yscale','log')
% xlim([-3 3])

%%
syms x R k1 k2 b1 b2
y = -(R*(log(log(1+exp(log(1+exp(b1*k1*x))+log(1+exp(b2*k2*x)))))) - log(1+exp(log(1+exp(b1*k1*x))+log(1+exp(b2*k2*x)))));

dydb1 = simplify(diff(y,b1));
dy2db12 = simplify(diff(dydb1,b1));
dydb2 = simplify(diff(y,b2));
dy2db22 = simplify(diff(dydb2,b2));
dy2db1b2 = simplify(diff(dydb1,b2));
dy2db2b1 = simplify(diff(dydb2,b1));

f = matlabFunction(dy2db12)
fdet = matlabFunction(dy2db12*dy2db22 - dy2db1b2*dy2db2b1)

%%
xvals = linspace(-5,5,500);
Rv = 0.5;
kv = 1;
k2v = 1;
b1v = 1;
b2v = 1;

yv = f(Rv,b1v,b2v,kv,k2v,xvals);
yvd = fdet(Rv,b1v,b2v,kv,k2v,xvals);

figure
plot(xvals,yvd)
set(gca,'yscale','log')
xlim([-3 3])

figure
plot(xvals,yv)
set(gca,'yscale','log')
xlim([-3 3])

%%
syms x R k b t
y = -(R*(log(log(1+exp(log(1+exp(b*k*(x-t))))))) - log(1+exp(log(1+exp(b*k*(x-t))))));

dydb = simplify(diff(y,b));
dy2db2 = simplify(diff(dydb,b));
dydt = simplify(diff(y,t));
dy2dt2 = simplify(diff(dydt,t));
dy2dbt = simplify(diff(dydb,t));
dy2dtb = simplify(diff(dydt,b));

f = matlabFunction(dy2db2)
fdet = matlabFunction(dy2db2*dy2dt2 - dy2dbt*dy2dtb)

%%
xvals = linspace(-5,5,500);
Rv = 0.5;
kv = 1;
bv = 1;
tv = 0;

yv = f(Rv,bv,kv,tv,xvals);
yvd = fdet(Rv,bv,kv,tv,xvals);

figure
plot(xvals,yvd)
set(gca,'yscale','log')
xlim([-3 3])

figure
plot(xvals,yv)
set(gca,'yscale','log')
xlim([-3 3])

%%
syms x1 x2 R k1 k2
y = -(R*(log(log(1+exp(log(1+exp(k1*x1+k2*x2)))))) - log(1+exp(log(1+exp(k1*x1+k2*x2)))));

dydk1 = simplify(diff(y,k1));
dy2dk12 = simplify(diff(dydk1,k1));
dydk2 = simplify(diff(y,k2));
dy2dk22 = simplify(diff(dydk2,k2));
dy2dk1k2 = simplify(diff(dydk1,k2));
dy2dk2k1 = simplify(diff(dydk2,k1));

f = matlabFunction(dy2dk12)
fdet = matlabFunction(dy2dk12*dy2dk22 - dy2dk1k2*dy2dk2k1)

%%
[x1vals,x2vals] = meshgrid(linspace(-3,3,100),linspace(-3,3,100));
Rv = 0.5;
k1v = 1;
k2v = 1;

yv = f(Rv,k1v,k2v,x1vals,x2vals);
% yvd = fdet(Rv,b1v,b2v,kv,k2v,xvals);

figure
plot(xvals,yvd)
set(gca,'yscale','log')
xlim([-3 3])

figure
plot(xvals,yv)
set(gca,'yscale','log')
xlim([-3 3])

%%
syms x R k c
y = -(R*(log(log(1+exp(log(1+exp(k*x))))+c)) - log(1+exp(log(1+exp(k*x))))+c);

dydk = simplify(diff(y,k));
dy2dk2 = simplify(diff(dydk,k));
dydc = simplify(diff(y,c));
dy2dc2 = simplify(diff(dydc,c));
dy2dkc = simplify(diff(dydk,c));
dy2dck = simplify(diff(dydc,k));

f = matlabFunction(dy2dc2)
fdet = matlabFunction(dy2dc2*dy2dk2 - dy2dck*dy2dkc)

%%
xvals = linspace(-20,20,500);
Rv = 0.5;
kv = 1;
cv = -1;

yv = f(Rv,cv,kv,xvals);
yvd = fdet(Rv,cv,kv,xvals);

figure
plot(xvals,yvd)
set(gca,'yscale','log')
xlim([-3 3])

figure
plot(xvals,yv)
set(gca,'yscale','log')
xlim([-3 3])
