function yhat = my_sigmoid(param,x)
a = param(1);
b = param(2);
c = param(3);
d = param(4);
yhat = a*(1.0./(1.0+exp((b-x)/c)))+d;
