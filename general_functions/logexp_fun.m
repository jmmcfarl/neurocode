function lexp_out = logexp_fun(x,params)

internal = params(2)*(x + params(1));
too_big = find(internal > 50);
lexp = log(1+exp(internal));
lexp(too_big) = params(2)*(x(too_big)+params(1));

lexp_out = params(3)*lexp;

if length(params) == 4
    lexp_out = lexp_out + params(4);
end