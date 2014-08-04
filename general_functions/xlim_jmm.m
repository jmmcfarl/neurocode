function [] = xlim_jmm(xl)

if any(isnan(xl))
    xl = [0 1];
end
if range(xl) == 0
    xl = [0 1];
end
xlim(xl);