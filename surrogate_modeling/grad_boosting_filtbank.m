function [f,grad] = grad_boosting_filtbank(b_vec,y_tilda,filt_outs)

g_int = filt_outs*b_vec;
g = g_int;
g(g < 0) = 0;
f = sum((y_tilda - g).^2);

g_prime = zeros(size(g));
g_prime(g_int > 0) = 1;
grad = -2*((y_tilda-g).*g_prime)'*filt_outs;
grad = grad(:);