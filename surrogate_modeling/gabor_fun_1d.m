function y = gabor_fun_1d(x,x0,sigma,lambda,psi,a,b)

y = b.*exp(-((x-x0).^2./2./sigma.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi))+a;
