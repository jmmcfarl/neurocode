function pdf = inv_gauss_pdf(x,mu,lambda)

pdf = sqrt(lambda./(2*pi*x.^3)).*exp(-lambda*(x-mu).^2./(2*mu^2*x));

