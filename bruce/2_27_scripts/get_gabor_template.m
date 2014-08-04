function gabor = get_gabor_template(X,Y,x0,y0,theta,lambda,psi,bandwidth,ar)

theta = theta - pi/2;
xp = (X-x0)*cos(theta)+(Y-y0)*sin(theta);
yp = -(X-x0)*sin(theta)+(Y-y0)*cos(theta);

sigmax = bandwidth*lambda;
sigmay = sigmax*ar;
gabor = exp(-(xp.^2/2/sigmax^2 + yp.^2/2/sigmay^2)) .* cos(2*pi*xp/lambda+psi);
gabor = gabor - mean(gabor(:));