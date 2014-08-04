function cur_mask = get_pgabor_mask(params,phase,msize)

x0 = params(1);
y0 = params(2);
theta = params(3);
lambda = params(4);
sigma = params(5);
gamma = params(6);

xax = -floor(msize(2)/2):floor(msize(2)/2);
yax = -floor(msize(1)/2):floor(msize(1)/2);
xax(msize(2)+1:end) = [];
yax(msize(1)+1:end) = [];
[X,Y] = meshgrid(xax,yax);

xp = (X-x0)*cos(theta)+(Y-y0)*sin(theta);
yp = -(X-x0)*sin(theta)+(Y-y0)*cos(theta);

sigmax = sigma;
sigmay = sigma;
cur_mask = exp(-(xp.^2/2/sigmax^2 + gamma*yp.^2/2/sigmay^2)) .* cos(2*pi*xp/lambda+phase);
cur_mask = cur_mask - mean(cur_mask(:));

end