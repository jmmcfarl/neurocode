function cur_mask = get_pgauss_mask_v2(X,Y,params)

x0 = params(1);
y0 = params(2);
theta = params(3);
lambda = params(4);
sigma = params(5);
gamma = params(6);

xp = (X-x0)*cos(theta)+(Y-y0)*sin(theta);
yp = -(X-x0)*sin(theta)+(Y-y0)*cos(theta);

env = exp(-(xp.^2/2/sigma^2 + gamma*yp.^2/2/sigma^2));
cur_mask = env/sum(env(:));

end