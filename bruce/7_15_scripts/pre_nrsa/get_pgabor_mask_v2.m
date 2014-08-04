function cur_mask = get_pgabor_mask_v2(X,Y,params,phase)

x0 = params(1);
y0 = params(2);
theta = params(3)-pi/2;
lambda = params(4);
sigma = params(5);
gamma = params(6);

xp = (X-x0)*cos(theta)+(Y-y0)*sin(theta);
yp = -(X-x0)*sin(theta)+(Y-y0)*cos(theta);

env = exp(-(xp.^2/2/sigma^2 + gamma*yp.^2/2/sigma^2));
env = env/sum(env(:));
osc =cos(2*pi*xp/lambda+phase);
dc = sum(env(:).*osc(:));
osc =cos(2*pi*xp/lambda+phase)-dc;
cur_mask = env.*osc;

cur_mask = cur_mask/norm(cur_mask(:));
% cur_mask = cur_mask - mean(cur_mask(:)); %should make this DC subtraction localized to the gabor

end