function [gabor_params,LL] = fit_gabor_params_v2(Xmat,Robs,msize,const_params,LB,UB)

options.Display = 'iter';
% options.Algorithm = 'sqp';
% options.MaxFunEvals = 1000;

K0 = const_params(4);
fixed_params = const_params;

% [params LL exitflag] = fminunc( @(K) get_pgabor_LL(K, Xmat, Robs, msize,hold_const,fixed_params), K0,options);
[params LL] = fmincon( @(K) get_pgabor_LL(K, Xmat, Robs, msize,fixed_params), ...
    K0,[],[],[],[],LB,UB,[],options);
gabor_params = init_params;
gabor_params(hold_const == 0) = params;

end

function LL = get_pgabor_LL(K,Xmat,Robs,msize,fixed_params)

params = fixed_params;
params(4) = K;
phase = atan2(params(9),params(8));

cur_mask1 = get_gabor_mask(params(1:4),phase,msize);

mask1_out = Xmat*cur_mask1(:);
lin_out = params(8)*mask1_out;

g = lin_out + params(10);
too_large = find(g > 100);
expg = exp(g);
r = log(1+expg);
r(too_large) = g(too_large);

r(r < 1e-20) = 1e-20;

LL = sum(Robs.*log(r)-r);
Nspks = sum(Robs);
LL = -LL/Nspks;

residual = (Robs./r + 1).*(expg./(1+expg));
residual(too_large) = (Robs(too_large)./r(too_large)+1);

deriv_mask = get_gabor_mask_deriv(params(1:4),phase,msize);
deriv_maskout = Xmat*deriv_mask(:);
LLgrad = sum(residual.*deriv_maskout);
LLgrad = -LLgrad/Nspks;


end

function cur_mask = get_gabor_mask(params,phase,msize)

x0 = params(1);
y0 = params(2);
theta = params(3);
lambda = params(4);

xax = -floor(msize(2)/2):floor(msize(2)/2);
yax = -floor(msize(1)/2):floor(msize(1)/2);
[X,Y] = meshgrid(xax,yax);

xp = (X-x0)*cos(theta)+(Y-y0)*sin(theta);
yp = -(X-x0)*sin(theta)+(Y-y0)*cos(theta);

cur_mask = exp(-((xp.^2 + yp.^2)*2/lambda^2)) .* cos(2*pi*xp/lambda+phase);

end

function cur_mask = get_gabor_mask_deriv(params,phase,msize)

x0 = params(1);
y0 = params(2);
theta = params(3);
lambda = params(4);

xax = -floor(msize(2)/2):floor(msize(2)/2);
yax = -floor(msize(1)/2):floor(msize(1)/2);
[X,Y] = meshgrid(xax,yax);

xp = (X-x0)*cos(theta)+(Y-y0)*sin(theta);
yp = -(X-x0)*sin(theta)+(Y-y0)*cos(theta);

cur_mask = exp(-((xp.^2 + yp.^2)*2/lambda^2)) .* cos(2*pi*xp/lambda+phase);

cur_mask = 2*cur_mask.*(X.^2+Y.^2)/lambda^3;
cur_mask = cur_mask + 2*pi*xp/lambda^2*exp(-((xp.^2 + yp.^2)*2/lambda^2)) .* sin(2*pi*xp/lambda+phase);

end
