function [gabor_params,LL] = fit_gabor_energy_mod_analytic(XX,YY,init_params,Xmat,Robs,hold_const,LB,UB)

options = optimset('GradObj','on');
options.Display = 'iter';
% options.Algorithm = 'sqp';

K0 = init_params(hold_const == 0);
fixed_params = init_params(hold_const == 1);

[params LL] = fmincon( @(K) get_pgabor_LL_analytic(K, Xmat, Robs,XX,YY,hold_const,fixed_params), ...
    K0,[],[],[],[],LB(hold_const==0),UB(hold_const==0),[],options);
gabor_params = init_params;
gabor_params(hold_const == 0) = params;

end

function [LL,grad] = get_pgabor_LL_analytic(K,Xmat,Robs,XX,YY,hold_const,fixed_params)

Fsd = XX(1,2)-XX(1,1);

params = zeros(8,1);
params(hold_const==1) = fixed_params;
params(hold_const==0) = K;

x0 = params(1);
y0 = params(2);
theta = params(3);
lambda = params(4);
sigma = params(5);
gamma = params(6);
ew = params(7);
c = params(8);

%compute useful quantities
xp = (XX-x0)*cos(theta) + (YY - y0)*sin(theta);
yp = -(XX-x0)*sin(theta) + (YY - y0)*cos(theta);
genv = exp(-xp.^2/2/sigma^2 - gamma*yp.^2/2/sigma^2);
cfun = cos(2*pi*xp/lambda);
sfun = cos(2*pi*xp/lambda + pi/2);
kappa = 2*pi*sigma/lambda;
norm_const = sqrt(2*sqrt(gamma))*Fsd/sqrt(pi*sigma^2);
dc_offset = sqrt(2)*exp(-kappa^2);

cur_mask1 =  norm_const*(cfun-dc_offset).*genv;
cur_mask2 = norm_const*(sfun-dc_offset).*genv;

mask1_out = Xmat*cur_mask1(:);
mask2_out = Xmat*cur_mask2(:);

energy_out = ew*sqrt(mask1_out.^2+mask2_out.^2);

g = energy_out + c;
too_large = find(g > 100);
r = log(1+exp(g));
r(too_large) = g(too_large);

r(r < 1e-20) = 1e-20;

LL = sum(Robs.*log(r)-r);
Nspks = sum(Robs);
LL = -LL/Nspks;

residual = (Robs./r - 1).*exp(g)./(1 + exp(g));
residual(too_large) = (Robs(too_large)./r(too_large)-1);

residual_temp2 = residual*ew.*(mask1_out.^2+mask2_out.^2).^(-1/2);

grad = zeros(size(params));
grad(8) = sum(residual);
grad(7) = sum(residual.*sqrt(mask1_out.^2+mask2_out.^2));

if hold_const(1) == 0
    %calculate df/dx0
    df0_dx0 = norm_const*(cfun-dc_offset).*genv.*(xp/sigma^2*cos(theta) - gamma*yp/sigma^2*sin(theta)) + ...
        genv.*sfun.*2*pi/lambda*cos(theta);
    df1_dx0 = norm_const*(sfun-dc_offset).*genv.*(xp/sigma^2*cos(theta) - gamma*yp/sigma^2*sin(theta)) - ...
        genv.*cfun.*2*pi/lambda*cos(theta);
    df0_dx0_out = Xmat*df0_dx0(:);
    df1_dx0_out = Xmat*df1_dx0(:);
    grad(1) = sum(residual_temp2.*(mask1_out.*df0_dx0_out + mask2_out.*df1_dx0_out));
end

if hold_const(2) == 0
    df0_dy0 = norm_const*(cfun-dc_offset).*genv.*(xp/sigma^2*sin(theta) + gamma*yp/sigma^2*cos(theta)) + ...
        genv.*sfun.*2*pi/lambda*sin(theta);
    df1_dy0 = norm_const*(sfun-dc_offset).*genv.*(xp/sigma^2*sin(theta) + gamma*yp/sigma^2*cos(theta)) - ...
        genv.*cfun.*2*pi/lambda*sin(theta);
    df0_dy0_out = Xmat*df0_dy0(:);
    df1_dy0_out = Xmat*df1_dy0(:);
    grad(2) = sum(residual_temp2.*(mask1_out.*df0_dy0_out + mask2_out.*df1_dy0_out));
end

if hold_const(4) == 0
    df0_dlambda = -norm_const*genv.*sfun.*(2*pi*xp/lambda^2);
    df1_dlambda = norm_const*genv.*cfun.*(2*pi*xp/lambda^2);
    df0_dlambda_out = Xmat*df0_dlambda(:);
    df1_dlambda_out = Xmat*df1_dlambda(:);
    grad(4) = sum(residual_temp2.*(mask1_out.*df0_dlambda_out + mask2_out.*df1_dlambda_out));
end
if hold_const(5) == 0
   df0_dsigma = norm_const*cfun.*genv.*(xp.^2/sigma^3 + gamma*yp.^2/sigma^3);
   df1_dsigma = norm_const*sfun.*genv.*(xp.^2/sigma^3 + gamma*yp.^2/sigma^3);
   df0_dsigma_out = Xmat*df0_dsigma(:);
   df1_dsigma_out = Xmat*df1_dsigma(:);
    grad(5) = sum(residual_temp2.*(mask1_out.*df0_dsigma_out + mask2_out.*df1_dsigma_out));
end

grad = -grad/Nspks;
grad = grad(hold_const==0);

end

