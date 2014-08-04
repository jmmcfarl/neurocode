function [gabor_params,LL] = fit_temporal_gabor_params(init_params,Xmat,gain,offset,Robs,stim_ids,msize,hold_const,LB,UB)

options.Display = 'iter';
% options.Algorithm = 'interior-point';
% options.Algorithm = 'spq';
% options.Hessian = 'lbfgs';
options.MaxFunEvals = 1000;

K0 = init_params(hold_const == 0);
fixed_params = init_params(hold_const == 1);
K0 = K0;
[params LL] = fmincon( @(K) get_pgabor_LL(K, Xmat, gain,offset, Robs, stim_ids, msize,hold_const,fixed_params), ...
    K0,[],[],[],[],LB(hold_const==0),UB(hold_const==0),[],options);
gabor_params = init_params;
gabor_params(hold_const == 0) = params;

end

function LL = get_pgabor_LL(K,Xmat,gain,offset,Robs,stim_ids,msize,hold_const,fixed_params)

params = zeros(10,1);
params(hold_const==1) = fixed_params;
params(hold_const==0) = K;

cur_mask1 = get_pgabor_mask(params(1:6),0,msize);
cur_mask2 = get_pgabor_mask(params(1:6),pi/2,msize);

mask1_out = Xmat*cur_mask1(:);
mask2_out = Xmat*cur_mask2(:);

energy_out = params(7)*sqrt(mask1_out.^2+mask2_out.^2);
lin_out = params(8)*mask1_out + params(9)*mask2_out;
tot_out = energy_out + lin_out + params(10);
g = tot_out(stim_ids).*gain + offset;

too_large = find(g > 100);
r = log(1+exp(g));
r(too_large) = g(too_large);

r(r < 1e-20) = 1e-20;

LL = sum(Robs.*log(r)-r);
Nspks = sum(Robs);
LL = -LL/Nspks;

end

