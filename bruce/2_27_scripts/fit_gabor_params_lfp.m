function [gabor_params,LL,err_sigma] = fit_gabor_params_lfp(init_params,Xmat,obs_amps,msize,hold_const,LB,UB)

options.Display = 'off';
options.Algorithm = 'sqp';
% options.MaxFunEvals = 1000;

K0 = init_params(hold_const == 0);
fixed_params = init_params(hold_const == 1);

% [params LL exitflag] = fminunc( @(K) get_pgabor_LL(K, Xmat, Robs, msize,hold_const,fixed_params), K0,options);
[params LL] = fmincon( @(K) get_pgabor_LL(K, Xmat, obs_amps, msize,hold_const,fixed_params), ...
    K0,[],[],[],[],LB(hold_const==0),UB(hold_const==0),[],options);
gabor_params = init_params;
gabor_params(hold_const == 0) = params;

cur_mask1 = get_pgabor_mask(gabor_params(1:6),0,msize);
cur_mask2 = get_pgabor_mask(gabor_params(1:6),pi/2,msize);

mask1_out = Xmat*cur_mask1(:);
mask2_out = Xmat*cur_mask2(:);

energy_out = gabor_params(7)*sqrt(mask1_out.^2+mask2_out.^2);
g = energy_out + gabor_params(8);
err_sigma = std(g-obs_amps');
end

function LL = get_pgabor_LL(K,Xmat,obs_amps,msize,hold_const,fixed_params)
    
    params = zeros(8,1);
    params(hold_const==1) = fixed_params;
    params(hold_const==0) = K;

    cur_mask1 = get_pgabor_mask(params(1:6),0,msize);
    cur_mask2 = get_pgabor_mask(params(1:6),pi/2,msize);
    
    mask1_out = Xmat*cur_mask1(:);
    mask2_out = Xmat*cur_mask2(:);

    energy_out = params(7)*sqrt(mask1_out.^2+mask2_out.^2);   
    g = energy_out + params(8);

    LL = sum((g-obs_amps').^2);
    
end

