function [gabor_params,LL] = fit_gabor_params_lfp_withcon(init_params,Xmat,obs_amps,msize,hold_const,LB,UB)

options.Display = 'iter';
options.Algorithm = 'sqp';
% options.MaxFunEvals = 1000;

K0 = init_params(hold_const == 0);
fixed_params = init_params(hold_const == 1);

% [params LL exitflag] = fminunc( @(K) get_pgabor_LL(K, Xmat, Robs, msize,hold_const,fixed_params), K0,options);
[params LL] = fmincon( @(K) get_pgabor_LL(K, Xmat, obs_amps, msize,hold_const,fixed_params), ...
    K0,[],[],[],[],LB(hold_const==0),UB(hold_const==0),[],options);
gabor_params = init_params;
gabor_params(hold_const == 0) = params;
end

function LL = get_pgabor_LL(K,Xmat,obs_amps,msize,hold_const,fixed_params)
    
    params = zeros(9,1);
    params(hold_const==1) = fixed_params;
    params(hold_const==0) = K;

    cur_mask1 = get_pgabor_mask(params(1:6),0,msize);
    cur_mask2 = get_pgabor_mask(params(1:6),pi/2,msize);
    cur_gmask = get_pgauss_mask(params(1:6),msize);
    
    mask1_out = Xmat*cur_mask1(:);
    mask2_out = Xmat*cur_mask2(:);
    cont_out = sqrt(Xmat.^2*cur_gmask(:));
    
    energy_out = sqrt(mask1_out.^2+mask2_out.^2);   
    
    cont_out = zscore(cont_out);
    energy_out = zscore(energy_out);
    
    g = params(7)*energy_out + params(8)*cont_out + params(9);

    LL = sum((g-obs_amps').^2);
    
end

