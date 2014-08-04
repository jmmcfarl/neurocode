function [gabor_params,LL] = fit_gabor_energy_mod_dualeye(XX,YY,init_params,Xmat_l,Xmat_r,Robs,hold_const,LB,UB)

options.Display = 'off';
options.Algorithm = 'sqp';
% options.tolfun = 1e-6;
% options.tolx = 1e-9;
% options.Algorithm = 'interior-point';
% options.MaxFunEvals = 1000;

K0 = init_params(hold_const == 0);
fixed_params = init_params(hold_const == 1);

[params LL] = fmincon( @(K) get_pgabor_LL(K, Xmat_l,Xmat_r, Robs,XX,YY,hold_const,fixed_params), ...
    K0,[],[],[],[],LB(hold_const==0),UB(hold_const==0),[],options);
gabor_params = init_params;
gabor_params(hold_const == 0) = params;

end

function LL = get_pgabor_LL(K,Xmat_l,Xmat_r,Robs,XX,YY,hold_const,fixed_params)
    
    params = zeros(9,1);
    params(hold_const==1) = fixed_params;
    params(hold_const==0) = K;

    cur_mask1 = get_pgabor_mask_v2(XX,YY,params(1:6),0);
    cur_mask2 = get_pgabor_mask_v2(XX,YY,params(1:6),pi/2);
    
    mask1_out = Xmat_l*cur_mask1(:);
    mask2_out = Xmat_l*cur_mask2(:);
    energy_out_l = params(7)*sqrt(mask1_out.^2+mask2_out.^2);

    mask1_out = Xmat_r*cur_mask1(:);
    mask2_out = Xmat_r*cur_mask2(:);
    energy_out_r = params(8)*sqrt(mask1_out.^2+mask2_out.^2);
    
    g = energy_out_l + energy_out_r + params(9);
    too_large = find(g > 100);
    r = log(1+exp(g));
    r(too_large) = g(too_large);
    
    r(r < 1e-20) = 1e-20;

    LL = sum(Robs.*log(r)-r);
    Nspks = sum(Robs);
    LL = -LL/Nspks;

end

