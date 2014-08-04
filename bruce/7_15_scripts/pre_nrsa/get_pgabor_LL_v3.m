function LL = get_pgabor_LL_v3(K,Xmat,Robs,XX,YY,hold_const,fixed_params)
    
    params = zeros(10,1);
    params(hold_const==1) = fixed_params;
    params(hold_const==0) = K;

    cur_mask1 = get_pgabor_mask_v2(XX,YY,params(1:6),0);
    cur_mask2 = get_pgabor_mask_v2(XX,YY,params(1:6),pi/2);
    
    mask1_out = Xmat*cur_mask1(:);
    mask2_out = Xmat*cur_mask2(:);

    energy_out = params(7)*sqrt(mask1_out.^2+mask2_out.^2);
    lin_out = params(8)*mask1_out + params(9)*mask2_out;
    
    g = energy_out + lin_out + params(10);
    too_large = find(g > 100);
    r = log(1+exp(g));
    r(too_large) = g(too_large);
       
%     r = exp(g);
    
    r(r < 1e-20) = 1e-20;

    LL = sum(Robs.*log(r)-r);
    Nspks = sum(Robs);
    LL = -LL/Nspks;

end
