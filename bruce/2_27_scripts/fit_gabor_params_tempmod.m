function [gabor_params,LL] = fit_gabor_params_tempmod(init_params,temp_mod,Xmat,Tmat,stim_ids,Robs,msize,hold_const,LB,UB)

options.Display = 'iter';
options.Algorithm = 'sqp';
% options.MaxFunEvals = 1000;

K0 = init_params(hold_const == 0);
fixed_params = init_params(hold_const == 1);

% [params LL exitflag] = fminunc( @(K) get_pgabor_LL(K, Xmat, Robs, msize,hold_const,fixed_params), K0,options);
[params LL] = fmincon( @(K) get_pgabor_LL(K, temp_mod, Xmat, Tmat, stim_ids, Robs, msize,hold_const,fixed_params), ...
    K0,[],[],[],[],LB(hold_const==0),UB(hold_const==0),[],options);
gabor_params = init_params;
gabor_params(hold_const == 0) = params;

end

function LL = get_pgabor_LL(K,temp_mod,Xmat,Tmat, stim_ids, Robs,msize,hold_const,fixed_params)
    
    params = zeros(7,1);
    params(hold_const==1) = fixed_params;
    params(hold_const==0) = K;

    cur_mask1 = get_pgabor_mask(params(1:6),0,msize);
    cur_mask2 = get_pgabor_mask(params(1:6),pi/2,msize);
%     cur_mask3 = get_pgabor_mask(params(1:5),params(8),msize);
    
    mask1_out = Xmat*cur_mask1(:);
    mask2_out = Xmat*cur_mask2(:);
%     mask3_out = Xmat*cur_mask3(:);

    energy_out = params(7)*sqrt(mask1_out.^2+mask2_out.^2);
    lin_out = params(8)*mask1_out + params(9)*mask2_out;
    
    g_stim = energy_out + lin_out;
    g_stim_temp = g_stim(stim_ids);
    temp_stim = (Tmat*temp_mod.stim_kern').*g_stim_temp;
    
    temp_sac = Tmat*temp_mod.sac_kern';
    
    g = temp_sac + temp_stim + params(end);
    
    too_large = find(g > 100);
    r = log(1+exp(g));
    r(too_large) = g(too_large);
    
    r(r < 1e-20) = 1e-20;
    
    LL = sum(Robs.*log(r)-r);
    Nspks = sum(Robs);
    LL = -LL/Nspks;

end

