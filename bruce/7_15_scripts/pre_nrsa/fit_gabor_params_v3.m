function [gabor_params,LL] = fit_gabor_params_v3(XX,YY,init_params,Xmat,Robs,hold_const,LB,UB)

options.Display = 'off';
options.Algorithm = 'sqp';
% options.tolfun = 1e-6;
% options.tolx = 1e-9;
% options.Algorithm = 'interior-point';
% options.MaxFunEvals = 1000;

K0 = init_params(hold_const == 0);
fixed_params = init_params(hold_const == 1);

[params LL] = fmincon( @(K) get_pgabor_LL_v3(K, Xmat, Robs,XX,YY,hold_const,fixed_params), ...
    K0,[],[],[],[],LB(hold_const==0),UB(hold_const==0),[],options);
gabor_params = init_params;
gabor_params(hold_const == 0) = params;

end


