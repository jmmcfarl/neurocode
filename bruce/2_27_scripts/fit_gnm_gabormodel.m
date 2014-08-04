function gabor_model = fit_gnm_gabormodel(X,Robs)

n_kerns = size(X,2);
w_init = zeros(n_kerns,1);
c_init = 0;

initial_params = [w_init; c_init];
options.Display = 'off';
[params LL exitflag] = minFunc( @(K) gabormodel_LL(K, Robs, X), initial_params,options);

gabor_model.ws = params(1:end-1);
gabor_model.c = params(end);
gabor_model.LL = LL;

end

function [LL,LLgrad] = gabormodel_LL(K,Robs,X)

    n_kerns = size(X,2);
    c = K(end);
    w = K(1:end-1);
    
    g = X*w+c;
    too_large = find(g > 100);
    expg = exp(g);
    r = log(1+expg);
    r(too_large) = g(too_large);
    
    r(r < 1e-20) = 1e-20; %minimum predicted rate
    
    LL = sum(Robs.*log(r)-r);    

    residual = (Robs./r - 1) .* expg ./ (1+expg);
    residual(too_large) = (Robs(too_large)./r(too_large) - 1);
    
    %initialize LL gradient
    LLgrad = zeros(length(K),1);

    % Calculate derivatives with respect to constant
    LLgrad(end) = sum(residual);
  
    LLgrad(1:end-1) = sum(bsxfun(@times,X,residual));
    
    Nspks = sum(Robs);
    LL = -LL/Nspks;
    LLgrad = -LLgrad/Nspks;

end