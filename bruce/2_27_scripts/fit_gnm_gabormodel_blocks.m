function gabor_model = fit_gnm_gabormodel_blocks(X,Robs,blockids)

n_blocks = length(unique(blockids));
n_kerns = size(X,2);
w_init = zeros(n_kerns,n_blocks);
c_init = zeros(1,n_blocks);

initial_params = [w_init(:); c_init(:)];
options.Display = 'iter';
[params LL exitflag] = minFunc( @(K) gabormodel_blocks_LL(K, Robs, X,blockids), initial_params,options);

gabor_model.ws = reshape(params(1:n_kerns*n_blocks),n_kerns,n_blocks);
gabor_model.c = params((n_kerns*n_blocks+1):end);
gabor_model.LL = LL;

end

function [LL,LLgrad] = gabormodel_blocks_LL(K,Robs,X,blockids)

n_blocks = length(unique(blockids));
n_kerns = size(X,2);
w = K(1:n_kerns*n_blocks);
c = K((n_kerns*n_blocks+1):end);
w = reshape(w,n_kerns,n_blocks);

%initialize LL gradient
LLgradw = zeros(n_kerns,n_blocks);
LLgradc = zeros(n_blocks,1);
LL = 0;

for blockid = 1:n_blocks
    cur_set = find(blockids == blockid);
    g = X(cur_set,:)*w(:,blockid) + c(blockid);
    too_large = find(g > 100);
    expg = exp(g);
    r = log(1+expg);
    r(too_large) = g(too_large);
    
    r(r < 1e-20) = 1e-20; %minimum predicted rate
    
    cur_LL = sum(Robs(cur_set).*log(r)-r);
    Nspks = sum(Robs(cur_set));
    cur_LL = -cur_LL/Nspks;
    LL = LL + cur_LL;

    residual = (Robs(cur_set)./r - 1) .* expg ./ (1+expg);
    residual(too_large) = (Robs(cur_set(too_large))./r(too_large) - 1);
    
    % Calculate derivatives with respect to constant
    LLgradc(blockid) = sum(residual);
    LLgradw(:,blockid) = sum(bsxfun(@times,X(cur_set,:),residual));
    
    LLgradc(blockid) = -LLgradc(blockid)/Nspks;
    LLgradw(:,blockid) = -LLgradw(:,blockid)/Nspks;
end
LLgrad = zeros(size(K));
LLgrad(1:n_kerns*n_blocks) = LLgradw(:);
LLgrad((n_kerns*n_blocks+1):end) = LLgradc;

end