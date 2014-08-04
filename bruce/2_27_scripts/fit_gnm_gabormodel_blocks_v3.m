function [scale,offset,LL] = fit_gnm_gabormodel_blocks_v3(X,Robs,blockids)

n_blocks = length(unique(blockids));
w_init = zeros(n_blocks,1);
c_init = zeros(n_blocks,1);

initial_params = [w_init(:); c_init(:)];

options.Display = 'off';
[params LL] = minFunc( @(K) gabormodel_blocks_LL_internal(K, Robs, X, blockids), initial_params,options);

scale = params(1:n_blocks);
offset = params((n_blocks+1):end);

end

function [LL, LLgrad] = gabormodel_blocks_LL_internal(K,Robs,X,blockids)

n_blocks = length(unique(blockids));
w = K(1:n_blocks);
c = K((n_blocks+1):end);

% %initialize LL gradient
LLgradw = zeros(n_blocks,1);
LLgradc = zeros(n_blocks,1);
g = zeros(size(X,1),1);
for blockid = 1:n_blocks
    cur_set = find(blockids == blockid);
    g(cur_set) = X(cur_set,:)*w(blockid) + c(blockid);
end
too_large = find(g > 100);
expg = exp(g);
r = log(1+expg);
r(too_large) = g(too_large);
r(r < 1e-20) = 1e-20; %minimum predicted rate

Nspks = sum(Robs);
LL = sum(Robs.*log(r)-r);
LL = -LL/Nspks;

residual = (Robs./r - 1) .* expg ./ (1+expg);
residual(too_large) = (Robs(too_large)./r(too_large) - 1);
for blockid = 1:n_blocks
    cur_set = find(blockids == blockid);
    % Calculate derivatives with respect to constant
    LLgradc(blockid) = sum(residual(cur_set));
    LLgradw(blockid) = sum(X(cur_set).*residual(cur_set));
end
LLgradc = -LLgradc/Nspks;
LLgradw = -LLgradw/Nspks;

LLgrad = [LLgradw; LLgradc];
end