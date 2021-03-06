function gabor_model = fit_gnm_gabormodel_blocks_v2(X,Robs,blockids,w_init,c_init)

n_blocks = length(unique(blockids));
n_kerns = size(X,2);
w_init = repmat(w_init,1,n_blocks);
c_init = repmat(c_init,1,n_blocks);

initial_params = [w_init(:); c_init(:)];

options.Display = 'off';
[params LL] = minFunc( @(K) gabormodel_blocks_LL_internal(K, Robs, X,blockids), initial_params,options);

gabor_model.ws = reshape(params(1:n_kerns*n_blocks),n_kerns,n_blocks);
gabor_model.c = params((n_kerns*n_blocks+1):end);
gabor_model.LL = LL;

end

