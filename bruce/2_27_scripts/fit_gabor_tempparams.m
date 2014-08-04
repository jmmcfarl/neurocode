function [temp_params,LL] = fit_gabor_tempparams(init_params,gabor_params,Xmat,Robs,msize,blockids)

options.Display = 'iter';
options.Algorithm = 'sqp';

K0 = init_params(:);

[params LL] = fminunc( @(K) get_tempgabor_LL(K, gabor_params, Xmat, Robs, msize,blockids),K0,options);
temp_params = params;

end

function [LL,LLgrad] = get_tempgabor_LL(K, gabor_params, Xmat,Robs,msize,blockids)

n_blocks = length(unique(blockids));
kmat = reshape(K,n_blocks,4);

w_e = kmat(:,1);
w_s0 = kmat(:,2);
w_s90 = kmat(:,3);
c = kmat(:,4);

cur_mask1 = get_pgabor_mask(gabor_params,0,msize);
cur_mask2 = get_pgabor_mask(gabor_params,pi/2,msize);

mask1_out = Xmat*cur_mask1(:);
mask2_out = Xmat*cur_mask2(:);

%initialize LL gradient
LL = 0;
for i = 1:n_blocks
    cur_set = find(blockids == i);
    energy_out = sqrt(mask1_out(cur_set).^2+mask2_out(cur_set).^2);
    lin_out = w_s0(i)*mask1_out(cur_set) + w_s90(i)*mask2_out(cur_set);    
    g = lin_out + w_e(i)*energy_out + c(i);
    
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
    LLgrad_c(i) = -sum(residual)/Nspks;
    LLgrad_we(i) = -sum(bsxfun(@times,energy_out,residual))/Nspks;
    LLgrad_ws0(i) = -sum(bsxfun(@times,mask1_out(cur_set),residual))/Nspks;
    LLgrad_ws90(i) = -sum(bsxfun(@times,mask2_out(cur_set),residual))/Nspks;  
        
end

LLgrad = [LLgrad_we(:) LLgrad_ws0(:) LLgrad_ws90(:) LLgrad_c(:)];
LLgrad = LLgrad(:);
end

