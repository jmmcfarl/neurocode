function [beta,fin_LL] = smooth_regress_2d(response,pred,init_beta,smooth_lambda,sparse_lambda)

options.Display = 'iter';
options.optTol = 1;
options.maxIter = 500;
if sparse_lambda == 0
[bestk,LP] = minFunc(@(K) smooth_LL_internal(response,pred,K,smooth_lambda),init_beta,options);
else
    lambda = sparse_lambda*ones(size(init_beta));
 [bestk,LP] = L1General2_PSSas(@(K) smooth_LL_internal(response,pred,K,smooth_lambda),init_beta,lambda,options,0,100);   
end
beta = bestk;

pred_Y = pred*beta;
resid = response-pred_Y;
fin_LL = sum(resid.^2);
end

function [LL LLgrad] = smooth_LL_internal(Y,X,beta,lambda)

[n_pts,n_dims] = size(X);
side_len = sqrt(n_dims);
pred_Y = X*beta;
resid = Y-pred_Y;
LL = sum(resid.^2);

lapl_oper = [0 1 0;1 -4 1;0 1 0];
kern_mat = reshape(beta,side_len,side_len);
kern_mat = [zeros(side_len,1) kern_mat zeros(side_len,1)];
kern_mat = [zeros(1,side_len+2); kern_mat; zeros(1,side_len+2)];
lapl = conv2(kern_mat,lapl_oper,'same');
temp = lapl(2:end-1,2:end-1);
lapl_penalty = lambda*sum(temp(:).^2);
LL = LL + lapl_penalty;

LLgrad = -2*sum(bsxfun(@times,X,resid));

temp = zeros(size(lapl));
temp = conv2(lapl,lapl_oper,'same');
temp = temp(2:end-1,2:end-1);
LLgrad_pen = 2*lambda*temp(:);

LLgrad = LLgrad' + LLgrad_pen;

end

