function [beta,fin_LL] = smooth_regress_2d_doubkern(response,pred,init_beta,smooth_lambda,sparse_lambda,kern_inds)

options.Display = 'iter';
options.optTol = 1;
options.maxIter = 500;
if sparse_lambda == 0
    [bestk,LP] = minFunc(@(K) smooth_LL_internal(response,pred,K,smooth_lambda,kern_inds),init_beta,options);
else
    lambda = sparse_lambda*ones(size(init_beta));
    [bestk,LP] = L1General2_PSSas(@(K) smooth_LL_internal(response,pred,K,smooth_lambda,kern_inds),init_beta,lambda,options,0,100);
end
beta = bestk;

pred_Y = pred*beta;
resid = response-pred_Y;
fin_LL = sum(resid.^2);
end

function [LL LLgrad] = smooth_LL_internal(Y,X,beta,lambda,kern_inds)

[n_pts,tot_n_dims] = size(X);
kern1_inds = find(kern_inds==1);
kern2_inds = find(kern_inds==2);
n_dims1 = length(kern1_inds);
n_dims2 = length(kern2_inds);
side_len1 = sqrt(n_dims1);
side_len2 = sqrt(n_dims2);

pred_Y = X*beta;
resid = Y-pred_Y;
LL = sum(resid.^2);

lapl_oper = [0 1 0;1 -4 1;0 1 0];

kern_mat1 = reshape(beta(kern1_inds),side_len1,side_len1);
kern_mat1 = [zeros(side_len1,1) kern_mat1 zeros(side_len1,1)];
kern_mat1 = [zeros(1,side_len1+2); kern_mat1; zeros(1,side_len1+2)];
lapl1 = conv2(kern_mat1,lapl_oper,'same');
kern_mat2 = reshape(beta(kern2_inds),side_len2,side_len2);
kern_mat2 = [zeros(side_len2,1) kern_mat2 zeros(side_len2,1)];
kern_mat2 = [zeros(1,side_len2+2); kern_mat2; zeros(1,side_len2+2)];
lapl2 = conv2(kern_mat2,lapl_oper,'same');

temp = lapl1(2:end-1,2:end-1);
lapl_penalty1 = lambda*sum(temp(:).^2);
temp = lapl2(2:end-1,2:end-1);
lapl_penalty2 = lambda*sum(temp(:).^2);
LL = LL + lapl_penalty1 + lapl_penalty2;

LLgrad = -2*sum(bsxfun(@times,X,resid));

LLgrad_pen = zeros(size(LLgrad));

temp = zeros(size(lapl1));
temp = conv2(lapl1,lapl_oper,'same');
temp = temp(2:end-1,2:end-1);
LLgrad_pen(kern1_inds) = 2*lambda*temp(:);
temp = zeros(size(lapl2));
temp = conv2(lapl2,lapl_oper,'same');
temp = temp(2:end-1,2:end-1);
LLgrad_pen(kern2_inds) = 2*lambda*temp(:);

LLgrad = LLgrad' + LLgrad_pen';

end

