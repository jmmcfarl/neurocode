function [betas,LP] = smoothed_regression(Y,X,init_beta,pens,ranges,lambda)

options.Display = 'iter';
options.optTol = 1;
options.maxIter = 100;
lambda = ones(size(init_beta))*lambda;
% lambda(end-1:end) = 0;
% [bestk,LP] = minFunc(@(K) smooth_LL_internal(Y,X,K,pens,ranges),init_beta,options);
[bestk,LP] = L1General2_PSSas(@(K) smooth_LL_internal(Y,X,K,pens,ranges),init_beta,lambda,options,0,100);
betas = bestk;

end

function [LL LLgrad] = smooth_LL_internal(Y,X,beta,pens,ranges)

n_pts = length(Y);
pred_Y = X*beta;
resid = Y-pred_Y;
LL = sum(resid.^2);
LLgrad = -2*sum(bsxfun(@times,X,resid));
grad_pen = zeros(size(LLgrad));
for i = 1:size(ranges,1)
    range = ranges(i,1):ranges(i,2);
    chunk = beta(range)'; %'
    
%     chunk2 = [2*chunk(1)-chunk(2) chunk 2*chunk(end)-chunk(end-1)];
    chunk2 = [0 0 chunk 0 0];
    
    smooth_pen = pens(i) * sum((2*chunk2(2:end-1) - chunk2(1:end-2) - chunk2(3:end)).^2);
    LL = LL + smooth_pen;
    
    grad_pen(range) = grad_pen(range) - 2*pens(i) * (6*chunk2(3:end-2) - 4*chunk2(2:end-3) - 4*chunk2(4:end-1) + chunk2(1:end-4)+chunk2(5:end));
end
LLgrad = LLgrad - grad_pen;
LLgrad = LLgrad';
end