function [LL, LLgrad] = gabormodel_blocks_LL_internal(K,Robs,X,blockids)

n_blocks = length(unique(blockids));
n_kerns = size(X,2);
w = K(1:n_kerns*n_blocks);
c = K((n_kerns*n_blocks+1):end);
w = reshape(w,n_kerns,n_blocks);

% %initialize LL gradient
LLgradw = zeros(n_kerns,n_blocks);
LLgradc = zeros(n_blocks,1);
% LL = 0;

% for blockid = 1:n_blocks
%     cur_set = find(blockids == blockid);
%     g = X(cur_set,:)*w(:,blockid) + c(blockid);
%     too_large = find(g > 100);
%     expg = exp(g);
%     r = log(1+expg);
%     r(too_large) = g(too_large);
%
%     r(r < 1e-20) = 1e-20; %minimum predicted rate
%
%     cur_LL = sum(Robs(cur_set).*log(r)-r);
%     Nspks = sum(Robs(cur_set));
%     cur_LL = -cur_LL/Nspks;
%     LL = LL + cur_LL;
%
%     residual = (Robs(cur_set)./r - 1) .* expg ./ (1+expg);
%     residual(too_large) = (Robs(cur_set(too_large))./r(too_large) - 1);
%
%     % Calculate derivatives with respect to constant
%     LLgradc(blockid) = sum(residual);
%     LLgradw(:,blockid) = sum(bsxfun(@times,X(cur_set,:),residual));
%
%     LLgradc(blockid) = -LLgradc(blockid)/Nspks;
%     LLgradw(:,blockid) = -LLgradw(:,blockid)/Nspks;
% end
% LLgrad = zeros(size(K));
% LLgrad(1:n_kerns*n_blocks) = LLgradw(:);
% LLgrad((n_kerns*n_blocks+1):end) = LLgradc;


g = zeros(size(X,1),1);
for blockid = 1:n_blocks
    cur_set = find(blockids == blockid);
    g(cur_set) = X(cur_set,:)*w(:,blockid) + c(blockid);
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
    LLgradw(:,blockid) = sum(bsxfun(@times,X(cur_set,:),residual(cur_set)));
end
LLgradc = -LLgradc/Nspks;
LLgradw = -LLgradw/Nspks;

LLgrad = [LLgradw(:); LLgradc];
end