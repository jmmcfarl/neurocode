function [LL, LLgrad] = phasedep_glmfit(params, Robs, Xstim,XPBx,XLin,L2_mat,reg_lambda,stim_params,PBx)
%

%% USEFUL PARAMETERS
filtLen = prod(stim_params.stim_dims);
nBars = squeeze(stim_params.stim_dims(2));
nPBs = length(PBx);
lin_dims = stim_params.lin_dims;

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
% theta = params(end); %offset
kmat = reshape(params(1:filtLen*nPBs),filtLen,nPBs);
KPLin = params((filtLen*nPBs+1):(filtLen*nPBs + nPBs));
KLin = params((filtLen*nPBs + nPBs+1):(filtLen*nPBs + nPBs+lin_dims));
% KLin = params((filtLen*nPBs + 1):(filtLen*nPBs + lin_dims+1));
filt_outs = Xstim*kmat;

% G = XPBx*Klin + theta;%initialize overall generating function G
G = XPBx*KPLin + XLin*KLin;%initialize overall generating function G
% G = XLin*KLin;%initialize overall generating function G

%%
G = G + sum(filt_outs.*XPBx,2);
%% Compute predicted firing rate
max_gbeta = 50; %to prevent numerical overflow
expg = exp(G);
too_large = (G > max_gbeta);
r = log(1+expg); %alpha*log(1+exp(gbeta))
r(too_large) = G(too_large); %log(1+exp(x)) ~ x in limit of large x

%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(r) < min_pred_rate
    r(r < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

%% COMPUTE LL and LL gradient
LL = sum(Robs.* log(r) - r); %up to an overall constant

%'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
residual = (Robs./r - 1) .* expg ./ (1+expg);
residual(too_large) = (Robs(too_large)./r(too_large) - 1);

%initialize LL gradient
LLgrad = zeros(length(params),1)';

% % Calculate derivatives with respect to constant term (theta)
% LLgrad(end) = sum(residual);

LLgrad((filtLen*nPBs+1):(filtLen*nPBs + nPBs)) = residual'*XPBx;
LLgrad((filtLen*nPBs + nPBs + 1):(filtLen*nPBs + nPBs+lin_dims)) = residual'*XLin;
% LLgrad((filtLen*nPBs + 1):(filtLen*nPBs + lin_dims+1)) = residual'*XLin;

%%
full_resid = bsxfun(@times,XPBx,residual);
cur_LLgradmat = Xstim'*full_resid;
LLgrad(1:filtLen*nPBs) = cur_LLgradmat(:);
LLgrad = LLgrad';
%% COMPUTE L2 PENALTIES AND ASSOCIATED CONTRIBUTIONS TO THE LL GRADIENT
L2_penalty = reg_lambda*sum((L2_mat * kmat(:)).^2);
L2_grad_pen = 2*reg_lambda*(L2_mat' * L2_mat * kmat(:));

LL = LL - L2_penalty;
LLgrad(1:filtLen*nPBs) = LLgrad(1:filtLen*nPBs) - L2_grad_pen;

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

