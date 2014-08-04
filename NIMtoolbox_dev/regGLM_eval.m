function [LL, penLL, pred_rate, G] = regGLM_eval(glm,Robs,Xmat)
%
%[LL, penLL, pred_rate, G] = regGLM_eval(glm,Robs,Xmat)
% INPUTS: 
% glm: model struct
% Robs: binned spike counts
% Xmat: predictor matrix
% OUTPUTS:
% LL: log likelihood
% penLL: penalized log-likelihood
% pred_rate: predicted spike rate (unitless)
% G: generating function
%
%% Key parameters
[NT,kLen] = size(Xmat);
reg_params = glm.reg_params;
%make sure Robs is a column vector
if size(Robs,2) > size(Robs,1)
    Robs = Robs';
end

%% CREATE L2 REGULARIZATION MATRICES
n_L2mats = length(reg_params);
L2_mats = [];
for ii = 1:n_L2mats
    L2_mats{ii} = generate_L2_mat(reg_params(ii),kLen);
end

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
G = Xmat*glm.K + glm.theta;

%% Compute predicted firing rate and LL
if strcmp(glm.NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*glm.NL_params(1); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    pred_rate = glm.NL_params(2)*log(1+expg); %alpha*log(1+exp(gbeta))
    pred_rate(too_large) = glm.NL_params(2)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
%     pred_rate = pred_rate + glm.NL_params(3);
elseif strcmp(glm.NL_type,'exp')
    expg = exp(G);
    pred_rate = expg;
else
    error('invalid spk nl');
end
%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(pred_rate) < min_pred_rate
    pred_rate(pred_rate < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

LL = sum(Robs.* log(pred_rate) - pred_rate); %up to an overall constant

%% COMPUTE L2 PENALTIES
if ~isempty(L2_mats)
    L2_penalty = 0;
    for ii = 1:length(L2_mats)
        L2_penalty = L2_penalty + glm.lambdas(ii)*sum((L2_mats{ii} * glm.K).^2);
    end
else
    L2_penalty = 0;
end

%penalized LL
penLL = LL - L2_penalty;

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
Nspks = sum(Robs);
LL = LL/Nspks;
penLL = penLL/Nspks;
