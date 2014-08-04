function [fitmod] = regGLM_fit(Xmat,Robs,reg_params,lambdas,NL_type,NL_params,silent,init_params,const_params)
%
%[fitmod] = regGLM_fit(Xmat,Robs,<reg_params>,<lambdas>,<NL_type>,<NL_params>,<silent>,<init_params>,<const_params>)
%
%INPUTS:
% Xmat: [TxD] Xmatrix of predictors
% Robs: binned spike count (column) vector
% reg_params: structure of parameters defining the type of regularization used (see "create_L2_params.m")
% lambdas: vector of regularization strengths (hyperparameters). Must be same length as the structure array reg_params
% NL_type: type of spiking NL ('logexp' or 'exp')
% NL_params: parameeters of spiking NL function 
% silent: 0 or 1 for silent or verbose output
% init_params: Dx1 vector of initial parameters (shouldn't matter since this is a convex problem GLM, but can be good to avoid computational issues)
% const_params: vector of indices of parameters to hold constant
% OUTPUTS: 
% fitmod: structure of model parameters

%% PROCESS INPUTS
if nargin < 3
    reg_params = [];
end
if nargin < 4
    lambdas = [];
end
if nargin < 5 || isempty(NL_type)
    NL_type = 'logexp';
end
if nargin < 6 || isempty(NL_params)
    if strcmp(NL_type,'logexp')
        NL_params = [1 1]; %[beta alpha]
    elseif strcmp(NL_type,'exp')
        NL_params = [];
    else
        error('Unsupported NL type');
    end
end
if nargin < 7 || isempty(silent)
    silent = 1;
end
if nargin < 8
    init_params = [];
end
if nargin < 9 
    const_params = [];
end

if length(lambdas) ~= length(reg_params)
    error('Mismatch between specified lambdas and reg_params');
end

%make sure Robs is a column vector
if size(Robs,2) > size(Robs,1)
    Robs = Robs';
end

[NT,kLen] = size(Xmat); %stimulus dimensions

%% INITIAL PARAMETERS
if isempty(init_params)
    K0 = randn(kLen,1)/kLen*0.1;
else
    K0 = init_params(:);
end
if length(K0) == size(Xmat,2)
    K0(end+1) = 0; %initial constant term
end
K0 = K0(:);
%% GENERATE REGULARIZATION MATRICES
n_L2mats = length(reg_params);
L2_mats = [];
for ii = 1:n_L2mats
    L2_mats{ii} = generate_L2_mat(reg_params(ii),kLen);
end

%% identify any constraints
use_con = 0;
Aeq = []; %initialize constraint matrices
if ~isempty(const_params)
    Aeq = zeros(length(const_params),length(K0));
    for ii = 1:length(const_params)
        Aeq(ii,const_params(ii)) = 1;
    end
    use_con = 1;
end
beq = K0(const_params);  

%%
optim_params.MaxFunEvals = 100*length(K0);
optim_params.MaxIter = 1e3;
optim_params.Display = 'off';
if silent == 0
    optim_params.Display = 'iter';
end

if use_con == 1
    optim_params.GradObj = 'on';
    optim_params.LargeScale = 'off';
    optim_params.Algorithm = 'active-set';
    
   [Kopt] = fmincon( @(K) regGLM_fit_internal(K, Xmat, Robs, L2_mats, lambdas, NL_type, NL_params), K0, [],[],Aeq,beq,[],[],[],optim_params);
    
else
    if exist('minFunc','file') == 2 %try to use Mark Schmidt's optimizer
        
        %if using Mark Schmidt's optimization, some differences in option parameters
        optim_params.optTol = 1e-5;
        optim_params.progTol = 1e-7;
        optim_params.Method = 'lbfgs';
        [Kopt] = minFunc( @(K) regGLM_fit_internal(K, Xmat, Robs, L2_mats, lambdas, NL_type, NL_params), K0, optim_params);
        
    else %if using Matlab Optim toolbox:
        %default optimization parameters
        optim_params.LargeScale = 'off';
        optim_params.TolFun = 1e-6;
        optim_params.TolX = 1e-6;
        optim_params.HessUpdate = 'bfgs';
        optim_params.GradObj = 'on';
        
        [Kopt] = fminunc( @(K) regGLM_fit_internal(K, Xmat, Robs, L2_mats, lambdas, NL_type, NL_params), K0, optim_params);
    end
end
%% PARSE MODEL FIT
fitmod.K = Kopt(1:end-1);
fitmod.theta = Kopt(end);
fitmod.NL_type = NL_type;
fitmod.NL_params = NL_params;
fitmod.reg_params = reg_params;
fitmod.lambdas = lambdas;

[LL, penLL] = regGLM_eval(fitmod,Robs,Xmat);

fitmod.LL = LL;
fitmod.penLL = penLL;

end

function [LL, LLgrad] = regGLM_fit_internal(params, Xmat, Robs, L2_mats, lambdas, NL_type, NL_params)
%

%% USEFUL PARAMETERS
[NT,kLen] = size(Xmat);
K = params(1:end-1);

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
G = Xmat*K + params(end);

%% Compute predicted firing rate
if strcmp(NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*NL_params(1); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    r = NL_params(2)*log(1+expg); %alpha*log(1+exp(gbeta))
    r(too_large) = NL_params(2)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
%     r = r + NL_params(3);
elseif strcmp(NL_type,'exp')
    expg = exp(G);
    r = expg;
else
    error('invalid spk nl');
end
%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(r) < min_pred_rate
    r(r < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

%% COMPUTE LL and LL gradient
LL = sum(Robs.* log(r) - r); %up to an overall constant

%'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
if strcmp(NL_type,'logexp')
    residual = NL_params(2)*NL_params(1)*(Robs./r - 1) .* expg ./ (1+expg);
    residual(too_large) = NL_params(2)*NL_params(1)*(Robs(too_large)./r(too_large) - 1);
elseif strcmp(NL_type,'exp')
    residual = Robs - r;
else
    error('Unsupported spiking NL')
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant term (theta)
LLgrad(end) = sum(residual);

% Calculate derivative with respect to linear term
LLgrad(1:end-1) = residual'*Xmat;

%% COMPUTE L2 PENALTIES AND ASSOCIATED CONTRIBUTIONS TO THE LL GRADIENT
L2_penalty = 0;
L2_grad_pen = zeros(kLen,1);
if ~isempty(L2_mats)
    for ii = 1:length(L2_mats)
        L2_penalty = L2_penalty + lambdas(ii)*sum((L2_mats{ii} * K).^2);
        L2_grad_pen = L2_grad_pen + 2*lambdas(ii)*(L2_mats{ii}' * L2_mats{ii} * K);
    end
end
LL = LL - L2_penalty;
LLgrad(1:end-1) = LLgrad(1:end-1) - L2_grad_pen;

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;
end