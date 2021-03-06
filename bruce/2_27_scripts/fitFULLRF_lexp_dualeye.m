function [fit1 exitflag] = fitFULLRF_lexp_dualeye(fit0, X, Robs, max_evals, method, silent,optTol,progTol)

options = [];
if silent == 1
    options.Display = 'off';
else
    options.Display = 'iter';
end
options.maxFunEvals = max_evals;
options.maxIter = max_evals;
options.Method = method;
if nargin < 7
    options.optTol = 1e-5;
else
    options.optTol = optTol;
end
if nargin < 8
    options.progTol = 1e-7;
else
   options.progTol = progTol; 
end

if nargin < 6
    silent = 1;
end

fsdim = fit0.mods(1).fsdim;
if strcmp(fit0.basis,'white')
    pkern_len = size(fit0.pix_conv_mat,2);
elseif strcmp(fit0.basis,'pix')
    pkern_len = size(X,2);
end
kern_t = pkern_len/fsdim;

stimlen = size(X,1);

nmods = length(fit0.mods);
initial_params = [];
for m = 1:nmods
    cur_pkern = fit0.mods(m).pix';
    initial_params = [initial_params cur_pkern]; %add STCBcoefs to initial param vector
end
initial_params(end+1) = fit0.const; %add constant offset term to params

lambda = [];
for i = 1:nmods
    lambda = [lambda fit0.mods(i).lambda_L1x'.*ones(1,pkern_len)];
end
lambda = [lambda 0]/sum(Robs);

if max(lambda) > 0
    [params LL] = L1General2_PSSas(@(K) FULLBF_lexp_LLinternal_dualeye(K, Robs, X, fit0),initial_params',lambda',options,0);
    exitflag = 1;
else
    [params LL exitflag] = minFunc( @(K) FULLBF_lexp_LLinternal_dualeye(K, Robs, X, fit0), initial_params', options );
end

%%%%%%%%%%%%%%%%%%%%%%
% if silent == 0
% %     opts = optimset('GradObj','on','Algorithm','active-set','Display','iter','MaxIter',200,'MaxFunEvals',10000,'TolFun',1e-7);
%     opts = optimset('GradObj','on','LargeScale','on','Display','iter','MaxIter',200,'MaxFunEvals',10000,'TolFun',1e-7);
% else
%     opts = optimset('GradObj','on','Algorithm','active-set','MaxIter',200,'Display','off','MaxFunEvals',10000,'TolFun',1e-7);
% end
% [params LL exitflag] = fminunc( @(K) FULLBF2d_LLinternal(K, Robs, X, fit0), initial_params', opts );



fit0.LP_seq = [fit0.LP_seq LL];
fit0.LP_ax = [fit0.LP_ax fit0.LP_ax(end)+1];

% Reassign variables
fit1 = fit0;
fit1.const = params(end);
for n = 1:nmods
    cur_kern = params((n-1)*pkern_len+(1:pkern_len));
    if strcmp(fit0.basis,'white')
        cur_wkern = cur_kern'*fit0.kern_conv_mat;
        fit1.mods(n).k = cur_wkern(:);
        fit1.mods(n).pix = cur_kern(:);
    elseif strcmp(fit0.basis,'pix')
        kern_out = X * cur_kern;
        fit1.mods(n).k = cur_kern(:);
        fit1.mods(n).pix = cur_kern(:);
    end
end
