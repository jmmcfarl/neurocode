function [fit1 exitflag] = fitbetas_lexp_rot(fit0, g_mat, Robs, max_evals, method, silent,optTol,progTol)

options = [];
if silent == 1
    options.display = 'off';
else
    options.display = 'iter';
end
options.maxFunEvals = max_evals;
options.Method = method;
if nargin < 7
    options.optTol = 1e-4;
else
    options.optTol = optTol;
end
if nargin < 8
    options.progTol = 1e-6;
else
   options.progTol = progTol; 
end

if nargin < 6
    silent = 1;
end

% fsdim = fit0.mods(1).fsdim;
% if strcmp(fit0.basis,'white')
%     pkern_len = size(fit0.pix_conv_mat,2);
% elseif strcmp(fit0.basis,'pix')
%     pkern_len = size(X,2);
% end
% kern_t = pkern_len/fsdim;

stimlen = size(g_mat,1);

nmods = length(fit0.mods);
initial_params = [];
for m = 1:nmods
    cur_beta = fit0.mods(m).beta';
    initial_params = [initial_params cur_beta]; %add STCBcoefs to initial param vector
end
initial_params(end+1) = fit0.const; %add constant offset term to params

% [params LL exitflag] = minFunc( @(K) beta_lexp_LLinternal(K, Robs, g_mat, fit0), initial_params', options );

opts = optimset('GradObj','off','Display','iter','LargeScale','off','MaxIter',100,'MaxFunEvals',200,'TolFun',1e-8);
% [params LL exitflag output grad hessian] = fminunc( @(K) beta_lexp_LLinternal(K, Robs, g_mat, fit0), initial_params', opts );
[params LL exitflag output grad] = fmincon(@(K) beta_lexp_LLinternal(K, Robs, g_mat, fit0), initial_params',[],[],[],[],[0 0],[],[],opts);
% [f,df] = rot_lexp_LLinternal(initial_params', Robs, kern_output, fit0)

fit0.LP_seq = [fit0.LP_seq LL];
fit0.LP_ax = [fit0.LP_ax fit0.LP_ax(end)+1];

% Reassign variables
fit1 = fit0;
fit1.const = params(end);
for n = 1:nmods
    fit1.mods(n).beta = params(n);
end
