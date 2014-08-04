function [fit1 exitflag] = fitFULLRF_lexp_rot(fit0, kern_output, Robs, max_evals, method, silent,optTol,progTol)

options = [];
if silent == 1
    options.display = 'off';
else
    options.display = 'iter';
end
options.maxFunEvals = max_evals;
options.Method = method;
if nargin < 7
    options.optTol = 1e-6;
else
    options.optTol = optTol;
end
if nargin < 8
    options.progTol = 1e-8;
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

stimlen = size(kern_output,1);

nmods = length(fit0.mods);
initial_params = [];
for m = 1:nmods
    NLx = fit0.mods(m).nlx;
    NL = fit0.mods(m).nly;
    
    if strcmp(fit0.mods(m).nltype,'uncon')
        %compute derivative of non-linearity
        fpr = zeros(1,length(NLx)-1);
        for n = 1:length(fpr)
            fpr(n) = (NL(n+1)-NL(n))/(NLx(n+1)-NLx(n));
        end
        fprimes{m} = fpr;
    else
        fprimes{m} = [];
    end
    cur_stc_kern = fit0.mods(m).STCcf';
    initial_params = [initial_params cur_stc_kern]; %add STCBcoefs to initial param vector
end
initial_params(end+1) = fit0.const; %add constant offset term to params

[params LL exitflag] = minFunc( @(K) rot_lexp_LLinternal(K, Robs, kern_output, fit0, fprimes), initial_params', options );

% opts = optimset('GradObj','on','Display','iter','LargeScale','off','MaxIter',500,'MaxFunEvals',1000,'TolFun',1e-7,'TolX',1e-9);
% [params LL exitflag output grad hessian] = fminunc( @(K) rot_lexp_LLinternal(K, Robs, kern_output, fit0), initial_params', opts );



% [f,df] = rot_lexp_LLinternal(initial_params', Robs, kern_output, fit0)

fit0.LP_seq = [fit0.LP_seq LL];
fit0.LP_ax = [fit0.LP_ax fit0.LP_ax(end)+1];

% Reassign variables
fit1 = fit0;
fit1.const = params(end);
NSTCdims = size(kern_output,2);
for n = 1:nmods
    fit1.mods(n).STCcf = params((n-1)*NSTCdims + (1:NSTCdims));
    fit1.mods(n).k = fit1.STCbasis*fit1.mods(n).STCcf;
end
