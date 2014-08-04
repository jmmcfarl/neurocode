function [fit1 exitflag] = fitGNM_filters_helper(fit0, X, Robs, max_evals, method, silent, optTol, progTol)

%parse fitting parameters
if nargin < 6
    silent = 1;
end
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

%stimulus parameters
[stimlen klen] = size(X);
kern_t = klen/fit0.stim_params.fsdim;

%compute initial fit parameters
nmods = length(fit0.mods);
initial_params = [];
for m = 1:nmods
    cur_kern = fit0.mods(m).k';
    initial_params = [initial_params cur_kern]; %add coefs to initial param vector
    
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
end
% initial_params(end+1) = fit0.const; %add constant offset term to params
initial_params(end+1) = fit0.spk_theta; %add constant offset term to params

lambda = [];
for i = 1:nmods
    lambda = [lambda fit0.mods(i).lambda_L1x'.*ones(1,klen)];
end
lambda = [lambda 0]/sum(Robs);

if max(lambda) > 0
    [params LL] = L1General2_PSSas(@(K) fitGNM_filters_internal(K, Robs, X, fit0,fprimes),initial_params',lambda',options,silent);
    exitflag = 1;
else
    [params LL exitflag] = minFunc( @(K) fitGNM_filters_internal(K, Robs, X, fit0,fprimes), initial_params', options );
end

fit0.LP_seq = [fit0.LP_seq LL];
% Reassign variables
fit1 = fit0;
% fit1.const = params(end);
fit1.spk_theta = params(end);
for n = 1:nmods
    cur_kern = params((n-1)*klen+(1:klen));
    kern_out = X * cur_kern;
    fit1.mods(n).k = cur_kern(:);
end
