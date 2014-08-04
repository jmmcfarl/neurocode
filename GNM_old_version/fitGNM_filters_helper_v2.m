function [fit1 exitflag] = fitGNM_filters_helper_v2(fit0, X, linX, Robs, max_evals, method, silent, optTol, progTol,is_con)

%parse fitting parameters
if nargin < 7
    silent = 1;
end
options = [];
% if silent == 1
    options.Display = 'off';
% else
    options.Display = 'iter';
% end
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
if nmods == 0
    fprimes = [];
end
initial_params = [initial_params fit0.linK'];
initial_params(end+1) = fit0.spk_theta; %add constant offset term to params

lambda = zeros(size(initial_params));
for i = 1:nmods
    cur_set = (i-1)*klen + (1:klen);
    lambda(cur_set) = fit0.mods(i).lambda_L1x'.*ones(1,klen);
end
lambda(nmods*klen+1:end-1) = fit0.lambda_L1_lin;
lambda = lambda/sum(Robs);

if is_con == 0
    %     Aeq = zeros(1,nmods*klen+size(linX,2)+1);
    %     Aeq(end) = 1;
    %     Beq = 0;
    %     [params LL exitflag] = fmincon( @(K) fitGNM_filters_internal_v2(K, Robs, X, linX, fit0,fprimes), initial_params',[],[],Aeq,Beq,[],[],[],options);
    %     [params LL exitflag] = fminunc( @(K) fitGNM_filters_internal_v2(K, Robs, X, linX, fit0,fprimes), initial_params',options);
    if max(lambda) > 0
    [params LL] = L1General2_PSSas(@(K) fitGNM_filters_internal_v2(K, Robs, X, linX,fit0,fprimes),initial_params',lambda',options,0);
    exitflag = 1;
    else
        [params LL exitflag] = minFunc( @(K) fitGNM_filters_internal_v2(K, Robs, X, linX,fit0,fprimes), initial_params', options );
        
    end
else
    Aeq1 = zeros(kern_t,2*klen+size(linX,2)+1);
    for i = 1:kern_t
        Aeq1(i,i) = 1;
        Aeq1(i,i+3*kern_t) = -1;
    end
    Aeq2 = zeros(kern_t,2*klen+size(linX,2)+1);
    for i = 1:kern_t
        Aeq2(i,i+kern_t) = 1;
        Aeq2(i,i+2*kern_t) = 1;
    end
    Aeq = [Aeq1; Aeq2];
    Aeq(end+1,end) = 1;
    Beq = zeros(2*kern_t+1,1);
    [params LL exitflag] = fmincon( @(K) fitGNM_filters_internal_v2(K, Robs, X, linX, fit0,fprimes), initial_params',[],[],Aeq,Beq,[],[],[],options);
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
fit1.linK = params((nmods*klen+1):end-1);