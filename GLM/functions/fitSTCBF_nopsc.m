function fit1 = fitSTCBF_nopsc(fit0, kern_output, Robs)
%

tolF = 1e-5;
tolX = 1e-8;
cs = [];

Nmods = length(fit0.mods);
initial_params = [];

for m = 1:Nmods
    NLx = fit0.mods(m).nlx;
    NL = fit0.mods(m).nly;
    
    %compute derivative of non-linearity
    fpr = zeros(1,length(NLx)-1);
    for n = 1:length(fpr)
        fpr(n) = (NL(n+1)-NL(n))/(NLx(n+1)-NLx(n));
    end    
    fprimes{m} = fpr;
    
    initial_params = [initial_params; fit0.mods(m).STCcf]; %add STCBcoefs to initial param vector
end
initial_params(end+1) = fit0.const; %add constant offset term to params

silent = 0;
if silent == 0
    opts = optimset('Algorithm','active-set','GradObj','on','LargeScale','off','Display','iter','MaxIter',30,'TolFun',tolF,'TolX',tolX);
else
    opts = optimset('Algorithm','active-set','GradObj','on','LargeScale','off','Display','off','MaxIter',30,'TolFun',tolF,'TolX',tolX);    
end

lamrange = [];

%%%%%%%% THE FIT %%%%%%%%
%if you don't want exc and sup filters to mix
[params LL eflag] = fminunc( @(K) STCBF_LLinternal_nopsc(K,Robs,kern_output,fit0, lamrange,fprimes), initial_params, opts );

% options.Method = 'cg';
% [params LL eflag] = minFunc( @(K) STCBF_LLinternal_nopsc(K,Robs,kern_output,fit0, lamrange,fprimes), initial_params,options);


% Reassign variables
fit1 = fit0;
% fit1.LP = LL;
% fit1.LL = adjustLL_STCBF(fit1,LL,sum(Robs));
% fit1.LL = LL - LLadjust( params(1:NK+NRP), lamrange, length(spks) );

NSTCdims = size(kern_output,2);
for n = 1:Nmods
    fit1.mods(n).STCcf = params((n-1)*NSTCdims + (1:NSTCdims));
    fit1.mods(n).k = fit1.STCbasis*fit1.mods(n).STCcf;
end

fit1 = normalizeRFs_STCB(fit1,kern_output);


% fprintf('Ent before: %.5f\n',total_kernel_entropy(fit0));
% fprintf('Ent after: %.5f\n',total_kernel_entropy(fit1));
% fprintf('Ent before: %.5f\n',average_kernel_std(fit0));
% fprintf('Ent after: %.5f\n',average_kernel_std(fit1));

% sdim = fit0.mods(1).SDIM;
% for i = 1:Nmods
%     old_std = kernel_std(fit0.STCbasis,fit0.mods(i).STCcf,sdim);
%     new_std = kernel_std(fit1.STCbasis,fit1.mods(i).STCcf,sdim);
%     fprintf('Std Before: %.5f   After: %.5f\n',old_std,new_std);
% end
