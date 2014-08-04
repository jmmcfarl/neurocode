function [fitp] = fit_GLM_sinphase_model(X, Robs, K0, silent, stim_params, reg_params,use_time,NL_type)
%

if nargin < 7
    use_time = 0;
end
if nargin < 8
    NL_type = 0;
end
[NT NPAR] = size(X);


%%%%%%%%%%%%%%%%
% MINIMIZATION %
%%%%%%%%%%%%%%%%
if silent == 0
    options.Display = 'iter';
else
    options.Display = 'off';
end
options.MaxIter = 1000;

if use_time == 0
    if NL_type == 0
        [bestk] = minFunc(@(K) LLelog_sinphase(K,X, Robs,stim_params,reg_params),K0,options);
        [LL,grad] = LLelog_sinphase(bestk,X,Robs,stim_params,reg_params);
    else
        [bestk] = minFunc(@(K) LLexp_sinphase(K,X, Robs,stim_params,reg_params),K0,options);
        [LL,grad] = LLexp_sinphase(bestk,X,Robs,stim_params,reg_params);
    end
else
    if NL_type == 0
        [bestk] = minFunc(@(K) LLelog_sinphase_timemod(K,X,Robs,stim_params,reg_params),K0,options);
        [LL,grad] = LLelog_sinphase_timemod(bestk,X,Robs,stim_params,reg_params);
    else
        [bestk] = minFunc(@(K) LLexp_sinphase_timemod(K,X,Robs,stim_params,reg_params),K0,options);
        [LL,grad] = LLexp_sinphase_timemod(bestk,X,Robs,stim_params,reg_params);
    end
end


%% Calculate final magnitude of gradient as evidence of why fit stopped
    % grad = sqrt(sum(grad.^2));
    
    %% Make structure with final results
    fitp.k = bestk;
    % fitp.LP = LL;
    % fitp.LL = LL - LLadjust(bestk,lamrange,length(spksN),lamrange2,llist);
