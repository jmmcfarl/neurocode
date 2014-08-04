function [fitp] = fit_GLM_phase_model( X,Robs, K0, silent, stim_params, reg_params,use_time,NL_type)
%

if nargin < 7
    use_time = 0;
end
if nargin < 8
    NL_type = 0;
end
% if nargin < 9
%     sinmod = 0;
% end
[NT NPAR] = size(X);


%%%%%%%%%%%%%%%%
% MINIMIZATION %
%%%%%%%%%%%%%%%%
if silent == 0
    options.Display = 'iter';
else
    options.Display = 'off';
end

sinmod = 0;
if use_time==1
%     if stim_params(2)==1 & length(K0) == stim_params(1)*2+stim_params(3)+1
    if stim_params(1)*stim_params(2)*2+stim_params(3)+1 == length(K0) & stim_params(1) >0
        fprintf('using sinmod\n')
        sinmod=1;
    elseif stim_params(1)*stim_params(2)*2+2*stim_params(3)+1 == length(K0) & stim_params(1) >0
        fprintf('using sinmod\n')
        sinmod=1;        
    end
end

lambda = [];
lambda = [lambda reg_params.l1_ind*ones(1,stim_params(1)*stim_params(2))];
if stim_params(end) == 2
    lambda = [lambda reg_params.l1_dep*ones(1,stim_params(1)*stim_params(2))];
end
if use_time == 1
    if reg_params.l1t_ind > 0
        lambda = [lambda reg_params.l1t_ind*ones(1,2*stim_params(3))];
    else
        lambda = [lambda zeros(1,2*stim_params(3))];
    end
    if length(lambda) > NPAR
        lambda = lambda(1:NPAR-1);
    end
end
if length(lambda) < NPAR
    lambda = [lambda 0];
end
lambda = [lambda 0]/sum(Robs);

if max(lambda) > 0
    if use_time == 0
        if NL_type == 0
            [bestk LL] = L1General2_PSSas(@(K) LLelog_phasemod(K,X, Robs,stim_params,reg_params),K0,lambda',options,0);
        else
            [bestk LL] = L1General2_PSSas(@(K) LLexp_phasemod(K,X, Robs,stim_params,reg_params),K0,lambda',options,0);
        end
    else
        if NL_type == 0
            [bestk LL] = L1General2_PSSas(@(K) LLelog_phase_timemod(K,X,Robs,stim_params,reg_params),K0,lambda',options,0);
        else
            [bestk LL] = L1General2_PSSas(@(K) LLexp_phase_timemod(K,X,Robs,stim_params,reg_params),K0,lambda',options,0);
        end
    end
else
    if use_time == 0
        if NL_type == 0
            [bestk] = minFunc(@(K) LLelog_phasemod(K,X, Robs,stim_params,reg_params),K0,options);
            [LL,grad] = LLelog_phasemod(bestk,X,Robs,stim_params,reg_params);
        else
            [bestk] = minFunc(@(K) LLexp_phasemod(K,X, Robs,stim_params,reg_params),K0,options);
            [LL,grad] = LLexp_phasemod(bestk,X,Robs,stim_params,reg_params);
        end
    elseif sinmod==0
        if NL_type == 0
            [bestk] = minFunc(@(K) LLelog_phase_timemod(K,X,Robs,stim_params,reg_params),K0,options);
            [LL,grad] = LLelog_phase_timemod(bestk,X,Robs,stim_params,reg_params);
        else
            [bestk] = minFunc(@(K) LLexp_phase_timemod(K,X,Robs,stim_params,reg_params),K0,options);
            [LL,grad] = LLexp_phase_timemod(bestk,X,Robs,stim_params,reg_params);
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
end
    
    
    %% Calculate final magnitude of gradient as evidence of why fit stopped
    % grad = sqrt(sum(grad.^2));
    
    %% Make structure with final results
    fitp.k = bestk;
    % fitp.LP = LL;
    % fitp.LL = LL - LLadjust(bestk,lamrange,length(spksN),lamrange2,llist);
