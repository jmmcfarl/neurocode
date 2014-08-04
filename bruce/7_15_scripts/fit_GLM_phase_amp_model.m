function [fitp] = fit_GLM_phase_amp_model( X, Robs, K0, silent, stim_params, reg_params,use_time)
%

if nargin < 7
    use_time = 0;
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

lambda = [];
lambda = [lambda reg_params.pl1_ind*ones(1,stim_params(1)*stim_params(3))];
lambda = [lambda reg_params.pl1_dep*ones(1,stim_params(1)*stim_params(3))];
lambda = [lambda reg_params.al1_ind*ones(1,stim_params(2)*stim_params(3))];
lambda = [lambda reg_params.al1_dep*ones(1,stim_params(2)*stim_params(3))];
if length(lambda) < NPAR
    lambda = [lambda 0];
end
lambda = [lambda 0]/sum(Robs);

if max(lambda) > 0
    [bestk LL] = L1General2_PSSas(@(K) LLelog_phase_ampmod(K,X,Robs,stim_params,reg_params),K0,lambda',options,0);
else
    [bestk] = minFunc(@(K) LLelog_phase_ampmod(K,X,Robs,stim_params,reg_params),K0,options);
    [LL,grad] = LLelog_phase_ampmod(bestk,X,Robs,stim_params,reg_params);
end


%% Calculate final magnitude of gradient as evidence of why fit stopped
% grad = sqrt(sum(grad.^2));

%% Make structure with final results
fitp.k = bestk;
% fitp.LP = LL;
% fitp.LL = LL - LLadjust(bestk,lamrange,length(spksN),lamrange2,llist);
