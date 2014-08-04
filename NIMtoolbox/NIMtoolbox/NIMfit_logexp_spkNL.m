function nim_out = NIMfit_logexp_spkNL(nim,Robs,Xstim,XLin,display,hold_const)
%
% nim_out = NIMfit_logexp_spkNL(nim,Robs,Xstim,<XLin>,<display>,<hold_const>)
%
% Fit 3 parameters of a alpha/beta*log(1+exp(beta*(x-theta))) model for the spiking nonlinearity.
% INPUTS:
%     nim: Model structure containing the current estimates of the spiking nonlinearity parameters (spk_alpha, spk_beta, and spk_theta)
%     Robs: vector of binned spike counts
%     Xstim: time-embedded stimulus mat
%     <Xlin>: Matrix specifying additional linear predictors
%     <isplay>: (0 or 1) to determine whether to monitor convergence and display plot of fit
%     <hold_const>: vector specifying which of the three parameters to hold constant [alpha beta theta]
%
% OUTPUTS:
%     nim_out: updated model structure with fit spk NL parameters

%%
if nargin < 4
    XLin = [];
end
if nargin < 5
    display = 0;
end
if nargin < 6
    hold_const = [];
end

%make sure Robs is a column vector
if size(Robs,2) > size(Robs,1)
    Robs = Robs';
end

NT = size(Xstim,1);

initial_params = nim.spk_NL_params;

%INITIALIZE CONSTRAINTS (keep theta from getting too big, and alpha and
%beta are non-negative)
LB = [-1e3 1e-3 1e-3];
UB = [1e3 1e3 1e3];
Aeq = [];
Beq = [];
for i = 1:length(hold_const)
    Aeq = [Aeq; zeros(1,3)];
    Aeq(end,hold_const(i)) = 1;
    Beq = [Beq; initial_params(hold_const(i))];
end

[~, ~, ~, G] = NIMmodel_eval(nim,Robs,Xstim,XLin);
G = G - nim.spk_NL_params(1);

opts.Display = 'off';
opts.GradObj = 'off';
opts.Algorithm = 'active-set';
fit_params = fmincon(@(K) spkNL_internal_LL(K,G,Robs), initial_params,[],[],Aeq,Beq,LB,UB,[],opts);

nim_out = nim;
nim_out.spk_NL_params = fit_params;

[LL, penLL] = NIMmodel_eval(nim_out,Robs,Xstim,XLin);
nim_out.LL_seq = cat(1,nim_out.LL_seq,LL);
nim_out.penLL_seq = cat(1,nim_out.penLL_seq,penLL);
nim_out.opt_history = cat(1,nim_out.opt_history,{'spkNL'});

%% FOR PLOT OF FIT
if display    
    dist_bins = 200;
    nonpar_bins = 100;
    bin_edges = my_prctile(G,linspace(0.05,99.95,nonpar_bins));
    bin_centers = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
    [n,x] = hist(G,dist_bins);
    nonpar = nan(nonpar_bins-1,1);
    for i = 1:nonpar_bins-1
        cur_set = find(G >= bin_edges(i) & G < bin_edges(i+1));
        if ~isempty(cur_set)
            nonpar(i) = mean(Robs(cur_set));
        end
    end
    
    pred = fit_params(3)*log(1 + exp(fit_params(2)*(bin_centers + fit_params(1))));
    init = initial_params(3)*log(1 + exp(initial_params(2)*(bin_centers + initial_params(1))));
    figure
    plot(bin_centers,nonpar/nim.stim_params.dt,'.-');
    hold on
    plot(bin_centers,pred/nim.stim_params.dt,'r-')
    plot(bin_centers,init/nim.stim_params.dt,'-','color',[0.2 0.8 0.2])
    hold on
    yl = ylim();
    plot(x,n/max(n)*yl(2),'k')
    legend('Nonparametric','Fit','Initial','Gen Dist')
    xlim(bin_edges([1 end]));
    xlabel('Generating signal','fontsize',12)
    ylabel('Firing rate','fontsize',12)
end


end

%% internal LL evaluation
%using the version without user-supplied gradient but can switch back to
%using the trust-region algo with gradient by uncommenting code below

function [LL grad] = spkNL_internal_LL(params,G,Robs)
% function [LL] = spkNL_internal_LL(params,G,Robs)

nspks = sum(Robs);

internal = params(2)*(G + params(1));
too_big = find(internal > 50);
lexp = log(1+exp(internal));
lexp(too_big) = params(2)*(G(too_big)+params(1));

r = params(3)*lexp;

r(r < 1e-50) = 1e-50;
LL = -sum(Robs.*log(r)-r)/nspks;

% %for gradient calculation: (faster to not supply gradient it seems)
% residual = (Robs./r - 1);
% 
% fract = exp(params(2)*(G+params(1)))./(1 + exp(params(2)*(G+params(1))));
% fract(too_big) = 1;
% 
% grad(1) = params(3)*params(2)*residual'*fract;
% grad(2) = params(3)*residual'*(fract.*(G+params(1)));
% grad(3) = residual'*r/params(3);
% 
% grad = -grad/nspks;

end