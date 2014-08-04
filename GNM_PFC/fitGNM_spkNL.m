function gnm = fitGNM_spkNL(gnm,g,spkbs,display,hold_const)

%Fit 3 parameters of a alpha/beta*log(1+exp(beta*(x-theta))) model for the spiking nonlinearity.
% INPUTS:
%     gnm: Model structure containing the current estimates of the spiking nonlinearity parameters (spk_alpha, spk_beta, and spk_theta)
%     g: vector of generating function values
%     spkbs: vector of spike bin indices
%     display: (0 or 1) to determine whether to monitor convergence and display plot of fit
%     hold_const: vector specifying which of the three parameters to hold constant [alpha beta theta]
%
% OUTPUTS:
%     gnm: updated model structure with fit spk NL parameters

if nargin < 4
    display = 0;
end
if nargin < 5
    hold_const = [];
end

stimlen = length(g);

%compute binned spike vector
Robs = zeros(1,stimlen);
ftable = tabulate(spkbs);
Robs(ftable(:,1)) = ftable(:,2);

%create default values if these fields aren't already present
if ~isfield(gnm,'spk_alpha')
    gnm.spk_alpha = 1;
end
if ~isfield(gnm,'spk_beta')
    gnm.spk_beta = 1;
end
if ~isfield(gnm,'spk_theta')
    gnm.spk_theta = 0;
end
%INITIAL PARAMETERS
init_params(1) = gnm.spk_alpha;
init_params(2) = gnm.spk_beta;
init_params(3) = 0;

orig_theta = gnm.spk_theta;

%INITIALIZE CONSTRAINTS
LB = [1e-3 1e-3 -1e3];
UB = [1e3 1e3 1e3];
Aeq = [];
Beq = [];
for i = 1:length(hold_const)
    Aeq = [Aeq; zeros(1,3)];
    Aeq(end,hold_const(i)) = 1;
    Beq = [Beq; init_params(i)];
end

%MODEL FITTING
if display == 1
    options.Display = 'iter';
else
    options.Display = 'none';
end
fit_params = fmincon(@(K) spkNL_internal_LL(K,g,Robs), init_params,[],[],Aeq,Beq,LB,UB,[],options);

%store fit values
gnm.spk_alpha = fit_params(1);
gnm.spk_beta = fit_params(2);
gnm.spk_theta = orig_theta + fit_params(3);

%% FOR PLOT OF FIT
if display    
    dist_bins = 200;
    nonpar_bins = 100;
    bin_edges = prctile(g,linspace(0.1,99.9,nonpar_bins));
    bin_centers = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
    [n,x] = hist(g,dist_bins);
    nonpar = nan(nonpar_bins-1,1);
    for i = 1:nonpar_bins-1
        cur_set = find(g >= bin_edges(i) & g < bin_edges(i+1));
        if ~isempty(cur_set)
            nonpar(i) = mean(Robs(cur_set));
        end
    end
    
    pred = fit_params(1)*log(1 + exp(fit_params(2)*(bin_centers + fit_params(3))));
    init = init_params(1)*log(1 + exp(init_params(2)*(bin_centers + init_params(3))));
    figure
    plot(bin_centers,nonpar,'.-');
    hold on
    plot(bin_centers,pred,'r.-')
    plot(bin_centers,init,'k.-')
    hold on
    yl = ylim();
    plot(x,n/max(n)*yl(2),'g')
    legend('Nonparametric','Predicted','Initial','Gen Dist')
end


end

%% internal LL evaluation
function [LL grad] = spkNL_internal_LL(params,g,Robs)

n_spks = sum(Robs);

internal = params(2)*(g + params(3));
too_big = find(internal > 50);
lexp = log(1+exp(internal));
lexp(too_big) = params(2)*(g(too_big)+params(3));

r = params(1)*lexp;

r(r < 1e-50) = 1e-50;
LL = -sum(Robs'.*log(r)-r)/n_spks;

fract = exp(params(2)*(g+params(3)))./(1 + exp(params(2)*(g+params(3))));
fract(too_big) = 1;

% grad(1) = 1/params(2)*sum(lexp);
% grad(2) = -params(1)/params(2)^2*sum(lexp) + params(1)/params(2)*sum(fract.*(g+params(3)));
% grad(3) = -params(1)*sum(fract);
grad(1) = sum(lexp);
grad(2) = params(1)*sum(fract.*(g+params(3)));
grad(3) = params(1)*params(2)*sum(fract);

grad = -grad/n_spks;

end