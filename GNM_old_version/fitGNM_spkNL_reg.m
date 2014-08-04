function gnm = fitGNM_spkNL_reg(gnm,X,spkbs,n_pts_perdim)

if nargin < 4
    n_pts_perdim = 50;
end

options.Display = 'none';

[stimlen] = size(X,1);

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
%INITIAL PARAMETERS, assuming spk NL of form: F[x] = alpha*log(1+exp(beta*(x-theta))
init_params(1) = gnm.spk_alpha;
init_params(2) = gnm.spk_beta;
init_params(3) = gnm.spk_theta;

%INITIALIZE CONSTRAINTS
LB = [1e-3 1e-3 -1e3];
UB = [1e3 1e3 1e3];
Aeq = [];
Beq = [];

kmat = get_k_mat(gnm);
kern_outs = X*kmat;
nmods = size(kern_outs,2);
for i = 1:nmods
    kern_outs(kern_outs(:,i) < 0,i) = 0;
    kern_outs(:,i) = kern_outs(:,i)*gnm.mods(i).w;
end
g = sum(kern_outs,2);

%FIT UNCONSTRAINED PARAMETERS
fit_params = fmincon(@(K) spkNL_internal_LL(K,g,Robs), init_params,[],[],Aeq,Beq,LB,UB,[],options);

%% NOW CONSTRUCT LOOKUP TABLE OF APPROXIMATE SCALE CONVERSIONS FOR ALPHA AND BETA
%convert scale conversions into penalty terms
%NOTE: these penalty terms already are divided by nspks
[nll, pnll, lpen] = getLL_GNM(gnm,X,spkbs,'none');
l2_pens = lpen.l2x + lpen.grad + lpen.lapl + lpen.laplXT;
l1_pens = lpen.l1x;
if max(l2_pens) + max(l1_pens) > 0
    
    search_fac = 2; %search over range from 1/search_fac*init_val to search_fac*init_val
    [alpha_t,beta_t] = meshgrid(linspace(1/search_fac*fit_params(1),search_fac*fit_params(1),n_pts_perdim),...
        linspace(1/search_fac*fit_params(2),search_fac*fit_params(2),n_pts_perdim));
    lut = nan(n_pts_perdim,n_pts_perdim,nmods);
    for i = 1:n_pts_perdim
        for j = 1:n_pts_perdim
            internal_1 = beta_t(i,j)*g;
            internal_2 = alpha_t(i,j)*log(1+exp(internal_1));
            internal_2(internal_1 > 50) = alpha_t(i,j)*beta_t(i,j)*g(internal_1 > 50);
            response = log(-1 + exp(internal_2));
            response(internal_2 > 50) = internal_2(internal_2 > 50);
            %         response = log(-1 + exp(alpha_t(i,j)*(log(1+exp(beta_t(i,j)*g)))));
            %         lut(i,j,:) = regress(response,kern_outs);
            lut(i,j,:) = kern_outs\response;
        end
    end
    
    l1_pens = l1_pens'; l2_pens = l2_pens';
    l1_pens = shiftdim(l1_pens,-1);
    l2_pens = shiftdim(l2_pens,-1);
    lut_l1 = bsxfun(@times,lut,l1_pens);
    lut_l2 = bsxfun(@times,lut.^2,l2_pens);
    
    lut_pen = lut_l1 + lut_l2;
    %%
    options.Display = 'iter';
    LB = [1/search_fac*fit_params(1) 1/search_fac*fit_params(2) -1e3];
    UB = [search_fac*fit_params(1) search_fac*fit_params(2) 1e3];
    fit_params_reg = fmincon(@(K) spkNL_internal_LL_reg(K,g,Robs,lut_pen,alpha_t,beta_t),...
        fit_params,[],[],Aeq,Beq,LB,UB,[],options);
    
    if ~any(isnan(fit_params_reg))
        gnm.spk_alpha = fit_params_reg(1);
        gnm.spk_beta = fit_params_reg(2);
        gnm.spk_theta = fit_params_reg(3);
    end
else
    gnm.spk_alpha = fit_params(1);
    gnm.spk_beta = fit_params(2);
    gnm.spk_theta = fit_params(3);
end
% %now rescale filters
% for i = 1:nmods
%     resc_value = interp2(alpha_t,beta_t,squeeze(lut(:,:,i)),fit_params_reg(1),fit_params_reg(2));
%     gnm.mods(i).k = gnm.mods(i).k*resc_value;
% end

end

%% internal LL evaluation
function [LL grad] = spkNL_internal_LL(params,g,Robs)

n_spks = sum(Robs);

internal = params(2)*(g + params(3));
too_big = find(internal > 50);
lexp = log(1+exp(internal));
lexp(too_big) = params(2)*(g(too_big)+params(3));

r = params(1)*lexp;

r(r < 1e-10) = 1e-10;
LL = -sum(Robs'.*log(r)-r)/n_spks;

fract = exp(params(2)*(g+params(3)))./(1 + exp(params(2)*(g+params(3))));
fract(too_big) = 1;

% grad(1) = 1/params(2)*sum(lexp);
% grad(2) = -params(1)/params(2)^2*sum(lexp) + params(1)/params(2)*sum(fract.*(g-params(3)));
% grad(3) = -params(1)*sum(fract);
grad(1) = sum(lexp);
grad(2) = params(1)*sum(fract.*(g+params(3)));
grad(3) = params(1)*params(2)*sum(fract);

grad = -grad/n_spks;

end

%% internal LL evaluation with lookup table regularization
function LL = spkNL_internal_LL_reg(params,g,Robs,lut_pen,lut_alphas,lut_betas)

n_spks = sum(Robs);

internal = params(2)*(g + params(3));
too_big = find(internal > 50);
lexp = log(1+exp(internal));
lexp(too_big) = params(2)*(g(too_big)+params(3));

r = params(1)*lexp;

r(r < 1e-10) = 1e-10;
LL = -sum(Robs'.*log(r)-r)/n_spks;

nmods = size(lut_pen,3);
LP = 0;
for i = 1:nmods
    cur_pen = interp2(lut_alphas,lut_betas,squeeze(lut_pen(:,:,i)),params(1),params(2));
    LP = LP + cur_pen;
end

LL = LL + LP;

end
