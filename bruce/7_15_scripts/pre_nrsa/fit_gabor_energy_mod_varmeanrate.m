function [gabor_params,LL] = fit_gabor_energy_mod_varmeanrate(XX,YY,init_params,Xmat,Robs,expt_inds,hold_const,LB,UB,priors)

if nargin < 9
    for i = 1:6
        priors(i).type = [];
    end
end
options.Display = 'off';
options.Algorithm = 'sqp';
options.MaxFunEvals = 2000;

K0 = init_params(hold_const == 0);
fixed_params = init_params(hold_const == 1);

un_expts = unique(expt_inds);
n_expts = length(un_expts);
K0 = [K0 zeros(1,n_expts)];
LB = [LB repmat(-Inf,1,n_expts)];
UB = [UB repmat(Inf,1,n_expts)];

if hold_const(end) == 1
    fixed_params = [fixed_params repmat(fixed_params(end),1,n_expts)];
    hold_const = [hold_const ones(1,n_expts)];
else
    hold_const = [hold_const zeros(1,n_expts)];
end

[params LL] = fmincon( @(K) get_pgabor_LL(K, Xmat, Robs,expt_inds,XX,YY,hold_const,fixed_params,priors), ...
    K0,[],[],[],[],LB(hold_const==0),UB(hold_const==0),[],options);
gabor_params = init_params;
gabor_params(hold_const == 0) = params;

end

function LL = get_pgabor_LL(K,Xmat,Robs,expt_inds,XX,YY,hold_const,fixed_params,priors)
    
    params = zeros(size(hold_const));
    params(hold_const==1) = fixed_params;
    params(hold_const==0) = K;

    cur_mask1 = get_pgabor_mask_v2(XX,YY,params(1:6),0);
    cur_mask2 = get_pgabor_mask_v2(XX,YY,params(1:6),pi/2);
    
    mask1_out = Xmat*cur_mask1(:);
    mask2_out = Xmat*cur_mask2(:);

    energy_out = params(7)*sqrt(mask1_out.^2+mask2_out.^2);
    
    g = energy_out + params(expt_inds+7)';
    too_large = find(g > 100);
    r = log(1+exp(g));
    r(too_large) = g(too_large);
    
    r(r < 1e-20) = 1e-20;

    LL = sum(Robs.*log(r)-r);
    
    tot_lprior = 0;
    for i = 1:6
       if hold_const(i) == 0 && ~isempty(priors(i).type)
           if strcmp(priors(i).type,'gauss')
               cur_lprior = log(normpdf(params(i),priors(i).theta(1),priors(i).theta(2)));
           elseif strcmp(priors(i).type,'gam')
               cur_lprior = log(gampdf(params(i),priors(i).theta(1),priors(i).theta(2)));
           else
                error('unsupported prior')
           end
            
           tot_lprior = tot_lprior + cur_lprior;
       end
    end
    LL = LL + tot_lprior;
    
    Nspks = sum(Robs);
    LL = -LL/Nspks;

end

