function [LL,LLgrad] = LLgauss_phasemod(K,X,Yobs,stim_params,reg_params)

[NT,NPAR] = size(X);

k = K(1:end-1);
b = K(end);

% kx = X*k+Xc+b;
y = X*k+b;

% LL = sum(Robs.* log(r) - r);
LL = -sum((Yobs-y).^2);

% residual = (Robs./r - 1) .* ekx ./ (1+ekx);
% residual(too_large) = (Robs(too_large)./r(too_large) - 1);
residual = -2*(y-Yobs);

LLgrad = zeros(NPAR+1,1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);
LLgrad(1:end-1) = residual'*X;

%%
filt_klen = stim_params(1)*stim_params(2);
K_ind = reshape(K(1:filt_klen),stim_params);
if NPAR > filt_klen + 1
    include_dep = 1;
    K_dep = reshape(K((filt_klen+1):(2*filt_klen)),stim_params);
else
    include_dep = 0;
end

%% COMPUTE PHASE D1 PENS
% lapl_oper = [-1; 2; -1];

K_ind_derivs = nan(size(K_ind));
K_ind_derivs(2:end,:) = K_ind(2:end,:) - K_ind(1:end-1,:);
if reg_params.is_phase == 1
    K_ind_derivs(1,:) = K_ind(1,:) - K_ind(end,:);
else
    K_ind_derivs(1,:) = 0;
end
K_ind_pen = reg_params.dl1_ind*sum(K_ind_derivs(:).^2);

K_ind_lapl = nan(size(K_ind));
K_ind_lapl(2:end-1,:) = 2*K_ind(2:end-1,:) - K_ind(1:end-2,:)-K_ind(3:end,:);
if reg_params.is_phase==1
    K_ind_lapl(1,:) = 2*K_ind(1,:) - K_ind(end,:)-K_ind(2,:);
    K_ind_lapl(end,:) = 2*K_ind(end,:) - K_ind(1,:) - K_ind(end-1,:);
else
    K_ind_lapl(1,:) = 2*K_ind(1,:)-2*K_ind(2,:);
    K_ind_lapl(end,:) = 2*K_ind(end,:)-2*K_ind(end-1,:);
end

LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.dl1_ind*K_ind_lapl(:);

if include_dep
    K_dep_derivs = nan(size(K_dep));
    K_dep_derivs(2:end,:) = K_dep(2:end,:) - K_dep(1:end-1,:);
    if reg_params.is_phase == 1
        K_dep_derivs(1,:) = K_dep(1,:) - K_dep(end,:);
    else
        K_dep_derivs(1,:) = 0;
    end
    K_dep_pen = reg_params.dl1_dep*sum(K_dep_derivs(:).^2);
    
    K_dep_lapl = nan(size(K_dep));
    K_dep_lapl(2:end-1,:) = 2*K_dep(2:end-1,:) - K_dep(1:end-2,:)-K_dep(3:end,:);
    if reg_params.is_phase == 1
        K_dep_lapl(1,:) = 2*K_dep(1,:) - K_dep(end,:)-K_dep(2,:);
        K_dep_lapl(end,:) = 2*K_dep(end,:) - K_dep(1,:) - K_dep(end-1,:);
    else
        K_dep_lapl(1,:) = 2*K_dep(1,:) - 2*K_dep(2,:);
        K_dep_lapl(end,:) = 2*K_dep(end,:) - 2*K_dep(end-1,:);
    end
    LLgrad(filt_klen+1:2*filt_klen) = LLgrad((filt_klen+1):2*filt_klen) - 2*reg_params.dl1_dep*K_dep_lapl(:);
else
    K_dep_pen = 0;
end

%% COMPUTE PHASE D2 PENS
K_ind_d2_pen = reg_params.dl2_ind*sum(K_ind_lapl(:).^2);

K_ind_d2grad = zeros(size(K_ind));
K_ind_d2grad(2:end-1,:) = 2*K_ind_lapl(2:end-1,:) - K_ind_lapl(1:end-2,:)-K_ind_lapl(3:end,:);
if reg_params.is_phase==1
    K_ind_d2grad(1,:) = 2*K_ind_lapl(1,:) - K_ind_lapl(end,:)-K_ind_lapl(2,:);
    K_ind_d2grad(end,:) = 2*K_ind_lapl(end,:) - K_ind_lapl(1,:) - K_ind_lapl(end-1,:);
else
    K_ind_d2grad(1,:) = 2*K_ind_lapl(1,:)-2*K_ind_lapl(2,:);
    K_ind_d2grad(end,:) = 2*K_ind_lapl(end,:)-2*K_ind_lapl(end-1,:);
end
LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.dl2_ind*K_ind_d2grad(:);

if include_dep
    K_dep_d2_pen = reg_params.dl2_dep*sum(K_dep_lapl(:).^2);
    
    K_dep_d2grad = zeros(size(K_ind));
    K_dep_d2grad(2:end-1,:) = 2*K_dep_lapl(2:end-1,:) - K_dep_lapl(1:end-2,:)-K_dep_lapl(3:end,:);
    if reg_params.is_phase == 1
        K_dep_d2grad(1,:) = 2*K_dep_lapl(1,:) - K_dep_lapl(end,:)-K_dep_lapl(2,:);
        K_dep_d2grad(end,:) = 2*K_dep_lapl(end,:) - K_dep_lapl(1,:) - K_dep_lapl(end-1,:);
    else
        K_dep_d2grad(1,:) = 2*K_dep_lapl(1,:) - 2*K_dep_lapl(2,:);
        K_dep_d2grad(end,:) = 2*K_dep_lapl(end,:) - 2*K_dep_lapl(end-1,:);
    end
    LLgrad(filt_klen+1:2*filt_klen) = LLgrad((filt_klen+1):2*filt_klen) - 2*reg_params.dl2_dep*K_dep_lapl(:);
else
    K_dep_d2_pen = 0;
end


%% COMPUTE FREQ D2 PENS
if stim_params(2) > 1
    K_find_lapl = nan(size(K_ind));
    K_find_lapl(:,2:end-1) = 2*K_ind(:,2:end-1)-K_ind(:,1:end-2)-K_ind(:,3:end);
    K_find_lapl(:,1) = 2*K_ind(:,1) - 2*K_ind(:,2);
    K_find_lapl(:,end) = 2*K_ind(:,end)-2*K_ind(:,end-1);
    K_find_d2_pen = reg_params.dl2_freq_ind*sum(K_find_lapl(:).^2);
    
    K_find_d2grad = zeros(size(K_ind));
    K_find_d2grad(:,2:end-1) = 2*K_find_lapl(:,2:end-1)-K_find_lapl(:,1:end-2)-K_find_lapl(:,3:end);
    K_find_d2grad(:,1) = 2*K_find_lapl(:,1)-2*K_find_lapl(:,2);
    K_find_d2grad(:,end) = 2*K_find_lapl(:,end)-2*K_find_lapl(:,end-1);
    
    LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.dl2_freq_ind*K_find_d2grad(:);
    
    if include_dep
        K_fdep_lapl = nan(size(K_dep));
        K_fdep_lapl(:,2:end-1) = 2*K_dep(:,2:end-1)-K_dep(:,1:end-2)-K_dep(:,3:end);
        K_fdep_lapl(:,1) = 2*K_dep(:,1) - 2*K_dep(:,2);
        K_fdep_lapl(:,end) = 2*K_dep(:,end)-2*K_dep(:,end-1);
        K_fdep_d2_pen = reg_params.dl2_freq_dep*sum(K_fdep_lapl(:).^2);
        
        K_fdep_d2grad = zeros(size(K_dep));
        K_fdep_d2grad(:,2:end-1) = 2*K_fdep_lapl(:,2:end-1)-K_fdep_lapl(:,1:end-2)-K_fdep_lapl(:,3:end);
        K_fdep_d2grad(:,1) = 2*K_fdep_lapl(:,1)-2*K_fdep_lapl(:,2);
        K_fdep_d2grad(:,end) = 2*K_fdep_lapl(:,end)-2*K_fdep_lapl(:,end-1);
        
        LLgrad(filt_klen+1:2*filt_klen) = LLgrad((filt_klen+1):2*filt_klen) - 2*reg_params.dl2_freq_dep*K_fdep_d2grad(:);
        
    else
        K_fdep_d2_pen = 0;
    end
else
    K_find_d2_pen = 0;
    K_fdep_d2_pen = 0;
end
%%
LL = LL - K_ind_pen - K_dep_pen - K_ind_d2_pen - K_dep_d2_pen - K_find_d2_pen - K_fdep_d2_pen;

dlen = length(Yobs);
LL = -LL/dlen;
LLgrad = -LLgrad/dlen;

