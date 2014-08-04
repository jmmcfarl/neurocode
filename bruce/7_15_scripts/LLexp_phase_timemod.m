function [LL,LLgrad] = LLexp_phase_timemod(K,X,Robs,stim_params,reg_params)

[NT,NPAR] = size(X);
filt_klen = stim_params(1)*stim_params(2);
filt_tlen = stim_params(3);

k = K(1:end-1);
b = K(end);

% kx = X*k+Xc+b;
kx = X*k+b;
% kx = X_time*k_time + b;
% kx(phase_inds) = kx(phase_inds) + X_phase*k_phase;

too_large = find(kx > 100);
kx(too_large) = 100;
r = exp(kx);
% r = log(1+ekx);
r(too_large) = kx(too_large);
r(r < 1e-20) = 1e-20;

LL = sum(Robs.* log(r) - r);

% residual = (Robs./r - 1) .* ekx ./ (1+ekx);
% residual(too_large) = (Robs(too_large)./r(too_large) - 1);
residual = (Robs./r - 1) .* r;

LLgrad = zeros(NPAR+1,1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);
LLgrad(1:end-1) = residual'*X;
% LLgrad(phase_set) = residual(phase_inds)'*X_phase;
% LLgrad(time_set) = residual'*X_time;
%%
K_ind = reshape(K(1:filt_klen),stim_params(1:2));
% if NPAR > filt_klen + 1
%     include_dep = 1;
%     K_dep = reshape(K((filt_klen+1):2*filt_klen),stim_params);
% else
%     include_dep = 0;
% end


if stim_params(end) == 2
    Kt_ind = K((2*filt_klen+1):(2*filt_klen+filt_tlen));
    K_dep = reshape(K((filt_klen+1):2*filt_klen),stim_params(1:2));
    if NPAR > 2*filt_klen+filt_tlen+3
        Kt_dep_use = 1;
        Kt_dep = K((2*filt_klen+filt_tlen+1):(2*filt_klen+2*filt_tlen));
    else
        Kt_dep_use = 0;
    end
else
    Kt_ind = K((filt_klen+1):(filt_klen+filt_tlen));
    if NPAR > filt_klen + filt_tlen + 3
        Kt_dep_use = 1;
        Kt_dep = K((filt_klen+filt_tlen+1):(filt_klen+2*filt_tlen));
    else
        Kt_dep_use = 0;
    end
end
%% COMPUTE PHASE D1 PENS
% lapl_oper = [-1; 2; -1];

if stim_params(1) > 0
    K_ind_derivs = nan(size(K_ind));
    K_ind_derivs(2:end,:) = K_ind(2:end,:) - K_ind(1:end-1,:);
    if reg_params.is_phase==1
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
    K_ind_lapl(1,:) = 2*K_ind(1,:) - 2*K_ind(2,:);
    K_ind_lapl(end,:) = 2*K_ind(end,:) - 2*K_ind(end-1,:);        
    end
    LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.dl1_ind*K_ind_lapl(:);
    
    if stim_params(end)==2
        K_dep_derivs = nan(size(K_dep));
        K_dep_derivs(2:end,:) = K_dep(2:end,:) - K_dep(1:end-1,:);
        if reg_params.is_phase==1
        K_dep_derivs(1,:) = K_dep(1,:) - K_dep(end,:);
        else
            K_dep_derivs(1,:) = 0;
        end
        K_dep_pen = reg_params.dl1_ind*sum(K_dep_derivs(:).^2);
        
        K_dep_lapl = nan(size(K_dep));
        K_dep_lapl(2:end-1,:) = 2*K_dep(2:end-1,:) - K_dep(1:end-2,:)-K_dep(3:end,:);
        if reg_params.is_phase==1
        K_dep_lapl(1,:) = 2*K_dep(1,:) - K_dep(end,:)-K_dep(2,:);
        K_dep_lapl(end,:) = 2*K_dep(end,:) - K_dep(1,:) - K_dep(end-1,:);
        else
         K_dep_lapl(1,:) = 2*K_dep(1,:) - 2*K_dep(2,:);
        K_dep_lapl(end,:) = 2*K_dep(end,:) - 2*K_dep(end-1,:);
       end
        LLgrad(filt_klen+1:(2*filt_klen)) = LLgrad((filt_klen+1):(2*filt_klen)) - 2*reg_params.dl1_dep*K_dep_lapl(:);
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
    K_ind_d2grad(1,:) = 2*K_ind_lapl(1,:) - 2*K_ind_lapl(2,:);
    K_ind_d2grad(end,:) = 2*K_ind_lapl(end,:) - 2*K_ind_lapl(end-1,:);
    end
    LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.dl2_ind*K_ind_d2grad(:);
    
    if stim_params(end)==2
        K_dep_d2_pen = reg_params.dl2_dep*sum(K_dep_lapl(:).^2);
        
        K_dep_d2grad = zeros(size(K_ind));
        K_dep_d2grad(2:end-1,:) = 2*K_dep_lapl(2:end-1,:) - K_dep_lapl(1:end-2,:)-K_dep_lapl(3:end,:);
        if reg_params.is_phase==1
        K_dep_d2grad(1,:) = 2*K_dep_lapl(1,:) - K_dep_lapl(end,:)-K_dep_lapl(2,:);
        K_dep_d2grad(end,:) = 2*K_dep_lapl(end,:) - K_dep_lapl(1,:) - K_dep_lapl(end-1,:);
        else
        K_dep_d2grad(1,:) = 2*K_dep_lapl(1,:) - 2*K_dep_lapl(2,:);
        K_dep_d2grad(end,:) = 2*K_dep_lapl(end,:) - 2*K_dep_lapl(end-1,:);
        end
        LLgrad(filt_klen+1:(2*filt_klen)) = LLgrad((filt_klen+1):(2*filt_klen)) - 2*reg_params.dl2_dep*K_dep_lapl(:);
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
        
        if stim_params(end)==2
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
else
    K_ind_pen = 0;
    K_dep_pen = 0;
    K_ind_d2_pen = 0;
    K_dep_d2_pen = 0;
    K_find_d2_pen = 0;
    K_fdep_d2_pen = 0;
    
end

%% COMPUTE TIME D2 PENS
K_tind_lapl = nan(size(Kt_ind));
K_tind_lapl(2:end-1) = 2*Kt_ind(2:end-1) - Kt_ind(1:end-2) - Kt_ind(3:end);
K_tind_lapl(1) = 2*Kt_ind(1)-2*Kt_ind(2);
K_tind_lapl(end) = 2*Kt_ind(end)-2*Kt_ind(end-1);
K_tind_d2_pen = reg_params.dl2_time_ind*sum(K_tind_lapl(:).^2);

K_tind_d2grad = zeros(size(Kt_ind));
K_tind_d2grad(2:end-1) = 2*K_tind_lapl(2:end-1) - K_tind_lapl(1:end-2) - K_tind_lapl(3:end);
K_tind_d2grad(1) = 2*K_tind_lapl(1) - 2*K_tind_lapl(2);
K_tind_d2grad(end) = 2*K_tind_lapl(end) - 2*K_tind_lapl(end-1);

if stim_params(end)==2
    LLgrad((2*filt_klen+1):(2*filt_klen+filt_tlen)) = LLgrad((2*filt_klen+1):(2*filt_klen+filt_tlen)) - ...
        2*reg_params.dl2_time_ind*K_tind_d2grad(:);
else
    LLgrad((filt_klen+1):(filt_klen+filt_tlen)) = LLgrad((filt_klen+1):(filt_klen+filt_tlen)) - ...
        2*reg_params.dl2_time_ind*K_tind_d2grad(:);
end

if Kt_dep_use == 1
    K_tdep_lapl = nan(size(Kt_dep));
    K_tdep_lapl(2:end-1) = 2*Kt_dep(2:end-1) - Kt_dep(1:end-2) - Kt_dep(3:end);
    K_tdep_lapl(1) = 2*Kt_dep(1)-2*Kt_dep(2);
    K_tdep_lapl(end) = 2*Kt_dep(end)-2*Kt_dep(end-1);
    K_tdep_d2_pen = reg_params.dl2_time_dep*sum(K_tdep_lapl(:).^2);
    
    K_tdep_d2grad = zeros(size(Kt_dep));
    K_tdep_d2grad(2:end-1) = 2*K_tdep_lapl(2:end-1) - K_tdep_lapl(1:end-2) - K_tdep_lapl(3:end);
    K_tdep_d2grad(1) = 2*K_tdep_lapl(1) - 2*K_tdep_lapl(2);
    K_tdep_d2grad(end) = 2*K_tdep_lapl(end) - 2*K_tdep_lapl(end-1);
    
    if stim_params(end)==2
        LLgrad((2*filt_klen+filt_tlen+1):(2*filt_klen+2*filt_tlen)) = LLgrad((2*filt_klen+filt_tlen+1):(2*filt_klen+2*filt_tlen)) - ...
            2*reg_params.dl2_time_dep*K_tdep_d2grad(:);
    else
        LLgrad((filt_klen+filt_tlen+1):(filt_klen+2*filt_tlen)) = LLgrad((filt_klen+filt_tlen+1):(filt_klen+2*filt_tlen)) - ...
            2*reg_params.dl2_time_dep*K_tdep_d2grad(:);
    end
else
    K_tdep_d2_pen =  0;
end
%%
LL = LL - K_ind_pen - K_dep_pen - K_ind_d2_pen - K_dep_d2_pen - K_find_d2_pen - ...
    K_fdep_d2_pen - K_tind_d2_pen - K_tdep_d2_pen;

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

