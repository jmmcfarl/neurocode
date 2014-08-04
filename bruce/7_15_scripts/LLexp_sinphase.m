function [LL,LLgrad] = LLexp_sinphase(K,X,Robs,stim_params,reg_params)

[NT,NPAR] = size(X);

k = K(1:end-1);
b = K(end);

% kx = X*k+Xc+b;
kx = X*k+b;
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

%%
filt_klen = stim_params(1)*stim_params(2);
K_sin = reshape(K(1:filt_klen),stim_params(1:2));
K_cos = reshape(K((filt_klen+1):(2*filt_klen)),stim_params(1:2));

%% COMPUTE D2/DF2 PENS
K_sin_lapl = zeros(size(K_sin));
K_sin_lapl(2:end-1,:) = 2*K_sin(2:end-1,:) - K_sin(1:end-2,:)-K_sin(3:end,:);

% K_sin_lapl(1,:) = 2*K_sin(1,:)-2*K_sin(2,:);
% K_sin_lapl(end,:) = 2*K_sin(end,:)-2*K_sin(end-1,:);

K_cos_lapl = zeros(size(K_sin));
K_cos_lapl(2:end-1,:) = 2*K_cos(2:end-1,:) - K_cos(1:end-2,:)-K_cos(3:end,:);
% K_cos_lapl(1,:) = 2*K_cos(1,:)-2*K_cos(2,:);
% K_cos_lapl(end,:) = 2*K_cos(end,:)-2*K_cos(end-1,:);
% K_cos_lapl([1 end],:) = 0;


K_sin_d2_pen = reg_params.dl2_freq*sum(K_sin_lapl(:).^2);
K_cos_d2_pen = reg_params.dl2_freq*sum(K_cos_lapl(:).^2);

K_sin_d2grad = zeros(size(K_sin));
K_sin_d2grad(2:end-1,:) = 2*K_sin_lapl(2:end-1,:) - K_sin_lapl(1:end-2,:)-K_sin_lapl(3:end,:);
% K_sin_d2grad(1,:) = 2*K_sin_lapl(1,:)-2*K_sin_lapl(2,:);
% K_sin_d2grad(end,:) = 2*K_sin_lapl(end,:)-2*K_sin_lapl(end-1,:);
% K_sin_d2grad([1 end],:) = 0;
K_sin_d2grad(1,:) = -K_sin_lapl(2,:);
K_sin_d2grad(end,:) = -K_sin_lapl(end-1,:);

K_cos_d2grad = zeros(size(K_cos));
K_cos_d2grad(2:end-1,:) = 2*K_cos_lapl(2:end-1,:) - K_cos_lapl(1:end-2,:)-K_cos_lapl(3:end,:);
% K_cos_d2grad(1,:) = 2*K_cos_d2grad(1,:) - 2*K_cos_d2grad(2,:);
% K_cos_d2grad(end,:) = 2*K_cos_d2grad(end,:) - 2*K_cos_d2grad(end-1,:);
% K_cos_d2grad([1 end],:) = 0;
K_cos_d2grad(1,:) = -K_cos_lapl(2,:);
K_cos_d2grad(end,:) = -K_cos_lapl(end-1,:);

LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.dl2_freq*K_sin_d2grad(:);   
LLgrad(filt_klen+1:2*filt_klen) = LLgrad((filt_klen+1):2*filt_klen) - 2*reg_params.dl2_freq*K_cos_d2grad(:);

%% COMPUTE D/DF PENS
K_sin_deriv = zeros(size(K_sin));
K_sin_deriv(2:end,:) = K_sin(2:end,:) - K_sin(1:end-1,:);

K_cos_deriv = zeros(size(K_cos));
K_cos_deriv(2:end,:) = K_cos(2:end,:) - K_cos(1:end-1,:);

K_sin_d_pen = reg_params.dl_freq*sum(K_sin_deriv(:).^2);
K_cos_d_pen = reg_params.dl_freq*sum(K_cos_deriv(:).^2);

K_sin_dgrad = zeros(size(K_sin));
K_sin_dgrad(2:end-1,:) = 2*K_sin(2:end-1,:) - K_sin(1:end-2,:) - K_sin(3:end,:);
K_sin_dgrad(1,:) = K_sin(1,:) - K_sin(2,:);
K_sin_dgrad(end,:) = K_sin(end,:) - K_sin(end-1,:);

K_cos_dgrad = zeros(size(K_cos));
K_cos_dgrad(2:end-1,:) = 2*K_cos(2:end-1,:) - K_cos(1:end-2,:) - K_cos(3:end,:);
K_cos_dgrad(1,:) = K_cos(1,:) - K_cos(2,:);
K_cos_dgrad(end,:) = K_cos(end,:) - K_cos(end-1,:);

LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.dl_freq*K_sin_dgrad(:);   
LLgrad(filt_klen+1:2*filt_klen) = LLgrad((filt_klen+1):2*filt_klen) - 2*reg_params.dl_freq*K_cos_dgrad(:);

%% COMPUTE CH D2 PENS
if stim_params(2) > 1
    K_fsin_lapl = zeros(size(K_sin));
    K_fsin_lapl(:,2:end-1) = 2*K_sin(:,2:end-1)-K_sin(:,1:end-2)-K_sin(:,3:end);
    %     K_fsin_lapl(:,1) = 2*K_sin(:,1) - 2*K_sin(:,2);
    %     K_fsin_lapl(:,end) = 2*K_sin(:,end)-2*K_sin(:,end-1);
%     K_fsin_lapl(:,[1 end]) = 0;
    K_fsin_d2_pen = reg_params.dl2_ch*sum(K_fsin_lapl(:).^2);
    
    K_fsin_d2grad = zeros(size(K_sin));
    K_fsin_d2grad(:,2:end-1) = 2*K_fsin_lapl(:,2:end-1)-K_fsin_lapl(:,1:end-2)-K_fsin_lapl(:,3:end);
    %     K_fsin_d2grad(:,1) = 2*K_fsin_lapl(:,1)-2*K_fsin_lapl(:,2);
    %     K_fsin_d2grad(:,end) = 2*K_fsin_lapl(:,end)-2*K_fsin_lapl(:,end-1);
%     K_fsin_d2grad(:,[1 end]) = 0;
    K_fsin_d2grad(:,1) = -K_fsin_lapl(:,2);
    K_fsin_d2grad(:,end) = - K_fsin_lapl(:,end-1);
    LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.dl2_ch*K_fsin_d2grad(:);
    
    K_fcos_lapl = zeros(size(K_cos));
    K_fcos_lapl(:,2:end-1) = 2*K_cos(:,2:end-1)-K_cos(:,1:end-2)-K_cos(:,3:end);
%     K_fcos_lapl(:,1) = 2*K_cos(:,1) - 2*K_cos(:,2);
%     K_fcos_lapl(:,end) = 2*K_cos(:,end)-2*K_cos(:,end-1);
%     K_fcos_lapl(:,[1 end]) = 0;
    K_fcos_d2_pen = reg_params.dl2_ch*sum(K_fcos_lapl(:).^2);
    
    K_fcos_d2grad = zeros(size(K_cos));
    K_fcos_d2grad(:,2:end-1) = 2*K_fcos_lapl(:,2:end-1)-K_fcos_lapl(:,1:end-2)-K_fcos_lapl(:,3:end);
%     K_fcos_d2grad(:,1) = 2*K_fcos_lapl(:,1)-2*K_fcos_lapl(:,2);
%     K_fcos_d2grad(:,end) = 2*K_fcos_lapl(:,end)-2*K_fcos_lapl(:,end-1);
%     K_fcos_d2grad(:,[1 end]) = 0;
    K_fcos_d2grad(:,1) = -K_fcos_lapl(:,2);
    K_fcos_d2grad(:,end) = -K_fcos_lapl(:,end-1);

    LLgrad(filt_klen+1:2*filt_klen) = LLgrad((filt_klen+1):2*filt_klen) - 2*reg_params.dl2_ch*K_fcos_d2grad(:);
    
else
    K_fsin_d2_pen = 0; K_fcos_d2_pen = 0;
end

%% CH D PENS
if stim_params(2) > 1
    
    K_fsin_deriv = zeros(size(K_sin));
    K_fsin_deriv(:,2:end) = K_sin(:,2:end) - K_sin(:,1:end-1);
        
    K_fcos_deriv = zeros(size(K_cos));
    K_fcos_deriv(:,2:end) = K_cos(:,2:end) - K_cos(:,1:end-1);
    
    K_fsin_d_pen = reg_params.dl_ch*sum(K_fsin_deriv(:).^2);
    K_fcos_d_pen = reg_params.dl_ch*sum(K_fcos_deriv(:).^2);
    
    K_fsin_dgrad = zeros(size(K_sin));
K_fsin_dgrad(:,2:end-1) = 2*K_sin(:,2:end-1) - K_sin(:,1:end-2) - K_sin(:,3:end);
K_fsin_dgrad(:,1) = K_sin(:,1) - K_sin(:,2);
K_fsin_dgrad(:,end) = K_sin(:,end) - K_sin(:,end-1);
    
    K_fcos_dgrad = zeros(size(K_cos));
K_fcos_dgrad(:,2:end-1) = 2*K_cos(:,2:end-1) - K_cos(:,1:end-2) - K_cos(:,3:end);
K_fcos_dgrad(:,1) = K_cos(:,1) - K_cos(:,2);
K_fcos_dgrad(:,end) = K_cos(:,end) - K_cos(:,end-1);
    
    LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.dl_ch*K_fsin_dgrad(:);
    LLgrad(filt_klen+1:2*filt_klen) = LLgrad((filt_klen+1):2*filt_klen) - 2*reg_params.dl_ch*K_fcos_dgrad(:);
    
else
    K_fsin_d_pen = 0; K_fcos_d_pen = 0;
end


%%
K2_sin_pen = reg_params.d2_phase*sum(K_sin(:).^2);
K2_cos_pen = reg_params.d2_phase*sum(K_cos(:).^2);
LLgrad(1:filt_klen) = LLgrad(1:filt_klen) - 2*reg_params.d2_phase*K_sin(:);
LLgrad(filt_klen+1:2*filt_klen) = LLgrad((filt_klen+1):2*filt_klen) - 2*reg_params.d2_phase*K_cos(:);


%%
LL = LL - K_sin_d2_pen - K_cos_d2_pen - K_sin_d_pen - K_cos_d_pen - ...
    K_fsin_d2_pen - K_fcos_d2_pen -  K_fsin_d_pen - K_fcos_d_pen -  K2_sin_pen - K2_cos_pen;

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

