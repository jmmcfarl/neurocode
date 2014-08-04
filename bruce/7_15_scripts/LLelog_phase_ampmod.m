function [LL,LLgrad] = LLelog_phase_ampmod(K,X,Robs,stim_params,reg_params)

[NT,NPAR] = size(X);

k = K(1:end-1);
b = K(end);

kx = X*k+b;
too_large = find(kx > 50);
ekx = exp(kx);
r = log(1+ekx);
r(too_large) = kx(too_large);
r(r < 1e-20) = 1e-20;

LL = sum(Robs.* log(r) - r);

residual = (Robs./r - 1) .* ekx ./ (1+ekx);
residual(too_large) = (Robs(too_large)./r(too_large) - 1);

LLgrad = zeros(NPAR+1,1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);
LLgrad(1:end-1) = residual'*X;

%%
pfilt_klen = stim_params(1)*stim_params(3);
afilt_klen = stim_params(2)*stim_params(3);

pk_ind_inds = 1:pfilt_klen;
pk_dep_inds = (pfilt_klen+1):(2*pfilt_klen);
ak_ind_inds = (2*pfilt_klen+1):(2*pfilt_klen+afilt_klen);
ak_dep_inds = (2*pfilt_klen+afilt_klen+1):(2*pfilt_klen+2*afilt_klen);

pK_ind = reshape(K(pk_ind_inds),stim_params([1 3]));
pK_dep = reshape(K(pk_dep_inds),stim_params([1 3]));
aK_ind = reshape(K(ak_ind_inds),stim_params([2 3]));
aK_dep = reshape(K(ak_dep_inds),stim_params([2 3]));


%% COMPUTE PHASE D1 PENS
pK_ind_derivs = nan(size(pK_ind));
pK_ind_derivs(2:end,:) = pK_ind(2:end,:) - pK_ind(1:end-1,:);
pK_ind_derivs(1,:) = pK_ind(1,:) - pK_ind(end,:);
pK_ind_pen = reg_params.pdl1_ind*sum(pK_ind_derivs(:).^2);

pK_ind_lapl = nan(size(pK_ind));
pK_ind_lapl(2:end-1,:) = 2*pK_ind(2:end-1,:) - pK_ind(1:end-2,:)-pK_ind(3:end,:);
pK_ind_lapl(1,:) = 2*pK_ind(1,:) - pK_ind(end,:)-pK_ind(2,:);
pK_ind_lapl(end,:) = 2*pK_ind(end,:) - pK_ind(1,:) - pK_ind(end-1,:);

LLgrad(pk_ind_inds) = LLgrad(pk_ind_inds) - 2*reg_params.pdl1_ind*pK_ind_lapl(:);

pK_dep_derivs = nan(size(pK_dep));
pK_dep_derivs(2:end,:) = pK_dep(2:end,:) - pK_dep(1:end-1,:);
pK_dep_derivs(1,:) = pK_dep(1,:) - pK_dep(end,:);
pK_dep_pen = reg_params.pdl1_dep*sum(pK_dep_derivs(:).^2);

pK_dep_lapl = nan(size(pK_dep));
pK_dep_lapl(2:end-1,:) = 2*pK_dep(2:end-1,:) - pK_dep(1:end-2,:)-pK_dep(3:end,:);
pK_dep_lapl(1,:) = 2*pK_dep(1,:) - pK_dep(end,:)-pK_dep(2,:);
pK_dep_lapl(end,:) = 2*pK_dep(end,:) - pK_dep(1,:) - pK_dep(end-1,:);
LLgrad(pk_dep_inds) = LLgrad(pk_dep_inds) - 2*reg_params.pdl1_dep*pK_dep_lapl(:);

%% COMPUTE PHASE D2 PENS
pK_ind_d2_pen = reg_params.pdl2_ind*sum(pK_ind_lapl(:).^2);

pK_ind_d2grad = zeros(size(pK_ind));
pK_ind_d2grad(2:end-1,:) = 2*pK_ind_lapl(2:end-1,:) - pK_ind_lapl(1:end-2,:)-pK_ind_lapl(3:end,:);
pK_ind_d2grad(1,:) = 2*pK_ind_lapl(1,:) - pK_ind_lapl(end,:)-pK_ind_lapl(2,:);
pK_ind_d2grad(end,:) = 2*pK_ind_lapl(end,:) - pK_ind_lapl(1,:) - pK_ind_lapl(end-1,:);
LLgrad(pk_ind_inds) = LLgrad(pk_ind_inds) - 2*reg_params.pdl2_ind*pK_ind_d2grad(:);

pK_dep_d2_pen = reg_params.pdl2_dep*sum(pK_dep_lapl(:).^2);

pK_dep_d2grad = zeros(size(pK_ind));
pK_dep_d2grad(2:end-1,:) = 2*pK_dep_lapl(2:end-1,:) - pK_dep_lapl(1:end-2,:)-pK_dep_lapl(3:end,:);
pK_dep_d2grad(1,:) = 2*pK_dep_lapl(1,:) - pK_dep_lapl(end,:)-pK_dep_lapl(2,:);
pK_dep_d2grad(end,:) = 2*pK_dep_lapl(end,:) - pK_dep_lapl(1,:) - pK_dep_lapl(end-1,:);
LLgrad(pk_dep_inds) = LLgrad(pk_dep_inds) - 2*reg_params.pdl2_dep*pK_dep_lapl(:);


%% COMPUTE FREQ D2 PENS
pK_find_lapl = nan(size(pK_ind));
pK_find_lapl(:,2:end-1) = 2*pK_ind(:,2:end-1)-pK_ind(:,1:end-2)-pK_ind(:,3:end);
pK_find_lapl(:,1) = 2*pK_ind(:,1) - 2*pK_ind(:,2);
pK_find_lapl(:,end) = 2*pK_ind(:,end)-2*pK_ind(:,end-1);
pK_find_d2_pen = reg_params.pdl2_freq_ind*sum(pK_find_lapl(:).^2);

pK_find_d2grad = zeros(size(pK_ind));
pK_find_d2grad(:,2:end-1) = 2*pK_find_lapl(:,2:end-1)-pK_find_lapl(:,1:end-2)-pK_find_lapl(:,3:end);
pK_find_d2grad(:,1) = 2*pK_find_lapl(:,1)-2*pK_find_lapl(:,2);
pK_find_d2grad(:,end) = 2*pK_find_lapl(:,end)-2*pK_find_lapl(:,end-1);

LLgrad(pk_ind_inds) = LLgrad(pk_ind_inds) - 2*reg_params.pdl2_freq_ind*pK_find_d2grad(:);

pK_fdep_lapl = nan(size(pK_dep));
pK_fdep_lapl(:,2:end-1) = 2*pK_dep(:,2:end-1)-pK_dep(:,1:end-2)-pK_dep(:,3:end);
pK_fdep_lapl(:,1) = 2*pK_dep(:,1) - 2*pK_dep(:,2);
pK_fdep_lapl(:,end) = 2*pK_dep(:,end)-2*pK_dep(:,end-1);
pK_fdep_d2_pen = reg_params.pdl2_freq_dep*sum(pK_fdep_lapl(:).^2);

pK_fdep_d2grad = zeros(size(pK_dep));
pK_fdep_d2grad(:,2:end-1) = 2*pK_fdep_lapl(:,2:end-1)-pK_fdep_lapl(:,1:end-2)-pK_fdep_lapl(:,3:end);
pK_fdep_d2grad(:,1) = 2*pK_fdep_lapl(:,1)-2*pK_fdep_lapl(:,2);
pK_fdep_d2grad(:,end) = 2*pK_fdep_lapl(:,end)-2*pK_fdep_lapl(:,end-1);

LLgrad(pk_dep_inds) = LLgrad(pk_dep_inds) - 2*reg_params.pdl2_freq_dep*pK_fdep_d2grad(:);
   
%% COMPUTE AMP D1 PENS
aK_ind_derivs = nan(size(aK_ind));
aK_ind_derivs(2:end,:) = aK_ind(2:end,:) - aK_ind(1:end-1,:);
aK_ind_derivs(1,:) = 0;
aK_ind_pen = reg_params.adl1_ind*sum(aK_ind_derivs(:).^2);

aK_ind_lapl = nan(size(aK_ind));
aK_ind_lapl(2:end-1,:) = 2*aK_ind(2:end-1,:) - aK_ind(1:end-2,:)-aK_ind(3:end,:);
aK_ind_lapl(1,:) = 2*aK_ind(1,:)-2*aK_ind(2,:);
aK_ind_lapl(end,:) = 2*aK_ind(end,:)-2*aK_ind(end-1,:);

LLgrad(ak_ind_inds) = LLgrad(ak_ind_inds) - 2*reg_params.adl1_ind*aK_ind_lapl(:);

aK_dep_derivs = nan(size(aK_dep));
aK_dep_derivs(2:end,:) = aK_dep(2:end,:) - aK_dep(1:end-1,:);
aK_dep_derivs(1,:) = 0;
aK_dep_pen = reg_params.adl1_ind*sum(aK_dep_derivs(:).^2);

aK_dep_lapl = nan(size(aK_dep));
aK_dep_lapl(2:end-1,:) = 2*aK_dep(2:end-1,:) - aK_dep(1:end-2,:)-aK_dep(3:end,:);
aK_dep_lapl(1,:) = 2*aK_dep(1,:) - 2*aK_dep(2,:);
aK_dep_lapl(end,:) = 2*aK_dep(end,:) - 2*aK_dep(end-1,:);
LLgrad(ak_dep_inds) = LLgrad(ak_dep_inds) - 2*reg_params.adl1_dep*aK_dep_lapl(:);

%% COMPUTE AMP D2 PENS
aK_ind_d2_pen = reg_params.adl2_ind*sum(aK_ind_lapl(:).^2);

aK_ind_d2grad = zeros(size(aK_ind));
aK_ind_d2grad(2:end-1,:) = 2*aK_ind_lapl(2:end-1,:) - aK_ind_lapl(1:end-2,:)-aK_ind_lapl(3:end,:);
aK_ind_d2grad(1,:) = 2*aK_ind_lapl(1,:)-2*aK_ind_lapl(2,:);
aK_ind_d2grad(end,:) = 2*aK_ind_lapl(end,:)-2*aK_ind_lapl(end-1,:);
LLgrad(ak_ind_inds) = LLgrad(ak_ind_inds) - 2*reg_params.adl2_ind*aK_ind_d2grad(:);

aK_dep_d2_pen = reg_params.adl2_dep*sum(aK_dep_lapl(:).^2);

aK_dep_d2grad = zeros(size(aK_ind));
aK_dep_d2grad(2:end-1,:) = 2*aK_dep_lapl(2:end-1,:) - aK_dep_lapl(1:end-2,:)-aK_dep_lapl(3:end,:);
aK_dep_d2grad(1,:) = 2*aK_dep_lapl(1,:) - 2*aK_dep_lapl(2,:);
aK_dep_d2grad(end,:) = 2*aK_dep_lapl(end,:) - 2*aK_dep_lapl(end-1,:);
LLgrad(ak_dep_inds) = LLgrad(ak_dep_inds) - 2*reg_params.adl2_dep*aK_dep_lapl(:);


%% COMPUTE AMP FREQ D2 PENS
aK_find_lapl = nan(size(aK_ind));
aK_find_lapl(:,2:end-1) = 2*aK_ind(:,2:end-1)-aK_ind(:,1:end-2)-aK_ind(:,3:end);
aK_find_lapl(:,1) = 2*aK_ind(:,1) - 2*aK_ind(:,2);
aK_find_lapl(:,end) = 2*aK_ind(:,end)-2*aK_ind(:,end-1);
aK_find_d2_pen = reg_params.adl2_freq_ind*sum(aK_find_lapl(:).^2);

aK_find_d2grad = zeros(size(aK_ind));
aK_find_d2grad(:,2:end-1) = 2*aK_find_lapl(:,2:end-1)-aK_find_lapl(:,1:end-2)-aK_find_lapl(:,3:end);
aK_find_d2grad(:,1) = 2*aK_find_lapl(:,1)-2*aK_find_lapl(:,2);
aK_find_d2grad(:,end) = 2*aK_find_lapl(:,end)-2*aK_find_lapl(:,end-1);

LLgrad(ak_ind_inds) = LLgrad(ak_ind_inds) - 2*reg_params.adl2_freq_ind*aK_find_d2grad(:);

aK_fdep_lapl = nan(size(aK_dep));
aK_fdep_lapl(:,2:end-1) = 2*aK_dep(:,2:end-1)-aK_dep(:,1:end-2)-aK_dep(:,3:end);
aK_fdep_lapl(:,1) = 2*aK_dep(:,1) - 2*aK_dep(:,2);
aK_fdep_lapl(:,end) = 2*aK_dep(:,end)-2*aK_dep(:,end-1);
aK_fdep_d2_pen = reg_params.adl2_freq_dep*sum(aK_fdep_lapl(:).^2);

aK_fdep_d2grad = zeros(size(aK_dep));
aK_fdep_d2grad(:,2:end-1) = 2*aK_fdep_lapl(:,2:end-1)-aK_fdep_lapl(:,1:end-2)-aK_fdep_lapl(:,3:end);
aK_fdep_d2grad(:,1) = 2*aK_fdep_lapl(:,1)-2*aK_fdep_lapl(:,2);
aK_fdep_d2grad(:,end) = 2*aK_fdep_lapl(:,end)-2*aK_fdep_lapl(:,end-1);

LLgrad(ak_dep_inds) = LLgrad(ak_dep_inds) - 2*reg_params.adl2_freq_dep*aK_fdep_d2grad(:);


%%
LL = LL - pK_ind_pen - pK_dep_pen - pK_ind_d2_pen - pK_dep_d2_pen - pK_find_d2_pen - pK_fdep_d2_pen;
LL = LL - aK_ind_pen - aK_dep_pen - aK_ind_d2_pen - aK_dep_d2_pen - aK_find_d2_pen - aK_fdep_d2_pen;

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

