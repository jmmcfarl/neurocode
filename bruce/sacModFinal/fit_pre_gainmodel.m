function [preGainMod] = fit_pre_gainmodel(stim_mod,Robs,Xmat,Xsac_mat,poss_d2T,poss_L2,tr_inds,xv_inds)

cur_NT = length(Robs);
stim_dims = stim_mod.stim_params(1).stim_dims;
flen = stim_dims(1);
klen = prod(stim_dims);
n_lags = size(Xsac_mat,2);


%convert Xmat into TxLxD matrix
rXmat = reshape(Xmat,[cur_NT flen stim_dims(2)]);

%precompute the output of each filter using only time points at each
%latency to saccades
sac_filt_outs = nan(n_lags,length(stim_mod.mods),cur_NT);
all_Filts = [stim_mod.mods(1:end).filtK];
for pp = 1:n_lags
    sac_emb = create_time_embedding(Xsac_mat(:,pp),NMMcreate_stim_params(flen));
    sac_Xmat = reshape(bsxfun(@times,rXmat,sac_emb),[cur_NT klen]);
    sac_filt_outs(pp,:,:) = (sac_Xmat*all_Filts)';
end

%defaults for minFunc are 1e-5 and 1e-9
optim_params.optTol = 1e-5;
optim_params.progTol = 1e-8;
optim_params.Method = 'lbfgs';
optim_params.Display = 'off';

%initialize regularization params
rp = NMMcreate_reg_params('lambda_d2T',1,'boundary_conds',[0 0 0]);
nim = NMMinitialize_model(NMMcreate_stim_params(n_lags),1,{'lin'},rp);
L2_mats = create_L2_matrices_NMM( nim );

if length(poss_d2T) > 1 || length(poss_L2) > 1
    L2_xvLL = nan(length(poss_d2T),length(poss_L2));
    for ii = 1:length(poss_d2T)
        cur_lambda_d2T = poss_d2T(ii);
        for jj = 1:length(poss_L2)
            cur_lambda_L2 = poss_L2(jj);
            
            fprintf('Fitting pre-model d2T %d/%d, L2 %d/%d\n',ii,length(poss_d2T),jj,length(poss_L2));
            initial_params = [zeros(n_lags,1); zeros(n_lags,1); stim_mod.spk_NL_params(1)];
            %fit upstream filter and offset filter, given post-gain filter
            [params,mval] = minFunc( @(K) LLinternal_alphas_full(stim_mod, K, Robs(tr_inds), rXmat(tr_inds,:,:), Xsac_mat(tr_inds,:) ...
                , sac_filt_outs(:,:,tr_inds), L2_mats, cur_lambda_d2T,cur_lambda_L2), initial_params, optim_params);
            [~,~,L2_xvLL(ii,jj)] = LLinternal_alphas_full(stim_mod, params, Robs(xv_inds), rXmat(xv_inds,:,:),...
                Xsac_mat(xv_inds,:) , sac_filt_outs(:,:,xv_inds), L2_mats, cur_lambda_d2T,cur_lambda_L2);
            
        end
    end
    
    [~,optloc] = min(L2_xvLL(:));
    [optloc_x,optloc_y] = ind2sub([length(poss_d2T) length(poss_L2)],optloc);
    opt_d2T = poss_d2T(optloc_x);
    opt_L2 = poss_L2(optloc_y);
    
else
    opt_d2T = poss_d2T; opt_L2 = poss_L2; L2_xvLL = nan;
end

initial_params = [zeros(n_lags,1); zeros(n_lags,1); 1; stim_mod.spk_NL_params(1)];
%fit upstream filter and offset filter, given post-gain filter
[params,mval] = minFunc( @(K) LLinternal_alphas_full(stim_mod, K, Robs, rXmat, Xsac_mat, ...
    sac_filt_outs, L2_mats, opt_d2T,opt_L2), initial_params, optim_params);

preGainMod.stim_kernel = params(1:n_lags);
preGainMod.off_kernel = params((n_lags+1):2*n_lags);
preGainMod.theta = params(end);
preGainMod.beta = params(end-1);
preGainMod.stim_mod = stim_mod;
preGainMod.opt_d2T = opt_d2T;
preGainMod.opt_L2 = opt_L2;
preGainMod.fullxvLL = L2_xvLL;

[~,~,LL] = LLinternal_alphas_full(stim_mod, params, Robs, rXmat,...
    Xsac_mat , sac_filt_outs, L2_mats, opt_d2T,opt_L2);

preGainMod.LL = -LL;

end
%%
function [LL,LLgrad,LLraw] = LLinternal_alphas_full( mod, params, Robs, Xmat, Xsac_mat, sac_filt_outs, L2_mats,lambda,lambda_L2)

%F[ kX{I + gamma(tau)J(tau)} g(tau) + c(tau)]

NT = length(Robs);
n_lags = size(Xsac_mat,2);
flen = mod.stim_params(1).stim_dims(1);
klen = prod(mod.stim_params(1).stim_dims);

%create X * (gamma * J) (summed over tau). Same dimensions as Xmat.
%This is the additive component of Xmat with upstream filtering
sac_emb = create_time_embedding(Xsac_mat*params(1:n_lags),NMMcreate_stim_params(flen));
sac_Xmat = bsxfun(@times,Xmat,sac_emb);

%add Xmat components with and without upstream filtering
cur_Xmat = reshape(Xmat + sac_Xmat,[NT klen]);

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
Nmods = length(mod.mods);
all_Filts = [mod.mods(1:end).filtK];

theta = params(end); % offset
beta = params(end-1); %overall gain scaling
Gstim = 0;
gint = nan(length(Robs),Nmods);

for ii = 1:Nmods
    
    gint(:,ii) = cur_Xmat*all_Filts(:,ii);
    
    % Process subunit g's with upstream NLs
    if strcmp(mod.mods(ii).NLtype,'nonpar')
        fgint = piecelin_process(gint(:,ii),mod.mods(ii).NLy,mod.mods(ii).NLx);
    elseif strcmp(mod.mods(ii).NLtype,'quad')
        fgint = gint(:,ii).^2;
    elseif strcmp(mod.mods(ii).NLtype,'lin')
        fgint = gint(:,ii);
    elseif strcmp(mod.mods(ii).NLtype,'threshlin')
        fgint = gint(:,ii);
        fgint(fgint < 0) = 0;
    else
        error('Invalid internal NL');
    end
    
    Gstim = Gstim + fgint*mod.mods(ii).sign;
    
end

%add in sac offset term
G = theta + Gstim*beta + Xsac_mat*params((n_lags+1):(2*n_lags));
%% Compute predicted firing rate
if strcmp(mod.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*mod.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    r = mod.spk_NL_params(4) + mod.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    r(too_large) = mod.spk_NL_params(4) + mod.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
elseif strcmp(mod.spk_NL_type,'exp')
    expg = exp(G);
    r = expg;
elseif strcmp(mod.spk_NL_type,'linear')
    r = G;
else
    error('invalid spk nl');
end

% Enforce mimodum predicted firing rate to avoid nan LLs
if ~strcmp(mod.spk_NL_type,'linear')
    min_pred_rate = 1e-50;
    if min(r) < min_pred_rate
        r(r < min_pred_rate) = min_pred_rate; %mimodum predicted rate
    end
end
%% COMPUTE LL and LL gradient
if strcmp(mod.spk_NL_type,'linear') % use MSE as cost function
    Nspks = length(Robs);
    LL = -sum( (Robs - r).^2 );
else
    Nspks = sum(Robs);
    LL = sum(Robs.* log(r) - r); %up to an overall constant
    %'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
end

%'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
if strcmp(mod.spk_NL_type,'logexp')
    residual = mod.spk_NL_params(3)*mod.spk_NL_params(2)*(Robs./r - 1) .* expg ./ (1+expg);
    residual(too_large) = mod.spk_NL_params(3)*mod.spk_NL_params(2)*(Robs(too_large)./r(too_large) - 1);
elseif strcmp(mod.spk_NL_type,'exp')
    residual = Robs - r;
elseif strcmp(mod.spk_NL_type,'linear')
    residual = 2*(Robs - r);
else
    error('Unsupported spiking NL')
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant term (theta)
LLgrad(end) = sum(residual);
LLgrad(end-1) = sum(residual.*Gstim);
%%

gint = gint*beta; %gain scaling
for pp = 1:n_lags
    
    for ii = 1:Nmods
        
        if strcmp(mod.mods(ii).NLtype,'lin')
            LLgrad(pp) = LLgrad(pp) + (residual.*gint(:,ii)* mod.mods(ii).sign)' *squeeze(sac_filt_outs(pp,ii,:));
        else
            if strcmp(mod.mods(ii).NLtype,'nonpar')
                fpg = piececonst_process(gint(:,ii),fprimes{ii}, mod.mods(ii).NLx);
            elseif strcmp(mod.mods(ii).NLtype,'quad')
                fpg = 2*gint(:,ii);
            elseif strcmp(mod.mods(ii).NLtype,'threshlin')
                fpg = gint(:,ii) >= 0;
                
            else
                error('Unsupported NL type')
            end
            LLgrad(pp) = LLgrad(pp) + (fpg.*residual*mod.mods(ii).sign)' * squeeze(sac_filt_outs(pp,ii,:));
        end
    end
    
end

cur_set = (n_lags+1):2*n_lags;
LLgrad(cur_set) = residual'*Xsac_mat;

%%
smooth_penalty = 0;
LLgrad_pen = zeros(size(LLgrad));
if lambda > 0
    smooth_penalty = smooth_penalty + lambda*sum((L2_mats.L2_d2T{1} * params(1:n_lags)).^2);
    LLgrad_pen(1:n_lags) = 2*lambda*(L2_mats.L2_d2T{1}' * L2_mats.L2_d2T{1} * params(1:n_lags));
    
    cur_set = (n_lags+1):2*n_lags;
    smooth_penalty = smooth_penalty + lambda*sum((L2_mats.L2_d2T{1} * params(cur_set)).^2);
    LLgrad_pen(cur_set) = 2*lambda*(L2_mats.L2_d2T{1}' * L2_mats.L2_d2T{1} * params(cur_set));
end
if lambda_L2 > 0
    L2_penalty = sum(params(1:2*n_lags).^2)*lambda_L2;
    LLgrad_pen(1:2*n_lags) = LLgrad_pen(1:2*n_lags) + 2*lambda_L2*params(1:2*n_lags);
else
    L2_penalty = 0;
end
LLraw = LL;
LL = LL - smooth_penalty - L2_penalty;
LLgrad = LLgrad - LLgrad_pen;

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
LLraw = -LLraw/Nspks;
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;
end
