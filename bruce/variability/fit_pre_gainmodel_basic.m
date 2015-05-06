function [preGainMod] = fit_pre_gainmodel_basic(stim_mod,Robs,Xmat,Xsac_mat,d2T_off,d2T_gain)

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


%create a new linear model
temp_stim_params(1) = NMMcreate_stim_params(1);
temp_stim_params(2) = NMMcreate_stim_params(n_lags);
temp_mod = NMMinitialize_model(temp_stim_params,[1 1],{'lin','lin'},[],[1 2]);
temp_mod.mods(1).filtK = 1;
temp_mod.spk_NL_params = stim_mod.spk_NL_params;
X{2} = Xsac_mat;

null_prate = mean(Robs);
nullLL = sum(Robs.*log(ones(size(Robs))*null_prate) - ones(size(Robs))*null_prate)/sum(Robs);

%initialize regularization params
rp = NMMcreate_reg_params('lambda_d2T',1,'boundary_conds',[Inf 0 0]);
nim = NMMinitialize_model(NMMcreate_stim_params(n_lags),1,{'lin'},rp);
L2_mats = create_L2_matrices_NMM( nim );


initial_params = [zeros(n_lags,1); zeros(n_lags,1); stim_mod.spk_NL_params(1)];
%fit upstream filter and offset filter, given post-gain filter

[params,mval] = minFunc( @(K) LLinternal_alphas_full(stim_mod, K, Robs, rXmat, Xsac_mat...
    ,sac_filt_outs, L2_mats, d2T_off,d2T_gain), initial_params, optim_params);

sac_emb = create_time_embedding(Xsac_mat*params(1:n_lags),NMMcreate_stim_params(flen));
sac_Xmat = bsxfun(@times,rXmat,sac_emb);
%add Xmat components with and without upstream filtering
cur_Xmat = reshape(rXmat + sac_Xmat,[length(Robs) klen]);
%stim output (with upstream sac mod)
[~,~,~,stimG] = NMMmodel_eval(stim_mod,Robs,cur_Xmat);
stimG = stimG - stim_mod.spk_NL_params(1);
X{1} = stimG;
temp_mod.spk_NL_params(1) = params(end);
temp_mod.mods(2).filtK = params((n_lags+1):2*n_lags);
temp_mod = NMMfit_logexp_spkNL(temp_mod,Robs,X);

preGainMod.stim_kernel = params(1:n_lags);
preGainMod.off_kernel = params((n_lags+1):2*n_lags);
preGainMod.theta = temp_mod.spk_NL_params(1);
preGainMod.stim_mod = stim_mod;
preGainMod.stim_mod.spk_NL_params = temp_mod.spk_NL_params;

[LL,pred_rate] = eval_pre_gainmodel( preGainMod, Robs, Xmat, Xsac_mat);
preGainMod.ovInfo = mean(pred_rate/mean(pred_rate).*log2(pred_rate/mean(pred_rate)));
preGainMod.ovLLimp = (LL-nullLL)/log(2);
preGainMod.nullLL = nullLL;

end

%%
function [LL,LLgrad,LLraw] = LLinternal_alphas_full( mod, params, Robs, Xmat, Xsac_mat, sac_filt_outs, L2_mats,lambda_off,lambda_gain)

%F[ kX{I + gamma(tau)J(tau)} + c(tau)]

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
% beta = params(end-1); %overall gain scaling
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
% G = theta + Gstim*beta + Xsac_mat*params((n_lags+1):(2*n_lags));
G = theta + Gstim + Xsac_mat*params((n_lags+1):(2*n_lags));
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
% LLgrad(end-1) = sum(residual.*Gstim);
%%

for pp = 1:n_lags
    
    for ii = 1:Nmods
        
        if strcmp(mod.mods(ii).NLtype,'lin')
            %             LLgrad(pp) = LLgrad(pp) + beta*(residual * mod.mods(ii).sign)' *squeeze(sac_filt_outs(pp,ii,:));
            LLgrad(pp) = LLgrad(pp) + (residual * mod.mods(ii).sign)' *squeeze(sac_filt_outs(pp,ii,:));
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
            %             LLgrad(pp) = LLgrad(pp) + beta*(fpg.*residual*mod.mods(ii).sign)' * squeeze(sac_filt_outs(pp,ii,:));
            LLgrad(pp) = LLgrad(pp) + (fpg.*residual*mod.mods(ii).sign)' * squeeze(sac_filt_outs(pp,ii,:));
        end
    end
    
end

cur_set = (n_lags+1):2*n_lags;
LLgrad(cur_set) = residual'*Xsac_mat;

%%
smooth_penalty = 0;
LLgrad_pen = zeros(size(LLgrad));
if lambda_gain > 0 || lambda_off > 0
    smooth_penalty = smooth_penalty + lambda_gain*sum((L2_mats.L2_d2T{1} * params(1:n_lags)).^2);
    LLgrad_pen(1:n_lags) = 2*lambda_gain*(L2_mats.L2_d2T{1}' * L2_mats.L2_d2T{1} * params(1:n_lags));
    
    cur_set = (n_lags+1):2*n_lags;
    smooth_penalty = smooth_penalty + lambda_off*sum((L2_mats.L2_d2T{1} * params(cur_set)).^2);
    LLgrad_pen(cur_set) = 2*lambda_off*(L2_mats.L2_d2T{1}' * L2_mats.L2_d2T{1} * params(cur_set));
end
L2_penalty = 0;
LLraw = LL;
LL = LL - smooth_penalty - L2_penalty;
LLgrad = LLgrad - LLgrad_pen;

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
LLraw = -LLraw/Nspks;
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;
end

function [LL,LLgrad] = LLinternal_alphas_gonly( mod, params, Robs, Xmat, Xsac_mat,off_kernel,stim_kernel)

%F[ kX{I + gamma(tau)J(tau)} + c(tau)]

NT = length(Robs);
n_lags = size(Xsac_mat,2);
flen = mod.stim_params(1).stim_dims(1);
klen = prod(mod.stim_params(1).stim_dims);

%create X * (gamma * J) (summed over tau). Same dimensions as Xmat.
%This is the additive component of Xmat with upstream filtering
sac_emb = create_time_embedding(Xsac_mat*stim_kernel,NMMcreate_stim_params(flen));
sac_Xmat = bsxfun(@times,Xmat,sac_emb);

%add Xmat components with and without upstream filtering
cur_Xmat = reshape(Xmat + sac_Xmat,[NT klen]);

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
Nmods = length(mod.mods);
all_Filts = [mod.mods(1:end).filtK];

beta = params(1); %overall gain scaling
theta = params(2); % offset
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
G = theta + beta*Gstim + Xsac_mat*off_kernel;
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
LLgrad(1) = sum(residual.*Gstim);
LLgrad(2) = sum(residual);

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;
end

