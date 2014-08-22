function [sacGainMod,sacGainOnlyMod] = fit_prepost_gainmodel(stim_mod,Robs,Xmat,Xsac_mat,lambda_d2T,lambda_L2,init_gain_kernel,max_iter)

cur_NT = length(Robs);
stim_dims = stim_mod.stim_params(1).stim_dims;
flen = stim_dims(1);
klen = prod(stim_dims);
n_lags = size(Xsac_mat,2);

%if no initial post-gain kernel is specified, assume it's all zeros
if nargin < 7 || isempty(init_gain_kernel)
    init_gain_kernel = zeros(n_lags,1);
end
if nargin < 8
    max_iter = Inf;
end

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
silent = 1;
optim_params.optTol = 1e-5;
optim_params.progTol = 1e-8;
optim_params.Method = 'lbfgs';
optim_params.Display = 'off';

%initialize regularization params
rp = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'boundary_conds',[0 0 0]);
nim = NMMinitialize_model(NMMcreate_stim_params(n_lags),1,{'lin'},rp);
L2_mats = create_L2_matrices_NMM( nim );

initial_params = [zeros(n_lags,1); zeros(n_lags,1); stim_mod.spk_NL_params(1)];
cur_params = initial_params;
gain_kernel = init_gain_kernel;

cur_LLimp = -Inf;
cur_LL = 1e5;
imp_thresh = -1e-6;
it = 1;
while cur_LLimp <= imp_thresh && it <= max_iter
    
    %fit upstream filter and offset filter, given post-gain filter
    fprintf('Iteration %d\n',it);
    [params{it},mval(it)] = minFunc( @(K) LLinternal_alphas_full(stim_mod, K, Robs, rXmat, Xsac_mat ...
        ,gain_kernel, sac_filt_outs, L2_mats, lambda_d2T,lambda_L2), cur_params, optim_params);
    cur_params = params{it};
    cur_LLimp = mval(it) - cur_LL;
    cur_LL = mval(it);
    
    %after one iteration save the model (this has pre-gain but no
    %post-gain)
    if it == 1
        sacGainOnlyMod.gain_kernel = zeros(n_lags,1);
        sacGainOnlyMod.stim_kernel = cur_params(1:n_lags);
        sacGainOnlyMod.off_kernel = cur_params((n_lags+1):2*n_lags);
        sacGainOnlyMod.theta = cur_params(end);
        sacGainOnlyMod.stim_mod = stim_mod;
        sacGainOnlyMod.LL = eval_sacgain_mod(sacGainOnlyMod,Robs,Xmat,Xsac_mat);
    end
    
    %create upstream filtered Xmat
    sac_emb = create_time_embedding(Xsac_mat*cur_params(1:n_lags),NMMcreate_stim_params(flen));
    sac_Xmat = bsxfun(@times,rXmat,sac_emb);
    cur_Xmat = reshape(rXmat + sac_Xmat,[cur_NT klen]);
    
    %get output of model with this upstream filtered Xmat
    [LL, penLL, pred_rate, G] = NMMmodel_eval(stim_mod,Robs,cur_Xmat);
    g_tot = G - stim_mod.spk_NL_params(1);
    Xsac_tot = bsxfun(@times,Xsac_mat,g_tot);
    
    %now fit the post-gain filter and offset filter, given this upstream
    %filter
    tr_stim{1} = [g_tot];
    tr_stim{2} = Xsac_mat;
    tr_stim{3} = Xsac_tot;
    clear sac_stim_params
    sac_stim_params(1) = NMMcreate_stim_params(1);
    sac_stim_params(2) = NMMcreate_stim_params(size(Xsac_tot,2));
    sac_stim_params(3) = NMMcreate_stim_params(size(Xsac_tot,2));
    sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
    
    mod_signs = [1 1 1];
    Xtargets = [1 2 3];
    NL_types = {'lin','lin','lin'};
    temp_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    temp_mod.mods(1).filtK = 1;
    temp_mod.mods(2).filtK = cur_params((n_lags+1):2*n_lags);
    temp_mod.spk_NL_params(1) = cur_params(end);
    temp_mod.spk_NL_params(2:end) = stim_mod.spk_NL_params(2:end);
    temp_mod.mods(1).reg_params = NMMcreate_reg_params();
    temp_mod = NMMfit_filters(temp_mod,Robs,tr_stim,[],[1 2 3],silent); %fit only offset and post-gain filters (dont refit gain on G)
    gainLL(it) = temp_mod.LL_seq(end);
    
    gain_kernel = temp_mod.mods(3).filtK + temp_mod.mods(1).filtK-1;
    gain_kern_est(it,:) = gain_kernel;
    cur_params((n_lags+1):2*n_lags) = temp_mod.mods(2).filtK; %update estimate of offset filter
    cur_params(end) = temp_mod.spk_NL_params(1);
    it = it + 1;
end

sacGainMod.gain_kernel = gain_kernel;
sacGainMod.stim_kernel = cur_params(1:n_lags);
sacGainMod.off_kernel = cur_params((n_lags+1):2*n_lags);
sacGainMod.theta = cur_params(end);
sacGainMod.stim_mod = stim_mod;
sacGainMod.gain_offset = temp_mod.mods(1).filtK;

sacGainMod.LL = eval_prepost_gainmodel( sacGainMod, Robs, Xmat, Xsac_mat);

end

%%
function [LL,LLgrad] = LLinternal_alphas_full( mod, params, Robs, Xmat, Xsac_mat, gain_kernel, sac_filt_outs, L2_mats,lambda,lambda_L2)

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

%compute multiplicative (post-filtering) gain output
gainfun = Xsac_mat*gain_kernel;
Gstim = Gstim + Gstim.*gainfun;

%add in sac offset term
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

%%

post_out = 1+gainfun;
for pp = 1:n_lags
    
    for ii = 1:Nmods
        
        if strcmp(mod.mods(ii).NLtype,'lin')
            LLgrad(pp) = LLgrad(pp) + (residual.*post_out.*gint(:,ii)* mod.mods(ii).sign)' *squeeze(sac_filt_outs(pp,ii,:));
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
            LLgrad(pp) = LLgrad(pp) + (post_out.*fpg.*residual*mod.mods(ii).sign)' * squeeze(sac_filt_outs(pp,ii,:));
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
LL = LL - smooth_penalty - L2_penalty;
LLgrad = LLgrad - LLgrad_pen;

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS

LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;
end
