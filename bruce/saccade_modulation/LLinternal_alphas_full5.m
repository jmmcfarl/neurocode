function [LL,LLgrad] = LLinternal_alphas_full5( mod, params, Robs, Xmat, Xsac_mat, L2_mats,lambda,lambda_L2)

%%
NT = length(Robs);
n_lags = size(Xsac_mat,2);
flen = mod.stim_params(1).stim_dims(1);
klen = prod(mod.stim_params(1).stim_dims);
sac_emb = create_time_embedding(Xsac_mat*params(1:n_lags),NMMcreate_stim_params(flen));
sac_Xmat = bsxfun(@times,Xmat,sac_emb);

cur_Xmat = reshape(Xmat + sac_Xmat,[NT klen]);
%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
Nmods = length(mod.mods);
all_Filts = [mod.mods(1:end).filtK];

theta = params(end); % offset
% theta = mod.spk_NL_params(1);
Gstim = 0; % initialize overall generating function G

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
Xsac_stim = bsxfun(@times,Xsac_mat,Gstim);
G = Gstim + theta;
G = G + Xsac_mat*params((n_lags+1):(2*n_lags)) + Xsac_stim*params((2*n_lags+1):(3*n_lags));
%% Compute predicted firing rate
if strcmp(mod.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*mod.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    r = mod.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    r(too_large) = mod.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
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
    
post_out = 1+Xsac_stim*params((2*n_lags+1):3*n_lags);
for pp = 1:n_lags
    sac_emb = create_time_embedding(Xsac_mat(:,pp),NMMcreate_stim_params(flen));
   sac_Xmat = reshape(bsxfun(@times,Xmat,sac_emb),[NT klen]);

    for ii = 1:Nmods
        
        if strcmp(mod.mods(ii).NLtype,'lin')
            LLgrad(pp) = LLgrad(pp) + (residual.*post_out.*gint(:,ii)* mod.mods(ii).sign)' *(sac_Xmat*all_Filts(:,ii));            
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
            LLgrad(pp) = LLgrad(pp) + (post_out.*fpg.*residual*mod.mods(ii).sign)' * (sac_Xmat*all_Filts(:,ii));
        end
    end

end

cur_set = (n_lags+1):2*n_lags;
LLgrad(cur_set) = residual'*Xsac_mat;
cur_set = (2*n_lags+1):3*n_lags;
LLgrad(cur_set) = residual'*Xsac_stim;

%%
smooth_penalty = 0;
LLgrad_pen = zeros(size(LLgrad));
if lambda > 0
    smooth_penalty = smooth_penalty + lambda*sum((L2_mats.L2_d2T{1} * params(1:n_lags)).^2);
    LLgrad_pen(1:n_lags) = 2*lambda*(L2_mats.L2_d2T{1}' * L2_mats.L2_d2T{1} * params(1:n_lags));
    
    cur_set = (n_lags+1):2*n_lags;
    smooth_penalty = smooth_penalty + lambda*sum((L2_mats.L2_d2T{1} * params(cur_set)).^2);
    LLgrad_pen(cur_set) = 2*lambda*(L2_mats.L2_d2T{1}' * L2_mats.L2_d2T{1} * params(cur_set));
    
    cur_set = (2*n_lags+1):3*n_lags;
    smooth_penalty = smooth_penalty + lambda*sum((L2_mats.L2_d2T{1} * params(cur_set)).^2);
    LLgrad_pen(cur_set) = 2*lambda*(L2_mats.L2_d2T{1}' * L2_mats.L2_d2T{1} * params(cur_set));
end
if lambda_L2 > 0
    L2_penalty = sum(params(1:3*n_lags).^2)*lambda_L2;
    LLgrad_pen(1:3*n_lags) = LLgrad_pen(1:3*n_lags) + 2*lambda_L2*params(1:3*n_lags);
else
    L2_penalty = 0;
end
LL = LL - smooth_penalty - L2_penalty;
LL = LL - smooth_penalty;
LLgrad = LLgrad - LLgrad_pen;

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS

LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;
end
