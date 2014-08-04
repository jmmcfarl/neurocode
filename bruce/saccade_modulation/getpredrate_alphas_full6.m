function [LL,pred_rate] = getpredrate_alphas_full6( mod, params, Robs, Xmat, Xsac_mat, gain_kernel)

NT = length(Robs);
n_lags = size(Xsac_mat,2);
flen = mod.stim_params(1).stim_dims(1);
klen = prod(mod.stim_params(1).stim_dims);
sac_emb = create_time_embedding(Xsac_mat*params(1:n_lags),NMMcreate_stim_params(flen));
sac_Xmat = bsxfun(@times,Xmat,sac_emb);

cur_Xmat = reshape(Xmat + sac_Xmat,[NT klen]);

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
Nmods = length(mod.mods);
all_Filts = [mod.mods(:).filtK];

theta = params(end); % offset
% theta = mod.spk_NL_params(1);
% G = theta; % initialize overall generating function G
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
Gstim = Gstim + Gstim.*(Xsac_mat*gain_kernel);
G = theta + Gstim + Xsac_mat*params((n_lags+1):(2*n_lags));
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
LL = LL/Nspks;
pred_rate = r;