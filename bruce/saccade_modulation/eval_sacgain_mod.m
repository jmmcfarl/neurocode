function [LL,pred_rate,Gstim,fgint] = eval_sacgain_mod( sacGainMod, Robs, Xmat, Xsac_mat)

NT = length(Robs);
stim_mod = sacGainMod.stim_mod;

if length(stim_mod.spk_NL_params) < 4
    stim_mod.spk_NL_params = [stim_mod.spk_NL_params 0];
end

n_lags = size(Xsac_mat,2);
flen = stim_mod.stim_params(1).stim_dims(1);
nPix = stim_mod.stim_params(1).stim_dims(2);
klen = prod(stim_mod.stim_params(1).stim_dims);

Xmat = reshape(Xmat,[NT flen nPix]);

sac_emb = create_time_embedding(Xsac_mat*sacGainMod.stim_kernel,NMMcreate_stim_params(flen));
sac_Xmat = bsxfun(@times,Xmat,sac_emb);

cur_Xmat = reshape(Xmat + sac_Xmat,[NT klen]);

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
Nmods = length(stim_mod.mods);
all_Filts = [stim_mod.mods(:).filtK];

theta = sacGainMod.theta; % offset
Gstim = 0;

gint = nan(length(Robs),Nmods);
fgint = nan(length(Robs),Nmods);
for ii = 1:Nmods
            
    gint(:,ii) = cur_Xmat*all_Filts(:,ii);
    
    % Process subunit g's with upstream NLs
    if strcmp(stim_mod.mods(ii).NLtype,'nonpar')
        fgint(:,ii) = piecelin_process(gint(:,ii),stim_mod.mods(ii).NLy,stim_mod.mods(ii).NLx);
    elseif strcmp(stim_mod.mods(ii).NLtype,'quad')
        fgint(:,ii) = gint(:,ii).^2;
    elseif strcmp(stim_mod.mods(ii).NLtype,'lin')
        fgint(:,ii) = gint(:,ii);
    elseif strcmp(stim_mod.mods(ii).NLtype,'threshlin')
        fgint(:,ii) = gint(:,ii);
        fgint(fgint(:,ii) < 0,ii) = 0;
    else
        error('Invalid internal NL');
    end
    
    Gstim = Gstim + fgint(:,ii)*stim_mod.mods(ii).sign;
    
end
Gstim = Gstim + Gstim.*(Xsac_mat*sacGainMod.gain_kernel);
G = theta + Gstim + Xsac_mat*sacGainMod.off_kernel;
%% Compute predicted firing rate
if strcmp(stim_mod.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*stim_mod.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    r = stim_mod.spk_NL_params(4) + stim_mod.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    r(too_large) = stim_mod.spk_NL_params(4) + stim_mod.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
elseif strcmp(stim_mod.spk_NL_type,'exp')
    expg = exp(G);
    r = expg;
elseif strcmp(stim_mod.spk_NL_type,'linear')
    r = G;    
else
    error('invalid spk nl');
end

% Enforce mimodum predicted firing rate to avoid nan LLs
if ~strcmp(stim_mod.spk_NL_type,'linear')
min_pred_rate = 1e-50;
if min(r) < min_pred_rate
    r(r < min_pred_rate) = min_pred_rate; %mimodum predicted rate
end
end
%% COMPUTE LL and LL gradient
if strcmp(stim_mod.spk_NL_type,'linear') % use MSE as cost function 
    Nspks = length(Robs);
    LL = -sum( (Robs - r).^2 );
else
    Nspks = sum(Robs);
    LL = sum(Robs.* log(r) - r); %up to an overall constant
    %'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
end
LL = LL/Nspks;
pred_rate = r;