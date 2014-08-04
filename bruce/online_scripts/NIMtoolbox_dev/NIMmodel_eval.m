function [LL, penLL, pred_rate, G, gint, fgint] = NIMmodel_eval(nim,Robs,Xstim,XLin)
%
% [LL, penLL, pred_rate, G, gint, fgint] = NIMmodel_eval(nim,Robs,Xstim,<XLin>)
%
% Evalutes the LL of the specified model, also returns other useful outputs
% of the model
%
% INPUTS:
%   nim: model structure
%   Robs: binned spikes
%   Xstim: time-embedded stimulus mat
%   <XLin>: Matrix specifying additional linear predictors
%
% OUTPUTS:
%   LL: log-likelihood (per spike) of the model
%   penLL: penalized log-likelihood (per spike)
%   pred_rate: predicted firing rate (unitless, so divide by dt to get Hz)
%   G: generating function (output of the model before the spk NL)
%   gint: TxNmods matrix of the output of each subunit's stimulus filter (before applying upstream NL)
%   fgint: TxNmods matrix of the output of each subunit (after applying upstream NL)


%%
if nargin < 4
    XLin = [];
end

%% Key parameters
NT = length(Robs);

if NT == 0               %%%% DAB ADDED
  NT = size(Xstim,1);	   %%%% DAB ADDED
end                      %%%% DAB ADDED
	
Nmods = length(nim.mods);
lin_dims = nim.stim_params.lin_dims;
spkhstlen = nim.spk_hist.spkhstlen;

%% CREATE L2 REGULARIZATION MATRICES
L2_mats = create_L2_matrices(nim);

%% CREATE SPIKE HISTORY MATRIX IF NEEDED
if spkhstlen > 0
    Xspkhst = create_spkhist_Xmat(Robs,nim.spk_hist.bin_edges);
else
    Xspkhst = [];
end

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
theta = nim.spk_NL_params(1); %offset
G = theta + zeros(NT,1); %initialize overall generating function G

Kmat = [nim.mods(:).filtK];
gint = Xstim*Kmat; %subunit generating functions

fgint = nan(NT,Nmods);
for n = 1:Nmods
    
    %process subunit g's with upstream NLs
    if strcmp(nim.mods(n).NLtype,'nonpar')
        fgint(:,n) = piecelin_process(gint(:,n),nim.mods(n).NLy,nim.mods(n).NLx);
    elseif strcmp(nim.mods(n).NLtype,'quad')
        fgint(:,n) = gint(:,n).^2;
    elseif strcmp(nim.mods(n).NLtype,'lin')
        fgint(:,n) = gint(:,n);
    elseif strcmp(nim.mods(n).NLtype,'threshlin')
        fgint(:,n) = gint(:,n);
        fgint(fgint(:,n) < 0,n) = 0;
    else
        error('Invalid internal NL');
    end
    
    %multiply by weight and add to generating function
    G = G + fgint(:,n)*nim.mods(n).sign;
end

%add spike history contribution
if spkhstlen > 0
    G = G + Xspkhst*nim.spk_hist.coefs;
end
%add contribution from linear filter
if lin_dims > 0
    G = G + XLin*nim.kLin;
end

%% IF YOU JUST WANT TO COMPUTE G gint and fgint
if isempty(Robs)
    LL = nan;
    penLL = nan;
    pred_rate = nan;
    return
end

%% OTHERWISE COMPUTE LL AND PREDICTED RATE
%make sure Robs is a column vector
if size(Robs,2) > size(Robs,1)
    Robs = Robs';
end

%% Compute predicted firing rate and LL
if strcmp(nim.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*nim.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    pred_rate = nim.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    pred_rate(too_large) = nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
elseif strcmp(nim.spk_NL_type,'exp')
    expg = exp(G);
    pred_rate = expg;
else
    error('invalid spk nl');
end
%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(pred_rate) < min_pred_rate
    pred_rate(pred_rate < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

LL = sum(Robs.* log(pred_rate) - pred_rate); %up to an overall constant

%% COMPUTE L2 PENALTIES
smooth_penalty = zeros(Nmods,1);
deriv_penalty = zeros(Nmods,1);
ridge_penalty = zeros(Nmods,1);
sparse_penalty = zeros(Nmods,1);

for n = 1:Nmods
    cur_kern = nim.mods(n).filtK;
    if nim.mods(n).reg_params.lambda_dT > 0
        deriv_penalty(n) = deriv_penalty(n) + nim.mods(n).reg_params.lambda_dT*sum((L2_mats.L2_dT * cur_kern).^2);
    end
    if nim.mods(n).reg_params.lambda_d2T > 0
        smooth_penalty(n) = smooth_penalty(n) + nim.mods(n).reg_params.lambda_d2T*sum((L2_mats.L2_d2T * cur_kern).^2);
    end
%     if nim.stim_params.stim_dims(2) > 1 && nim.stim_params.stim_dims(3) == 1
    if nim.stim_params.stim_dims(2) > 1 %if spatial dims
        if nim.mods(n).reg_params.lambda_dX > 0
            deriv_penalty(n) = deriv_penalty(n) + nim.mods(n).reg_params.lambda_dX*sum((L2_mats.L2_dX * cur_kern).^2);
        end
        if nim.mods(n).reg_params.lambda_d2X > 0
            smooth_penalty(n) = smooth_penalty(n) + nim.mods(n).reg_params.lambda_d2X*sum((L2_mats.L2_d2X * cur_kern).^2);
        end
        if nim.mods(n).reg_params.lambda_d2XT > 0
            smooth_penalty(n) = smooth_penalty(n) + nim.mods(n).reg_params.lambda_d2XT*sum((L2_mats.L2_d2XT * cur_kern).^2);
        end
%     elseif nim.stim_params.stim_dims(3) > 1
%         error('No support for 2-spatial dims yet');
    end
    if nim.mods(n).reg_params.lambda_L2 > 0
        ridge_penalty(n) = nim.mods(n).reg_params.lambda_L2*(cur_kern' * cur_kern);
    end
    if nim.mods(n).reg_params.lambda_L1 > 0
       sparse_penalty(n) = nim.mods(n).reg_params.lambda_L1*sum(abs(cur_kern)); 
    end
end

%penalized LL
penLL = LL - sum(smooth_penalty) - sum(ridge_penalty) - sum(deriv_penalty) - sum(sparse_penalty);

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
Nspks = sum(Robs);
LL = LL/Nspks;
penLL = penLL/Nspks;
