function [LL, penLL, pred_rate, G, nullLL] = NMMmodel_eval_EM( nim, Robs, Xstims, gamma,shifts, Gmults, regmat_custom)
%
% Usage: [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMmodel_eval( nim, Robs, Xstims, <Gmults>,<regmat_custom> )
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
%   nullLL: log-likelihood (per spike) of the null model (constant firing rate)


%% Set parameters
Nmods = length(nim.mods);
spkhstlen = nim.spk_hist.spkhstlen;

if (nargin < 6) || (length(Gmults) < Nmods)
  Gmults{Nmods} = [];
end
if nargin < 7
    regmat_custom = [];
end

% Process Xstims (in case multiple Xstims)
if ~iscell(Xstims)
	tmp = Xstims;
	clear Xstims
	Xstims{1} = tmp;
end

%% Key parameters
NT = length(Robs);

if NT == 0    
  NT = size(Xstims{1},1);
end

%% CREATE L2 REGULARIZATION MATRICES
L2_mats = create_L2_matrices_NMM(nim);
L2_mats.custom = regmat_custom;

%% CREATE SPIKE HISTORY MATRIX IF NEEDED
if spkhstlen > 0
  Xspkhst = create_spkhist_Xmat( Robs, nim.spk_hist.bin_edges );
else
  Xspkhst = [];
end

%% COMPUTE SHIFT MATRICES
n_shifts = length(shifts);
sdims = nim.stim_params(1).stim_dims;
L = prod(sdims);
for n = 1:n_shifts
    shift_mat{n} = spdiags( ones(L,1), -shifts(n)*sdims(1), L, L );
end

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
theta = nim.spk_NL_params(1); % offset
G = theta + zeros(NT,1); % initialize overall generating function G
G = repmat(G,1,n_shifts);

%Kmat = [nim.mods(:).filtK];
%gint = Xstim*Kmat; % subunit generating functions
Xtargets = [nim.mods(:).Xtarget];
stim_targets = find(Xtargets == 1);

% gint = nan(NT,Nmods);
% fgint = nan(NT,Nmods);
for n = 1:Nmods
    
    if ismember(n,stim_targets)
        ks = zeros(length(nim.mods(n).filtK),n_shifts);
        for ii = 1:n_shifts
            ks(:,ii) = shift_mat{ii} * nim.mods(n).filtK;
        end
        
        gint = Xstims{nim.mods(n).Xtarget} * ks;
    else
        gint = repmat(Xstims{nim.mods(n).Xtarget} * nim.mods(n).filtK,1,n_shifts);
    end
    
  % Process subunit g's with upstream NLs
	if strcmp(nim.mods(n).NLtype,'nonpar')
		fgint = piecelin_process( gint, nim.mods(n).NLy, nim.mods(n).NLx );
	elseif strcmp(nim.mods(n).NLtype,'quad')
		fgint = gint.^2;
	elseif strcmp(nim.mods(n).NLtype,'lin')
		fgint = gint;
	elseif strcmp(nim.mods(n).NLtype,'threshlin')
		fgint = gint;
		fgint(fgint < 0,n) = 0;
	else
		error('Invalid internal NL');
	end
    
	% Multiply by weight (and multiplier, if appl) and add to generating function
	if isempty(Gmults{n})
		G = G + fgint * nim.mods(n).sign;
	else
		G = G + (fgint.*Gmults{n}) * nim.mods(n).sign;
	end
	
end

%add spike history contribution
if spkhstlen > 0
	G = G + Xspkhst*nim.spk_hist.coefs;
end

%% OTHERWISE COMPUTE LL AND PREDICTED RATE
% Make sure Robs is a column vector
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
elseif strcmp(nim.spk_NL_type,'linear')
    pred_rate = G;    
else
    error('invalid spk nl');
end
%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(pred_rate) < min_pred_rate
    pred_rate(pred_rate < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

%% IF YOU JUST WANT TO COMPUTE G, gint, fgint and pred_rate
if isempty(Robs)
    LL = nan;
    penLL = nan;
    nullLL = nan;
    return
end

if strcmp(nim.spk_NL_type,'linear') % use MSE as cost function 
    Nspks = length(Robs);
    LL = - bsxfun(@minus,pred_rate,Robs).^2 ;
else
    Nspks = sum(Robs);
    LL = bsxfun(@times,log(pred_rate),Robs) - pred_rate;
end
LL = sum(sum(LL.*gamma));

%% COMPUTE L2 PENALTIES
smooth_penalty = zeros(Nmods,1);
deriv_penalty = zeros(Nmods,1);
ridge_penalty = zeros(Nmods,1);
sparse_penalty = zeros(Nmods,1);
custom_penalty = zeros(Nmods,1);
for n = 1:Nmods
    cur_kern = nim.mods(n).filtK;
    if nim.mods(n).reg_params.lambda_dT > 0
        deriv_penalty(n) = deriv_penalty(n) + nim.mods(n).reg_params.lambda_dT*sum((L2_mats.L2_dT{nim.mods(n).Xtarget} * cur_kern).^2);
    end
    if nim.mods(n).reg_params.lambda_d2T > 0
        smooth_penalty(n) = smooth_penalty(n) + nim.mods(n).reg_params.lambda_d2T*sum((L2_mats.L2_d2T{nim.mods(n).Xtarget} * cur_kern).^2);
    end
%     if nim.stim_params.stim_dims(2) > 1 && nim.stim_params.stim_dims(3) == 1
    if nim.stim_params(nim.mods(n).Xtarget).stim_dims(2) > 1 %if spatial dims
        if nim.mods(n).reg_params.lambda_dX > 0
            deriv_penalty(n) = deriv_penalty(n) + nim.mods(n).reg_params.lambda_dX*sum((L2_mats.L2_dX{nim.mods(n).Xtarget} * cur_kern).^2);
        end
        if nim.mods(n).reg_params.lambda_d2X > 0
            smooth_penalty(n) = smooth_penalty(n) + nim.mods(n).reg_params.lambda_d2X*sum((L2_mats.L2_d2X{nim.mods(n).Xtarget} * cur_kern).^2);
        end
        if nim.mods(n).reg_params.lambda_d2XT > 0
            smooth_penalty(n) = smooth_penalty(n) + nim.mods(n).reg_params.lambda_d2XT*sum((L2_mats.L2_d2XT{nim.mods(n).Xtarget} * cur_kern).^2);
        end
    end
    if nim.mods(n).reg_params.lambda_L2 > 0
        ridge_penalty(n) = nim.mods(n).reg_params.lambda_L2*(cur_kern' * cur_kern);
    end
    if nim.mods(n).reg_params.lambda_L1 > 0
        sparse_penalty(n) = nim.mods(n).reg_params.lambda_L1*sum(abs(cur_kern));
    end
    %for custom regularization
    if nim.mods(n).reg_params.lambda_custom > 0
			if isempty(L2_mats.custom)
				disp('Warning: penalty calculation is off because L2_mats.custom is not included.')
			else
        custom_penalty(n) = custom_penalty(n) + nim.mods(n).reg_params.lambda_custom*sum((L2_mats.custom * cur_kern).^2);
			end
    end
    
end

%penalized LL
penLL = LL - sum(smooth_penalty) - sum(ridge_penalty) - sum(deriv_penalty) - sum(sparse_penalty) - sum(custom_penalty);

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS

LL = LL/Nspks;
penLL = penLL/Nspks;

%% compute null LL
if nargout > 6
    avg_rate = mean(Robs);
    null_prate = ones(NT,1)*avg_rate;
    if strcmp(nim.spk_NL_type,'linear')
        nullLL = sum( (Robs - avg_rate).^2 );
    else
        nullLL = sum(Robs.* log(null_prate) - null_prate); %up to an overall constant
    end
    nullLL = nullLL/Nspks;
end