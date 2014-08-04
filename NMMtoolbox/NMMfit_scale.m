function nim_out = NMMfit_scale( nim, Robs, Xstims, Gmults, Uindx, silent, desired_optim_params )
%
% Usage: nim_out = NMMfit_scale( nim, Robs, Xstims, <Gmults{}>, <Uindx>, <silent>, <desired_optim_params> )
%
% Optimizes the upstream NLs (in terms of tent-basis functions) (plus extra linear terms if desired) for
% given stimulus filters
%
% INPUTS:
%       nim: model structure
%       Robs: binned spikes
%       Xstim: time-embedded stimulus mat
%       <XLin>: Matrix specifying additional linear predictors
%       <desired_targets>: Vector of indices specifying which subunits to optimize. 
%           (-1 = spk history filter, and -2 = extra linear filter) Default is to optimize all elements
%       <rescale_NLs>: set to 0 if you don't want to rescale upstream NLs after
%           estimating (otherwise set to 1)
%       <silent>: set to 0 if you want to turn on the optimization display
%       <desired_optim_params>: Struct of optimization parameters
% OUTPUTS:
%       nim_out: output model struct

%% USEFUL PARAMS
Nmods = length(nim.mods);
spkhstlen = nim.spk_hist.spkhstlen > 0;

%% Process Xstims (in case multiple Xstims)
if ~iscell(Xstims)
	tmp = Xstims;
	clear Xstims
	Xstims{1} = tmp;
end
if nargin < 6
  silent = 1;
end
if nargin < 7
  desired_optim_params = [];
end

NT = size(Xstims{1},1); %stimulus dimensions

%% PARSE INPUTS
if (nargin < 4) || (length(Gmults) < Nmods)
  Gmults{Nmods} = [];
end
if nargin < 5
  Uindx = [];
end

% Index X-matrices and Robs
RobsFULL = Robs;
if ~isempty(Uindx)
  for nn = 1:length(Xstims)
    Xstims{nn} = Xstims{nn}(Uindx,:);
  end
  Robs = RobsFULL(Uindx);
end

% Make sure Robs is a column vector
Robs = Robs(:);

%add spk NL constant if it isnt already there
if length(nim.spk_NL_params) < 4
	nim.spk_NL_params(4) = 0;
end

%% Initialize X-matrix and initial point
X = zeros(NT,Nmods+spkhstlen);  % each mod and spike history term 
initial_params = ones(Nmods+spkhstlen+1,1);
initial_params(end) = nim.spk_NL_params(1);

%% Put in spike-history effects (if applicable)
if spkhstlen > 0
    Xspkhst = create_spkhist_Xmat(Robs,nim.spk_hist.bin_edges);
		X(:,end) = Xspkhst * nim.spk_hist.coefs;
end

%% COMPUTE NET OUPTUT OF ALL SUBUNITS
%Kmat = [nim.mods(:).filtK];
%gint = Xstim*Kmat;
for imod = 1:Nmods
	
		fgint = Xstims{nim.mods(imod).Xtarget} * nim.mods(imod).filtK;
		
    % Process subunit g's with upstream NLs
    if strcmp(nim.mods(imod).NLtype,'nonpar')
        %fgint = piecelin_process(gint(:,imod),nim.mods(imod).NLy,nim.mods(imod).NLx);
				fgint = piecelin_process( fgint, nim.mods(imod).NLy ,nim.mods(imod).NLx );
    elseif strcmp(nim.mods(imod).NLtype,'quad')
        %fgint = gint(:,imod).^2;
				fgint = fgint.^2;
    elseif strcmp(nim.mods(imod).NLtype,'lin')
        %fgint = gint(:,imod);
    elseif strcmp(nim.mods(imod).NLtype,'threshlin')
        %fgint = gint(:,imod);
        fgint(fgint < 0) = 0;
    else
        error('Invalid internal NL');
		end
    
    % Multiply by weight (and multiplier, if appl) and add to generating function
		if isempty(Gmults{imod})
			X(:,imod) = fgint*nim.mods(imod).sign;
		else
			X(:,imod) = (fgint.*Gmults{imod}) * nim.mods(imod).sign;
		end
end

%% PROCESS CONSTRAINTS
LB = []; UB = []; A = []; Aeq = []; % Initialize constraint parameters


%% HANDLE OPTIMIZATION PARAMETERS
%default optimization parameters
if silent == 1
    optim_params.Display = 'off';
else
    optim_params.Display = 'iter';
end
optim_params.MaxFunEvals = 400;
optim_params.MaxIter = 1e3;
optim_params.TolFun = 1e-6;
optim_params.TolX = 1e-6;
optim_params.HessUpdate = 'bfgs';
optim_params.GradObj = 'on';
optim_params.Algorithm = 'Active-set';

%load in specified optimization parameters
if ~isempty(desired_optim_params)
    spec_fields = fieldnames(desired_optim_params);
    for i = 1:length(spec_fields)
        optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
    end
end

%% Optimization 
[params] = fminunc( @(K) LLfit_scale_internal(nim, K, Robs, X), initial_params, optim_params);

if ~silent
	fprintf( 'Resulting gains:\n' )
	disp(sprintf('  %0.2f\n', params(1:end-1) ))
end

nim_out = nim;

% Set new values
nim_out.spk_NL_params(1) = params(end);

if spkhstlen 
	nim_out.spk_hist.coefs = nim_out.spk_hist.coefs * params(end-1);
end

for ii = 1:Nmods
    
    % Incorporate gains into filters
    if strcmp(nim.mods(ii).NLtype,'quad')
        nim_out.mods(ii).filtK = nim_out.mods(ii).filtK * sqrt(params(ii));
    else
        nim_out.mods(ii).filtK = nim_out.mods(ii).filtK * params(ii);
    end
    
	% If non-parmetric nonlinearity, need to scale x- and y- axes
	if strcmp(nim.mods(ii).NLtype,'nonpar')
		nim_out.mods(ii).NLx = nim_out.mods(ii).NLx * params(ii);
		nim_out.mods(ii).NLy = nim_out.mods(ii).NLy * params(ii);
	end
end

%% COMPUTE FINAL LL Values
[LL, penLL,~,G] = NMMmodel_eval( nim_out, Robs, Xstims, Gmults );
nim_out.LL_seq = cat(1,nim_out.LL_seq,LL);
nim_out.penLL_seq = cat(1,nim_out.penLL_seq,penLL);
nim_out.opt_history = cat(1,nim_out.opt_history,{'scale'});

end



%%%%%% INTERNAL FUNCTIONS %%%%%%

function [LL, LLgrad] = LLfit_scale_internal(nim, params, Robs, X )
%
% Internal function for computing LL and LLgradient with respect to the
% tent-basis coeffs

%% Useful params
NT = size(X,1);

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
G = params(end) + X*params(1:end-1);

%% Compute predicted firing rate
if strcmp(nim.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*nim.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    r = nim.spk_NL_params(4) + nim.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    r(too_large) = nim.spk_NL_params(4) + nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
elseif strcmp(nim.spk_NL_type,'exp')
    expg = exp(G);
    r = expg;
elseif strcmp(nim.spk_NL_type,'linear')
    r = G;
else
    error('invalid spk nl');
end
%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(r) < min_pred_rate
    r(r < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

%% COMPUTE LL and LL gradient

if strcmp(nim.spk_NL_type,'linear') % use MSE as cost function 
    LL = -sum( (Robs - r).^2 );
else
    LL = sum(Robs.* log(r) - r); %up to an overall constant
    %'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
end

%'residual' = dLL/dr
if strcmp(nim.spk_NL_type,'logexp')
    residual = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs./r - 1) .* expg ./ (1+expg);
    residual(too_large) = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs(too_large)./r(too_large) - 1);
elseif strcmp(nim.spk_NL_type,'exp')
    residual = Robs - r;
elseif strcmp(nim.spk_NL_type,'linear')
    residual = 2*(Robs - r);
else
    error('Unsupported spiking NL')
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

LLgrad(1:end-1) = residual'*X;

% Calculate derivatives with respect to constant term (theta)
LLgrad(end) = sum(residual);

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

end

