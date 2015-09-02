function nim_out = NMMfit_filters( nim, Robs, Xstims, Gmults, Uindx, silent, desired_optim_params, regmat_custom,targets)
%
% Usage: nim_out = NMMfit_filters( nim, Robs, Xstims, <Gmults>, <Uindx>, <silent>, <desired_optim_params>, <regmat_custom>,<targets> )
%
% Optimizes the stimulus filters (plus extra linear terms if desired) for
% given upstream NLs
%
% INPUTS:
%       nim: model structure
%       Robs: binned spikes
%       Xstim: time-embedded stimulus mat
%       <XLin>: Matrix specifying additional linear predictors
%       <targets>: Vector of indices specifying which subunits to optimize.
%           (-1 = spk history filter, and -2 = offset only) Default is to optimize all elements
%       <silent>: 0 to display optimization iterations, and 1 to suppress them
%       <desired_optim_params>: Struct of optimization parameters
% OUTPUTS:
%       nim_out: output model struct

%% PROCESS INPUTS
Nmods = length(nim.mods);
if (nargin < 4) || (length(Gmults) < Nmods)
	Gmults{Nmods} = [];
end

if nargin < 5
	Uindx = [];
end

% Process Xstims (in case multiple Xstims)
if ~iscell(Xstims)
    tmp = Xstims;
    clear Xstims
    Xstims{1} = tmp;
end
if (nargin < 6) || isempty(silent)
	silent = 1;
end
if nargin < 7
	desired_optim_params = [];
end
if nargin < 8
	regmat_custom = [];
end
if nargin < 9
	targets = [];
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

for n = 1:Nmods
    [NT,filtLen] = size(Xstims{nim.mods(n).Xtarget}); %stimulus dimensions
    if filtLen ~= prod(nim.stim_params(nim.mods(n).Xtarget).stim_dims)
        error('Xstim dimensions dont match with stim_params')
    end
end

spkhstlen = nim.spk_hist.spkhstlen;

if max(targets) > Nmods %check input targets
	error('Invalid target specified');
end
if isempty(targets) %default is to optimize all model components
    targets = 1:Nmods;
    if spkhstlen > 0
        targets = [targets -1]; %optimize spike hist filter
    end
elseif targets == -2
    targets = [];
end

Ntargets = sum(targets > 0); %number of targeted subunits
non_targets = setdiff([1:Nmods -1 -2],targets); %elements of the model held constant

if ismember(-1,targets) && spkhstlen == 0
	error('No spk history term initialized')
end

% Add spk NL constant if it isnt already there
if length(nim.spk_NL_params) < 4
	nim.spk_NL_params(4) = 0;
end

%% PARSE INITIAL PARAMETERS
% Compute initial fit parameters
initial_params = [];
sign_con = [];
for imod = targets(targets > 0)
	cur_kern = nim.mods(imod).filtK';
    if isfield(nim.mods(imod),'Kcon')
	if nim.mods(imod).Kcon ~= 0
    sign_con(length(initial_params)+(1:length(cur_kern))) = nim.mods(imod).Kcon;
    end
  end
	initial_params = [initial_params; cur_kern']; % add coefs to initial param vector
end

% Add in spike history coefs
if ismember(-1,targets)
    initial_params = [initial_params; nim.spk_hist.coefs];
end

initial_params(end+1) = nim.spk_NL_params(1); % add constant offset
initial_params = initial_params(:);

%% COMPUTE L1 PENALTY IF APPLICABLE
lambda_L1 = zeros(size(initial_params));
cnt = 0;
for ii = 1:Ntargets
    filtLen = length(nim.mods(targets(ii)).filtK);
    cur_inds = (1:filtLen) + cnt;
    lambda_L1(cur_inds) = nim.mods(targets(ii)).reg_params.lambda_L1;
    cnt = cnt + filtLen;
end
lambda_L1 = lambda_L1/sum(Robs); % since we are dealing with LL/spk

%% PRECOMPUTE 'TENT-BASIS' DERIVATIVES OF UPSTREAM NLS IF NEEDED
if any(strcmp('nonpar',{nim.mods(targets(targets > 0)).NLtype}))
    for ii = 1:Ntargets
        if strcmp(nim.mods(targets(ii)).NLtype,'nonpar')
            NLx = nim.mods(targets(ii)).NLx;
            NL = nim.mods(targets(ii)).NLy;
            
            % Compute derivative of non-linearity
            fpr = zeros(1,length(NLx)-1);
            for n = 1:length(fpr)
                fpr(n) = (NL(n+1)-NL(n))/(NLx(n+1)-NLx(n));
            end
            fprimes{ii} = fpr;
        else
            fprimes{ii} = [];
        end
    end
else
    fprimes = [];
end

%% CREATE SPIKE HISTORY Xmat IF NEEDED
if nim.spk_hist.spkhstlen > 0
	Xspkhst = create_spkhist_Xmat( RobsFULL, nim.spk_hist.bin_edges );
  if ~isempty(Uindx)
    Xspkhst = Xspkhst(Uindx,:);
  end
else
	Xspkhst = [];
end

%% COMPUTE NET OUPTUT OF ALL NON-TARGET PREDICTORS
nt_gout = zeros(NT,1);
%Kmat = [nim.mods(:).filtK]; %filter matrix
%gint = Xstim*Kmat; %filter output of each subunit

for imod = non_targets(non_targets > 0) %for all subunits that aren't targeted
    
    fgint = Xstims{nim.mods(imod).Xtarget} * nim.mods(imod).filtK;
    
    % Process subunit g's with upstream NLs
    if strcmp(nim.mods(imod).NLtype,'nonpar')
        %fgint = piecelin_process(gint(:,imod),nim.mods(imod).NLy,nim.mods(imod).NLx);
        fgint = piecelin_process( fgint, nim.mods(imod).NLy, nim.mods(imod).NLx );
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
        nt_gout = nt_gout + fgint*nim.mods(imod).sign;
    else
        nt_gout = nt_gout + (fgint.*Gmults{imod}) * nim.mods(imod).sign;
    end
end

if ismember(-1,non_targets) && spkhstlen > 0
    nt_gout = nt_gout + Xspkhst*nim.spk_hist.coefs(:);
end

%% IDENTIFY ANY CONSTRAINTS
use_con = 0;
A = []; Aeq = []; %initialize constraint matrices
LB = -Inf*ones(size(initial_params));
UB = Inf*ones(size(initial_params));

% Constrain any of the filters to be positive or negative
if ~isempty(sign_con)
  LB(sign_con == 1) = 0;
  UB(sign_con == -1) = 0;
 
  use_con = 1;
end

if spkhstlen > 0 && ismember(-1,targets) %if optimizing spk history term
    %negative constraint on spk history coefs
    if nim.spk_hist.negCon == 1
        %spkhist_inds = (Ntargets*filtLen + 1):(Ntargets*filtLen + spkhstlen);
        spkhist_inds = cnt + (1:spkhstlen);
        UB(spkhist_inds) = 0;
        
        use_con = 1;
    end
end

beq = zeros(size(Aeq,1),1);
b = zeros(size(A,1),1);

%% GENERATE REGULARIZATION MATRICES
L2_mats = create_L2_matrices_NMM( nim );
L2_mats.custom = regmat_custom;

%%
if max(lambda_L1) > 0 && use_con == 1
    disp('Cant use L1 with constrained optimization, aborting constraints');
    use_con = 0;
end

%%
optim_params.MaxFunEvals = 100*length(initial_params);
optim_params.MaxIter = 1e3;
optim_params.Display = 'off';
if silent == 0
    optim_params.Display = 'iter';
end
if use_con == 0 %if no constraints
    
    %if using L1 reg
    if max(lambda_L1) > 0
        if exist('L1General2_PSSas','file') == 2
            optim_params.optTol = 1e-4;
            optim_params.progTol = 1e-8;
            if silent == 0
                optim_params.verbose = 2;
            else
                optim_params.verbose = 0;
            end
            % Load in specified optimization parameters
            if ~isempty(desired_optim_params)
                spec_fields = fieldnames(desired_optim_params);
                for i = 1:length(spec_fields)
                    optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
                end
            end
            
            [params] = L1General2_PSSas(@(K) LLfit_filters_internal(nim, K, Robs, Xstims,Xspkhst,Gmults,L2_mats,targets,nt_gout,fprimes),...
                initial_params,lambda_L1,optim_params);
            %             [params] = L1General2_PSSas(@(K) NIM_fit_filters_internal(nim, K, Robs, Xstim,Xspkhst,XLin,L2_mats,targets,nt_gout,fprimes),...
            %                 initial_params,lambda_L1);
        else
            error('Need to install Mark Schmidts L1 optimization toolbox for using L1');
        end
    else % if not using L1 reg
        
        if exist('minFunc','file') == 2 %try to use Mark Schmidt's optimizer
            % if using Mark Schmidt's optimization, some differences in option parameters
            optim_params.optTol = 1e-4;
            optim_params.progTol = 1e-8;
            optim_params.Method = 'lbfgs';
            if silent == 0
                optim_params.verbose = 2;
            else
                optim_params.verbose = 0;
            end
            %load in specified optimization parameters
            if ~isempty(desired_optim_params)
                spec_fields = fieldnames(desired_optim_params);
                for i = 1:length(spec_fields)
                    optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
                end
            end
            
            [params] = minFunc( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes ), initial_params, optim_params);
            
        else %if using Matlab Optim toolbox:
            
            % Default optimization parameters
            optim_params.LargeScale = 'off';
            optim_params.TolFun = 1e-6;
            optim_params.TolX = 1e-6;
            optim_params.HessUpdate = 'bfgs';
            optim_params.GradObj = 'on';
            
            %load in specified optimization parameters
            if ~isempty(desired_optim_params)
                spec_fields = fieldnames(desired_optim_params);
                for i = 1:length(spec_fields)
                    optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
                end
            end
            
            [params] = fminunc( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes ), initial_params, optim_params);
            
        end
    end
else %if there are constraints
    
    % Try to use Mark Schmidt constrained optimizer
    if exist('minConf_TMP','file') == 2 && isempty(A) && isempty(Aeq)
        % if using Mark Schmidt's optimization, some differences in option parameters
        optim_params.optTol = 1e-4;
        optim_params.progTol = 1e-6;
        if silent == 0
            optim_params.verbose = 2;
        else
            optim_params.verbose = 0;
        end
        [params] = minConf_TMP( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes ),...
            initial_params, LB,UB,optim_params);
    else
        % otherwise resort to matlab's
        optim_params.GradObj = 'on';
        optim_params.LargeScale = 'off';
        optim_params.Algorithm = 'active-set';
        optim_params.optTol = 1e-4;
        optim_params.progTol = 1e-6;
        [params] = fmincon( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes ),...
            initial_params, A,b,Aeq,beq,LB,UB,[],optim_params);
    end
    
end

%% PARSE MODEL FIT
nim_out = nim;
nim_out.spk_NL_params(1) = params(end);
if ismember(-1,targets)
  nim_out.spk_hist.coefs = params(cnt + (1:spkhstlen));
  %nim_out.spk_hist.coefs = params((Ntargets*filtLen+1):(Ntargets*filtLen+spkhstlen));
end
%if ismember(-2,targets)
%    nim_out.kLin = params((Ntargets*filtLen+spkhstlen+1):(Ntargets*filtLen+spkhstlen+lin_dims));
%end

cnt = 0;
for ii = 1:Ntargets
    filtLen = length(nim.mods(targets(ii)).filtK);
    cur_kern = params((1:filtLen) + cnt);
    nim_out.mods(targets(ii)).filtK = cur_kern(:);
    cnt = cnt + filtLen;
end

[LL, nullLL, ~, G, gint, fgint, penLL] = NMMeval_model( nim_out, Robs, Xstims, Gmults, [], regmat_custom );
nim_out.LL_seq = cat(1,nim_out.LL_seq,LL);
nim_out.penLL_seq = cat(1,nim_out.penLL_seq,penLL);
nim_out.opt_history = cat(1,nim_out.opt_history,{'filt'});

% Compute std dev of the output of each subunit
for n = 1:Nmods
    mod_norm = std(fgint(:,n));
    nim_out.mods(n).mod_norm = mod_norm;
end
end

%%%% INTERNAL FUNCTIONS %%%%%%%

function [LL, LLgrad] = LLfit_filters_internal(nim, params, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes)
%
% [LL, LLgrad] = LLfit_filters_internal(nim, params, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes)
%
% Internal function for computing LL and LLgradient with respect to the
% stimulus filters

%% USEFUL PARAMETERS
Ntargets = sum(targets > 0);
%lin_dims = nim.stim_params.lin_dims;
spkhstlen = nim.spk_hist.spkhstlen;

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
theta = params(end); % offset
G = theta + nt_gout; % initialize overall generating function G

%kmat = reshape(params(1:Ntargets*filtLen),filtLen,Ntargets);
%gint = Xstim*kmat;

gint = nan(length(Robs),Ntargets);

filtLen = zeros(Ntargets,1);  ks = cell(Ntargets,1);

% NKtot = 0;  
% for ii = 1:Ntargets
%     
%     tar = targets(ii);
%     
%     % Pull out (potentially different-sized) filters from params
%     filtLen(ii) = prod(nim.stim_params(nim.mods(tar).Xtarget).stim_dims);
%     ks{ii} = params(NKtot+(1:filtLen(ii)));
%     NKtot = NKtot + filtLen(ii);
%     
%     gint(:,ii) = Xstims{nim.mods(tar).Xtarget} * ks{ii};
%     
%     % Process subunit g's with upstream NLs
%     if strcmp(nim.mods(tar).NLtype,'nonpar')
%         fgint = piecelin_process(gint(:,ii),nim.mods(tar).NLy,nim.mods(tar).NLx);
%     elseif strcmp(nim.mods(tar).NLtype,'quad')
%         fgint = gint(:,ii).^2;
%     elseif strcmp(nim.mods(tar).NLtype,'lin')
%         fgint = gint(:,ii);
%     elseif strcmp(nim.mods(tar).NLtype,'threshlin')
%         fgint = gint(:,ii);
%         fgint(fgint < 0) = 0;
%     else
%         error('Invalid internal NL');
%     end
%     
%     % Multiply by weight (and multiplier, if appl) and add to generating function
%     if isempty(Gmults{tar})
%         G = G + fgint*nim.mods(tar).sign;
%     else
%         G = G + (fgint.*Gmults{tar}) * nim.mods(tar).sign;
%     end
%     
% end

Xtarg_set = [nim.mods(targets).Xtarget]; %vector of Xtargets for set of subunits being optimized
un_Xtargs = unique(Xtarg_set); %set of unique Xtargets
nXtargs = length(un_Xtargs); %number of uniuqe Xtargets
param_inds = cell(Ntargets,1); %this will store the index values of each subunit's filter coefs within the parameter vector
NKtot = 0;  %init counter
for jj = 1:Ntargets %loop over subunits that are being optimized
    filtLen(jj) = prod(nim.stim_params(nim.mods(targets(jj)).Xtarget).stim_dims); % Pull out (potentially different-sized) filters from params
    param_inds{jj} = NKtot + (1:filtLen(jj)); %set of param indices associated with this subunit's filters
    ks{jj} = params(param_inds{jj}); %store filter coefs
    NKtot = NKtot + filtLen(jj); %inc counter
end
for ii = 1:nXtargs %now loop over the unique Xtargets and compute the generating signals 
    cur_mods = find(Xtarg_set == un_Xtargs(ii)); %set of targeted subunits that act on this Xtarg
    gint(:,cur_mods) = Xstims{un_Xtargs(ii)} * cat(2,ks{cur_mods}); %get filter outputs
end   
mod_NL_types = {nim.mods(targets).NLtype}; %NL types for each targeted subunit
unique_NL_types = unique(mod_NL_types); %unique set of NL types being used
fgint = gint; %init subunit outputs by filter outputs
for ii = 1:length(unique_NL_types) %loop over unique subunit NL types
    cur_mods = find(strcmp(mod_NL_types,unique_NL_types{ii}));
    % Process subunit g's with upstream NLs
    if strcmp(unique_NL_types{ii},'nonpar')
        for jj = 1:length(cur_mods) %for nonpar NLs, loop over individual subunits and compute outputs
            fgint(:,cur_mods(jj)) = piecelin_process(gint(:,cur_mods(jj)),...
                nim.mods(targets(cur_mods(jj))).NLy,nim.mods(targets(cur_mods(jj))).NLx);
        end
    elseif strcmp(unique_NL_types{ii},'quad')
        fgint(:,cur_mods) = fgint(:,cur_mods).^2;
    elseif strcmp(unique_NL_types{ii},'lin')
        %do nothing
    elseif strcmp(unique_NL_types{ii},'threshlin')
        fgint(fgint(:,cur_mods) < 0,cur_mods) = 0; %threshold at 0
    else
        error('Invalid internal NL');
    end
end
 
mod_signs = [nim.mods(targets).sign]'; %signs of targeted subunits

% Multiply by weight (and multiplier, if appl) and add to generating function
if all(cellfun(@(x) isempty(x),Gmults))
    G = G + fgint*mod_signs;
else
   error('not implemented yet')
end

% Add contribution from spike history filter
if spkhstlen > 0 && ismember(-1,targets)
    G = G + Xspkhst*params(NKtot + (1:spkhstlen));
end

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

% Enforce minimum predicted firing rate to avoid nan LLs
if ~strcmp(nim.spk_NL_type,'linear')
    min_pred_rate = 1e-50;
    if min(r) < min_pred_rate
        r(r < min_pred_rate) = min_pred_rate; %minimum predicted rate
    end
end

%% COMPUTE LL and LL gradient
if strcmp(nim.spk_NL_type,'linear') % use MSE as cost function 
    Nspks = length(Robs);
    LL = -sum( (Robs - r).^2 );
else
    Nspks = sum(Robs);
    LL = sum(Robs.* log(r) - r); %up to an overall constant
    %'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
end

%'residual' = LL'*F'[] where F[.] is the spk NL and LL' is the derivative
%of the LL cost
if strcmp(nim.spk_NL_type,'logexp')
    residual = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs./r - 1) .* expg ./ (1+expg);
    residual(too_large) = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs(too_large)./r(too_large) - 1);
    residual(r == min_pred_rate) = 0; %for points beneath the lower threshold of the spk NL, take F'[.] = 0
elseif strcmp(nim.spk_NL_type,'exp')
    residual = Robs - r;
    residual(r == min_pred_rate) = 0;
elseif strcmp(nim.spk_NL_type,'linear')
    residual = 2*(Robs - r);    
else
    error('Unsupported spiking NL')
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant term (theta)
LLgrad(end) = sum(residual);

% Calculate derivative with respect to spk history filter
if spkhstlen > 0 && ismember(-1,targets)
    LLgrad(NKtot+(1:spkhstlen)) = residual'*Xspkhst;
end

%NOW COMPUTE LL grad WRT STIMULUS FILTERS

% placeholder = 0;
% for ii = 1:Ntargets
%     tar = targets(ii);
%     if strcmp(nim.mods(tar).NLtype,'lin')
%         % Check for multiplicative interactions
%         if isempty(Gmults{tar})
%             LLgrad(placeholder+(1:filtLen(ii))) = residual'*Xstims{nim.mods(tar).Xtarget} * nim.mods(tar).sign;
%         else
%             LLgrad(placeholder+(1:filtLen(ii))) = (Gmults{tar}.*residual')* Xstims{nim.mods(tar).Xtarget} * nim.mods(tar).sign;
%         end
%         
%     else
%         if strcmp(nim.mods(tar).NLtype,'nonpar')
%             fpg = piececonst_process(gint(:,ii),fprimes{ii}, nim.mods(tar).NLx);
%         elseif strcmp(nim.mods(tar).NLtype,'quad')
%             fpg = 2*gint(:,ii);
%         elseif strcmp(nim.mods(tar).NLtype,'threshlin')
%             fpg = gint(:,ii) >= 0;
%             
%         else
%             error('Unsupported NL type')
%         end
%         % Add for multiplicative interactions
%         if ~isempty(Gmults{tar})
%             fpg = fpg .* Gmults{tar};
%         end
%         LLgrad(placeholder+(1:filtLen(ii))) = (fpg.*residual)' * Xstims{nim.mods(tar).Xtarget} * nim.mods(tar).sign;
%     end
%     placeholder = placeholder + filtLen(ii);
% end

for ii = 1:un_Xtargs %loop over unique Xtargets
    cur_mod_inds = find(Xtarg_set == un_Xtargs(ii)); %set of subunits with this Xtarget
    cur_NL_types = mod_NL_types(cur_mod_inds); %current NL types
    cur_unique_NL_types = unique(cur_NL_types); %set of unique NL types
    
    if length(cur_mod_inds) == 1 && strcmp(cur_unique_NL_types,'lin') %if there's only a single linear subunit, faster calc
        LLgrad(param_inds{cur_mod_inds}) = residual'*Xstims{un_Xtargs(ii)} * nim.mods(cur_mod_inds).sign;
    else %otherwise, compute a matrix of derivatives fpg
        fpg = ones(length(residual),length(cur_mod_inds)); %initialize to linear NL derivative (all ones)
        for jj = 1:length(cur_unique_NL_types) %loop over unique NL types
            cur_mod_subinds = find(strcmp(cur_NL_types,cur_unique_NL_types{jj})); %indices of current subset of subunits
            if strcmp(cur_unique_NL_types{jj},'quad')
                fpg(:,cur_mod_subinds) = 2*gint(:,cur_mod_inds(cur_mod_subinds));
            elseif strcmp(cur_unique_NL_types{jj},'threshlin')
                fpg(:,cur_mod_subinds) = gint(:,cur_mod_inds(cur_mod_subinds)) >= 0;
            elseif strcmp(cur_unique_NL_types{jj},'nonpar') %if using nonpara NLs, loop over inidividual subunits to get derivatives
                for zz = 1:length(cur_mod_subinds)
                    fpg(:,cur_mod_subinds(zz)) = piececonst_process(gint(:,cur_mod_inds(cur_mod_subinds(zz))),...
                        fprimes{cur_mod_inds(cur_mod_subinds(zz))},nim.mods(cur_mod_inds(cur_mod_subinds(zz))).NLx);
                end
            end
        end
        target_params = cat(2,param_inds{cur_mod_inds}); %indices of filter coefs for current set of targeted subunits
        
        %LL grad is residual * f'(.) *X *w, computed in parallel for all
        %subunits targeting this Xtarg
        LLgrad(target_params) = bsxfun(@times,(bsxfun(@times,fpg,residual)'*Xstims{un_Xtargs(ii)}),mod_signs(cur_mod_inds))';
    end
end

%% COMPUTE L2 PENALTIES AND ASSOCIATED CONTRIBUTIONS TO THE LL GRADIENT
smooth_penalty = zeros(Ntargets,1);
deriv_penalty = zeros(Ntargets,1);
ridge_penalty = zeros(Ntargets,1);
custom_penalty = zeros(Ntargets,1);

LLgrad_pen = zeros(size(LLgrad));
placeholder = 0;
for ii = 1:Ntargets
    
    tar = targets(ii);
    
    % Temporal regularization
    if nim.mods(tar).reg_params.lambda_dT > 0
        deriv_penalty(ii) = deriv_penalty(ii) + nim.mods(tar).reg_params.lambda_dT*sum((L2_mats.L2_dT{nim.mods(tar).Xtarget} * ks{ii}).^2);
        cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_dT*(L2_mats.L2_dT{nim.mods(tar).Xtarget}' * L2_mats.L2_dT{nim.mods(tar).Xtarget} * ks{ii});
        LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
    end
    if nim.mods(tar).reg_params.lambda_d2T > 0
        smooth_penalty(ii) = smooth_penalty(ii) + nim.mods(tar).reg_params.lambda_d2T*sum((L2_mats.L2_d2T{nim.mods(tar).Xtarget} * ks{ii}).^2);
        cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_d2T*(L2_mats.L2_d2T{nim.mods(tar).Xtarget}' * L2_mats.L2_d2T{nim.mods(tar).Xtarget} * ks{ii});
        LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
    end
    
    % Spatial (and spatiotemporal) regularization
    if nim.stim_params(nim.mods(tar).Xtarget).stim_dims(2) > 1
        % 	if nPix(1) > 1 % for spatial stimuli
        
        if nim.mods(tar).reg_params.lambda_dX > 0
            deriv_penalty(ii) = deriv_penalty(ii) + nim.mods(tar).reg_params.lambda_dX*sum((L2_mats.L2_dX{nim.mods(tar).Xtarget} * ks{ii}).^2);
            cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_dX*(L2_mats.L2_dX{nim.mods(tar).Xtarget}' * L2_mats.L2_dX{nim.mods(tar).Xtarget} * ks{ii});
            LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
        end
        if nim.mods(tar).reg_params.lambda_d2X > 0
            smooth_penalty(ii) = smooth_penalty(ii) + nim.mods(tar).reg_params.lambda_d2X*sum((L2_mats.L2_d2X{nim.mods(tar).Xtarget} * ks{ii}).^2);
            cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_d2X*(L2_mats.L2_d2X{nim.mods(tar).Xtarget}' * L2_mats.L2_d2X{nim.mods(tar).Xtarget} * ks{ii});
            LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
        end
        if nim.mods(tar).reg_params.lambda_d2XT > 0
            smooth_penalty(ii) = smooth_penalty(ii) + nim.mods(tar).reg_params.lambda_d2XT*sum((L2_mats.L2_d2XT{nim.mods(tar).Xtarget} * ks{ii}).^2);
            cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_d2XT*(L2_mats.L2_d2XT{nim.mods(tar).Xtarget}' * L2_mats.L2_d2XT{nim.mods(tar).Xtarget} * ks{ii});
            LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
        end
    end
    
    %for custom regularization
    if nim.mods(tar).reg_params.lambda_custom > 0
        custom_penalty(ii) = custom_penalty(ii) + nim.mods(tar).reg_params.lambda_custom*sum((L2_mats.custom * ks{ii}).^2);
        cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_custom*(L2_mats.custom' * L2_mats.custom * ks{ii});
        LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
    end
    
    if nim.mods(tar).reg_params.lambda_L2 > 0
        ridge_penalty(ii) = nim.mods(tar).reg_params.lambda_L2*(ks{ii}' * ks{ii});
        LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + 2*nim.mods(tar).reg_params.lambda_L2*ks{ii};
    end
    
    placeholder = placeholder + filtLen(ii);
end

LL = LL - sum(smooth_penalty) - sum(ridge_penalty) - sum(deriv_penalty) - sum(custom_penalty);
LLgrad = LLgrad - LLgrad_pen;


%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS

LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

end
