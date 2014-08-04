function nim_out = NIMfit_filters_EM(nim,Robs,Xstim,XLin,gamma,targets,silent,desired_optim_params)
%
% nim_out = NIMfit_filters(nim,Robs,Xstim,<XLin>,<targets>,<silent>,<desired_optim_params>)
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
%           (-1 = spk history filter, and -2 = extra linear filter) Default is to optimize all elements
%       <silent>: 0 to display optimization iterations, and 1 to suppress them
%       <desired_optim_params>: Struct of optimization parameters
% OUTPUTS:
%       nim_out: output model struct

%% PROCESS INPUTS
if nargin < 4
    XLin = [];
end
if nargin < 6
    targets = [];
end
if nargin < 7
    silent = 1;
end
if nargin < 8
    desired_optim_params = [];
end

%make sure Robs is a column vector
if size(Robs,2) > size(Robs,1)
    Robs = Robs';
end

[NT,filtLen] = size(Xstim); %stimulus dimensions
if filtLen ~= prod(nim.stim_params.stim_dims)
    error('Xstim dimensions dont match with stim_params')
end
nmods = length(nim.mods);

lin_dims = nim.stim_params.lin_dims;
if lin_dims ~= size(XLin,2);
    error('Mismatch between XLin and lin_dims');
end

spkhstlen = nim.spk_hist.spkhstlen;

if max(targets) > nmods %check input targets
    error('Invalid target specified');
end
if isempty(targets) %default is to optimize all model components
    targets = 1:nmods;
    if spkhstlen > 0
        targets = [targets -1]; %optimize spike hist filter
    end
    if lin_dims > 0
        targets = [targets -2]; %optimize additional linear filter
    end
end
Ntargets = sum(targets > 0); %number of targeted subunits
non_targets = setdiff([1:nmods -1 -2],targets); %elements of the model held constant

if ismember(-1,targets) && spkhstlen == 0
    error('No spk history term initialized')
end
if ismember(-2,targets) && lin_dims == 0
    error('No extra linear filter optimized')
end

%% PARSE INITIAL PARAMETERS
%compute initial fit parameters
initial_params = [];
for imod = targets(targets > 0)
    cur_kern = nim.mods(imod).filtK';
    initial_params = [initial_params; cur_kern']; %add coefs to initial param vector
end

%add in spike history coefs
if ismember(-1,targets)
    initial_params = [initial_params; nim.spk_hist.coefs];
end

%add in extra linear filter
if ismember(-2,targets)
    initial_params = [initial_params; nim.kLin];
end

initial_params(end+1) = nim.spk_NL_params(1); %add constant offset

%% COMPUTE L1 PENALTY IF APPLICABLE
lambda_L1 = zeros(size(initial_params));
for ii = 1:Ntargets
    cur_inds = (ii-1)*filtLen + (1:filtLen);
    lambda_L1(cur_inds) = nim.mods(targets(ii)).reg_params.lambda_L1';
end
lambda_L1 = lambda_L1/sum(Robs); %since we are dealing with LL/spk

%% PRECOMPUTE 'TENT-BASIS' DERIVATIVES OF UPSTREAM NLS IF NEEDED
if any(strcmp('nonpar',{nim.mods(targets(targets > 0)).NLtype}))
    for ii = 1:Ntargets
        if strcmp(nim.mods(targets(ii)).NLtype,'nonpar')
            NLx = nim.mods(targets(ii)).NLx;
            NL = nim.mods(targets(ii)).NLy;
            
            %compute derivative of non-linearity
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
    Xspkhst = create_spkhist_Xmat(Robs,nim.spk_hist.bin_edges);
else
    Xspkhst = [];
end

%% COMPUTE NET OUPTUT OF ALL NON-TARGET PREDICTORS
nt_gout = zeros(NT,1);
Kmat = [nim.mods(:).filtK]; %filter matrix
gint = Xstim*Kmat; %filter output of each subunit
for imod = non_targets(non_targets > 0) %for all subunits that aren't targeted
    %process subunit g's with upstream NLs
    if strcmp(nim.mods(imod).NLtype,'nonpar')
        fgint = piecelin_process(gint(:,imod),nim.mods(imod).NLy,nim.mods(imod).NLx);
    elseif strcmp(nim.mods(imod).NLtype,'quad')
        fgint = gint(:,imod).^2;
    elseif strcmp(nim.mods(imod).NLtype,'lin')
        fgint = gint(:,imod);
    elseif strcmp(nim.mods(imod).NLtype,'threshlin')
        fgint = gint(:,imod);
        fgint(fgint < 0) = 0;
    else
        error('Invalid internal NL');
    end
    
    %multiply by weight and add to generating function
    nt_gout = nt_gout + fgint*nim.mods(imod).sign;
end
if ismember(-1,non_targets) && spkhstlen > 0
    nt_gout = nt_gout + Xspkhst*nim.spk_hist.coefs(:);
end
if ismember(-2,non_targets) && lin_dims > 0
    nt_gout = nt_gout + XLin*nim.kLin(:);
end

%% IDENTIFY ANY CONSTRAINTS
use_con = 0;
LB = [];UB = []; A = []; Aeq = []; %initialize constraint matrices
if spkhstlen > 0 && ismember(-1,targets) %if optimizing spk history term
    %negative constraint on spk history coefs
    if nim.spk_hist.negCon == 1
        spkhist_inds = (Ntargets*filtLen + 1):(Ntargets*filtLen + spkhstlen);
        LB = -Inf*ones(size(initial_params));
        UB = Inf*ones(size(initial_params));
        UB(spkhist_inds) = 0;
        
        use_con = 1;
    end
end
beq = zeros(size(Aeq,1),1);
b = zeros(size(A,1),1);

%% GENERATE REGULARIZATION MATRICES
L2_mats = create_L2_matrices(nim);

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
                optim_params.verbose = 1;
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
            
            [params] = L1General2_PSSas(@(K) LLfit_filters_internal(nim, K, Robs, Xstim,Xspkhst,XLin,L2_mats,targets,nt_gout,fprimes),...
                initial_params,lambda_L1,optim_params);
%             [params] = L1General2_PSSas(@(K) NIM_fit_filters_internal(nim, K, Robs, Xstim,Xspkhst,XLin,L2_mats,targets,nt_gout,fprimes),...
%                 initial_params,lambda_L1);
        else
            error('Need to install Mark Schmidts L1 optimization toolbox for using L1');
        end
    else
        %if not using L1 reg
        if exist('minFunc','file') == 2 %try to use Mark Schmidt's optimizer
            
            %if using Mark Schmidt's optimization, some differences in option parameters
            optim_params.optTol = 1e-4;
            optim_params.progTol = 1e-8;
            optim_params.Method = 'lbfgs';
            %load in specified optimization parameters
            if ~isempty(desired_optim_params)
                spec_fields = fieldnames(desired_optim_params);
                for i = 1:length(spec_fields)
                    optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
                end
            end
            
            [params] = minFunc( @(K) LLfit_filters_internal(nim, K, Robs, Xstim,Xspkhst,XLin,L2_mats,targets,nt_gout,fprimes), initial_params, optim_params);
            
        else %if using Matlab Optim toolbox:
            %default optimization parameters
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
            
            [params] = fminunc( @(K) LLfit_filters_internal(nim, K, Robs, Xstim,Xspkhst,XLin,L2_mats,targets,nt_gout,fprimes), initial_params, optim_params);
        end
    end
else %if there are constraints
    
    %try to use Mark Schmidt constrained optimizer
    if exist('minConf_TMP','file') == 2 && isempty(A) && isempty(Aeq)
        %if using Mark Schmidt's optimization, some differences in option parameters
        optim_params.optTol = 1e-4;
        optim_params.progTol = 1e-6;
        [params] = minConf_TMP( @(K) LLfit_filters_internal(nim, K, Robs, Xstim,Xspkhst,XLin,L2_mats,targets,nt_gout,fprimes),...
            initial_params, LB,UB,optim_params);
    else
        %otherwise resort to matlab's
        optim_params.GradObj = 'on';
        optim_params.LargeScale = 'off';
        optim_params.Algorithm = 'active-set';
        optim_params.optTol = 1e-4;
        optim_params.progTol = 1e-6;
        [params] = fmincon( @(K) LLfit_filters_internal(nim, K, Robs, Xstim,Xspkhst,XLin,L2_mats,targets,nt_gout,fprimes),...
            initial_params, A,b,Aeq,beq,LB,UB,[],optim_params);
    end
end

%% PARSE MODEL FIT
nim_out = nim;
nim_out.spk_NL_params(1) = params(end);
if ismember(-1,targets)
    nim_out.spk_hist.coefs = params((Ntargets*filtLen+1):(Ntargets*filtLen+spkhstlen));
end
if ismember(-2,targets)
    nim_out.kLin = params((Ntargets*filtLen+spkhstlen+1):(Ntargets*filtLen+spkhstlen+lin_dims));
end

for ii = 1:Ntargets
    cur_kern = params((ii-1)*filtLen+(1:filtLen));
    nim_out.mods(targets(ii)).filtK = cur_kern(:);
end

[LL, penLL, ~, G, gint, fgint] = NIMmodel_eval(nim_out,Robs,Xstim,XLin);
nim_out.LL_seq = cat(1,nim_out.LL_seq,LL);
nim_out.penLL_seq = cat(1,nim_out.penLL_seq,penLL);
nim_out.opt_history = cat(1,nim_out.opt_history,{'filt'});

% Compute std dev of the output of each subunit
for n = 1:nmods
    mod_norm = std(fgint(:,n));
    nim_out.mods(n).mod_norm = mod_norm;
end
end

%%%% INTERNAL FUNCTIONS %%%%%%%

function [LL, LLgrad] = LLfit_filters_internal(nim, params, Robs, Xstim, Xspkhst, XLin, L2_mats, targets, nt_gout, fprimes)
%
% [LL, LLgrad] = LLfit_filters_internal(nim, params, Robs, Xstim, Xspkhst, XLin, L2_mats, targets, nt_gout, fprimes)
%
% Internal function for computing LL and LLgradient with respect to the
% stimulus filters

%% USEFUL PARAMETERS
Ntargets = sum(targets > 0);
filtLen = prod(nim.stim_params.stim_dims);
lin_dims = nim.stim_params.lin_dims;
spkhstlen = nim.spk_hist.spkhstlen;
nPix = squeeze(nim.stim_params.stim_dims(2:3));

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
theta = params(end); %offset
G = theta + nt_gout; %initialize overall generating function G
kmat = reshape(params(1:Ntargets*filtLen),filtLen,Ntargets);
gint = Xstim*kmat;
for ii = 1:Ntargets
    
    %process subunit g's with upstream NLs
    if strcmp(nim.mods(targets(ii)).NLtype,'nonpar')
        fgint = piecelin_process(gint(:,ii),nim.mods(targets(ii)).NLy,nim.mods(targets(ii)).NLx);
    elseif strcmp(nim.mods(targets(ii)).NLtype,'quad')
        fgint = gint(:,ii).^2;
    elseif strcmp(nim.mods(targets(ii)).NLtype,'lin')
        fgint = gint(:,ii);
    elseif strcmp(nim.mods(targets(ii)).NLtype,'threshlin')
        fgint = gint(:,ii);
        fgint(fgint < 0) = 0;
    else
        error('Invalid internal NL');
    end
    
    %multiply by weight and add to generating function
    G = G + fgint*nim.mods(targets(ii)).sign;
end

%add contribution from spike history filter
if spkhstlen > 0 && ismember(-1,targets)
    G = G + Xspkhst*params((Ntargets*filtLen + 1):(Ntargets*filtLen + spkhstlen));
end
%add contribution from linear filter
if lin_dims > 0 && ismember(-2,targets)
    G = G + XLin*params((Ntargets*filtLen+spkhstlen+1):(Ntargets*filtLen+spkhstlen+lin_dims));
end

%% Compute predicted firing rate
if strcmp(nim.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*nim.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    r = nim.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    r(too_large) = nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
elseif strcmp(nim.spk_NL_type,'exp')
    expg = exp(G);
    r = expg;
else
    error('invalid spk nl');
end
%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(r) < min_pred_rate
    r(r < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

%% COMPUTE LL and LL gradient
LL = sum(Robs.* log(r) - r); %up to an overall constant

%'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
if strcmp(nim.spk_NL_type,'logexp')
    residual = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs./r - 1) .* expg ./ (1+expg);
    residual(too_large) = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs(too_large)./r(too_large) - 1);
elseif strcmp(nim.spk_NL_type,'exp')
    residual = Robs - r;
else
    error('Unsupported spiking NL')
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant term (theta)
LLgrad(end) = sum(residual);

% Calculate derivative with respect to spk history filter
if spkhstlen > 0 && ismember(-1,targets)
    LLgrad((Ntargets*filtLen+1):(Ntargets*filtLen+spkhstlen)) = residual'*Xspkhst;
end
% Calculate derivative with respect to linear term
if lin_dims > 0 && ismember(-2,targets)
    LLgrad((Ntargets*filtLen+spkhstlen+1):(Ntargets*filtLen+spkhstlen+lin_dims)) = residual'*XLin;
end

%NOW COMPUTE LL grad WRT STIMULUS FILTERS
% Calculate output of derivative module
chunk_size = 1000; %maximum chunk size for processing high-dimensional stim filters
if filtLen <= chunk_size
    use_chunking = 0;
else
    use_chunking = 1;
    NChunks = ceil(filtLen/chunk_size); %divide stim filters into this many chunks for piecewise processing
end
for ii = 1:Ntargets
    
    if strcmp(nim.mods(targets(ii)).NLtype,'lin')
        LLgrad(((ii-1)*filtLen+1):(ii*filtLen)) = residual'*Xstim*nim.mods(targets(ii)).sign;
        
    else
        if strcmp(nim.mods(targets(ii)).NLtype,'nonpar')
            fpg = piececonst_process(gint(:,ii),fprimes{ii}, nim.mods(targets(ii)).NLx);
        elseif strcmp(nim.mods(targets(ii)).NLtype,'quad')
            fpg = 2*gint(:,ii);
        elseif strcmp(nim.mods(targets(ii)).NLtype,'threshlin')
            fpg = gint(:,ii) >= 0;
        else
            error('Unsupported NL type')
        end
        LLgrad(((ii-1)*filtLen+1):(ii*filtLen)) = (fpg.*residual)' * Xstim * nim.mods(targets(ii)).sign;
    end    
end

%% COMPUTE L2 PENALTIES AND ASSOCIATED CONTRIBUTIONS TO THE LL GRADIENT
smooth_penalty = zeros(Ntargets,1);
deriv_penalty = zeros(Ntargets,1);
ridge_penalty = zeros(Ntargets,1);

LLgrad_pen = zeros(size(LLgrad));
for ii = 1:Ntargets
    
    %temporal regularization
    if nim.mods(targets(ii)).reg_params.lambda_dT > 0
        deriv_penalty(ii) = deriv_penalty(ii) + nim.mods(targets(ii)).reg_params.lambda_dT*sum((L2_mats.L2_dT * kmat(:,ii)).^2);
        cur_grad_pen = 2*nim.mods(targets(ii)).reg_params.lambda_dT*(L2_mats.L2_dT' * L2_mats.L2_dT * kmat(:,ii));
        LLgrad_pen((ii-1)*filtLen + (1:filtLen)) = LLgrad_pen((ii-1)*filtLen + (1:filtLen)) + cur_grad_pen;
    end
    if nim.mods(targets(ii)).reg_params.lambda_d2T > 0
        smooth_penalty(ii) = smooth_penalty(ii) + nim.mods(targets(ii)).reg_params.lambda_d2T*sum((L2_mats.L2_d2T * kmat(:,ii)).^2);
        cur_grad_pen = 2*nim.mods(targets(ii)).reg_params.lambda_d2T*(L2_mats.L2_d2T' * L2_mats.L2_d2T * kmat(:,ii));
        LLgrad_pen((ii-1)*filtLen + (1:filtLen)) = LLgrad_pen((ii-1)*filtLen + (1:filtLen)) + cur_grad_pen;
    end
    
    %spatial (and spatiotemporal) regularization
    if nPix(1) > 1 %for spatial stimuli
%         if nPix(2) == 1 %for 1d stim
            if nim.mods(targets(ii)).reg_params.lambda_dX > 0
                deriv_penalty(ii) = deriv_penalty(ii) + nim.mods(targets(ii)).reg_params.lambda_dX*sum((L2_mats.L2_dX * kmat(:,ii)).^2);
                cur_grad_pen = 2*nim.mods(targets(ii)).reg_params.lambda_dX*(L2_mats.L2_dX' * L2_mats.L2_dX * kmat(:,ii));
                LLgrad_pen((ii-1)*filtLen + (1:filtLen)) = LLgrad_pen((ii-1)*filtLen + (1:filtLen)) + cur_grad_pen;
            end
            if nim.mods(targets(ii)).reg_params.lambda_d2X > 0
                smooth_penalty(ii) = smooth_penalty(ii) + nim.mods(targets(ii)).reg_params.lambda_d2X*sum((L2_mats.L2_d2X * kmat(:,ii)).^2);
                cur_grad_pen = 2*nim.mods(targets(ii)).reg_params.lambda_d2X*(L2_mats.L2_d2X' * L2_mats.L2_d2X * kmat(:,ii));
                LLgrad_pen((ii-1)*filtLen + (1:filtLen)) = LLgrad_pen((ii-1)*filtLen + (1:filtLen)) + cur_grad_pen;
            end
            if nim.mods(targets(ii)).reg_params.lambda_d2XT > 0
                smooth_penalty(ii) = smooth_penalty(ii) + nim.mods(targets(ii)).reg_params.lambda_d2XT*sum((L2_mats.L2_d2XT * kmat(:,ii)).^2);
                cur_grad_pen = 2*nim.mods(targets(ii)).reg_params.lambda_d2XT*(L2_mats.L2_d2XT' * L2_mats.L2_d2XT * kmat(:,ii));
                LLgrad_pen((ii-1)*filtLen + (1:filtLen)) = LLgrad_pen((ii-1)*filtLen + (1:filtLen)) + cur_grad_pen;
            end
%         else %2d stim
%             error('Only support for up to one spatial dimension now');
%         end
    end
    
    if nim.mods(targets(ii)).reg_params.lambda_L2 > 0
        ridge_penalty(ii) = nim.mods(targets(ii)).reg_params.lambda_L2*(kmat(:,ii)' * kmat(:,ii));
        LLgrad_pen((ii-1)*filtLen+(1:filtLen)) = LLgrad_pen((ii-1)*filtLen+(1:filtLen)) + nim.mods(targets(ii)).reg_params.lambda_L2*kmat(:,ii);
    end
end

LL = LL - sum(smooth_penalty) - sum(ridge_penalty) - sum(deriv_penalty);
LLgrad = LLgrad - LLgrad_pen;


%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

end
