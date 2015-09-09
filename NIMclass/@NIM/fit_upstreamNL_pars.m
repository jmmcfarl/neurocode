function nim = fit_upstreamNL_pars(nim, Robs, Xstims, varargin)
%         nim = nim.fit_upstreamNL_pars(Robs, Xstims, <train_inds>, varargin)
%         Optimizes the parameters of the spkNL
%         INPUTS:
%            Robs: vector of response observations (e.g. spike counts)
%            Xstims: cell array of stimuli
%            <train_inds>: index values of data on which to fit the model [default to all indices in provided data]
%            optional flags:
%                ('gain_funs',gain_funs): matrix of multiplicative factors, one column for each subunit
%                ('optim_params',optim_params): struct of desired optimization parameters
%                'silent': include this flag to suppress the iterative optimization display
%                'hold_const': vector of parameter indices to hold constant
%        OUTPUTS:
%            nim: output model struct

% PROCESS INPUTS
Nsubs = length(nim.subunits); %number of subunits
fit_subs = 1:Nsubs; %defualt to fitting all subunits (plus -1 for spkHist filter)
gain_funs = []; %default has no gain_funs
train_inds = nan; %default nan means train on all data
optim_params = []; %default has no user-specified optimization parameters
fit_spk_hist = nim.spk_hist.spkhstlen > 0; %default is to fit the spkNL filter if it exists
silent = false; %default is to display optimization output
j = 1;
while j <= length(varargin)
    flag_name = varargin{j};
    if ~ischar(flag_name)%if not a flag, it must be train_inds
        train_inds = flag_name;
        j = j + 1; %only one argument here
    else
        switch lower(flag_name)
            case 'sub_inds'
                fit_subs = varargin{j+1};
                assert(all(ismember(fit_subs,1:Nsubs)),'invalid target subunits specified');
                j = j + 2;
            case 'gain_funs'
                gain_funs = varargin{j+1};
                j = j + 2;
            case 'optim_params'
                optim_params = varargin{j+1};
                assert(isstruct(optim_params),'optim_params must be a struct');
                j = j + 2;
            case 'silent'
                silent = true;
                j = j + 1;
            case 'hold_spkhist'
                fit_spk_hist = false;
                j = j + 1;
            otherwise
                error('Invalid input flag');
        end
    end
end
if size(Robs,2) > size(Robs,1); Robs = Robs'; end; %make Robs a column vector
nim.check_inputs(Robs,Xstims,train_inds,gain_funs); %make sure input format is correct

Nfit_subs = length(fit_subs); %number of targeted subunits
non_fit_subs = setdiff([1:Nsubs],fit_subs); %elements of the model held constant
spkhstlen = nim.spk_hist.spkhstlen; %length of spike history filter
if fit_spk_hist; assert(spkhstlen > 0,'no spike history term initialized!'); end;
if spkhstlen > 0 % create spike history Xmat IF NEEDED
    Xspkhst = create_spkhist_Xmat( Robs, nim.spk_hist.bin_edges);
else
    Xspkhst = [];
end
if ~isnan(train_inds) %if specifying a subset of indices to train model params
    for nn = 1:length(Xstims)
        Xstims{nn} = Xstims{nn}(train_inds,:); %grab the subset of indices for each stimulus element
    end
    Robs = Robs(Uindx);
    if ~isempty(Xspkhst); Xspkhst = Xspkhst(train_inds,:); end;
    if ~isempty(gain_funs); gain_funs = gain_funs(train_inds,:); end;
end

% PARSE INITIAL PARAMETERS
init_params = [];
for imod = fit_subs
    cur_params = nim.subunits(imod).NLparams;
    init_params = [init_params; cur_params]; % add coefs to initial param vector
end
Nfit_filt_params = length(init_params); %number of filter coefficients in param vector
% Add in spike history coefs
if fit_spk_hist
    init_params = [init_params; nim.spk_hist.coefs];
end

init_params(end+1) = nim.spkNL.theta; % add constant offset
[nontarg_g] = nim.process_stimulus(Xstims,non_fit_subs,gain_funs);
if ~fit_spk_hist && spkhstlen > 0 %add in spike history filter output, if we're not fitting it
    nontarg_g = nontarg_g + Xspkhst*nim.spk_hist.coefs(:);
end

[~, ~, gint] = nim.process_stimulus(Xstims,fit_subs,gain_funs);

%BOUND CONSTRAINTS
LB = -Inf*ones(size(init_params));
UB = Inf*ones(size(init_params));
use_con = 0;
%CONSTRAINT INIT

if use_con
optimizer = 'fmincon';
else
    optimizer = 'minFunc';
end
fit_opts = struct('fit_spk_hist', fit_spk_hist, 'fit_subs',fit_subs); %put any additional fitting options into this struct
%the function we want to optimize
opt_fun = @(K) internal_upstreamNL_par(nim,K,Robs,gint,Xspkhst,nontarg_g,gain_funs,fit_opts);
init_params = init_params';
optim_params = nim.set_optim_params(optimizer,optim_params,silent);
if ~silent; fprintf('Running optimization using %s\n\n',optimizer); end;
switch optimizer %run optimization
    case 'minFunc'
        [params] = minFunc(opt_fun, init_params, optim_params);
    case 'fminunc'
        [params] = fminunc(opt_fun, init_params, optim_params);
    case 'minConf_TMP'
        [params] = minConf_TMP(opt_fun, init_params, LB, UB, optim_params);
    case 'fmincon'
        [params] = fmincon(opt_fun, init_params, [], [], [], [], LB, UB, [], optim_params);
end

% opt_fun = @(K) internal_LL_spkNL(nim,K, Robs, G);
% params = fmincon(opt_fun, init_params, [], [], Aeq, Beq, LB, UB, [], optim_params);
[~,penGrad] = opt_fun(params);
% first_order_optim = max(abs(penGrad));
% if first_order_optim > nim.opt_check_FO
%     warning(sprintf('First-order optimality: %.3f, fit might not be converged!',first_order_optim));
% end
% 
% nim.spkNL.params = params(1:end-1);
% nim.spkNL.theta = params(end);
% 
% [LL,~,mod_internals,LL_data] = nim.eval_model(Robs,Xstims,'gain_funs',gain_funs);
% nim = nim.set_subunit_scales(mod_internals.fgint); %update filter scales
% cur_fit_details = struct('fit_type','spkNL','LL',LL,'filt_pen',LL_data.filt_pen,...
%     'NL_pen',LL_data.NL_pen,'FO_optim',first_order_optim);
% nim.fit_props = cur_fit_details;
% nim.fit_hist = cat(1,nim.fit_hist,cur_fit_details);
end

%%
function [LL, LLgrad] = internal_upstreamNL_par(nim,params,Robs,gint,Xspkhst,nontarg_g,gain_funs,fit_opts)
%computes the LL and its gradient for given set of spkNL parameters

fit_subs = fit_opts.fit_subs;
Nfit_subs = length(fit_subs); %number of targeted subs
theta = params(end); % offset
mod_weights = [nim.subunits(fit_subs).weight];

G = theta + nontarg_g; % initialize overall generating function G with the offset term and the contribution from nontarget subs

pcnt = 0;
fgint = nan(length(Robs),Nfit_subs);
for ii = 1:Nfit_subs %
    nparams = length(nim.subunits(fit_subs(ii)).NLparams);
    cur_params = params((1:nparams) + pcnt);
    pcnt = pcnt + nparams;
    nim.subunits(fit_subs(ii)).NLparams = cur_params;
    fgint(:,ii) = nim.subunits(fit_subs(ii)).apply_NL(gint(:,ii));
end

% Multiply by weight (and multiplier, if appl) and add to generating function
if isempty(gain_funs)
    G = G + fgint*mod_weights;
else
    G = G + (fgint.*gain_funs(:,fit_subs))*mod_weights;
end

% Add contribution from spike history filter
if fit_opts.fit_spk_hist
    G = G + Xspkhst*params(pcnt + (1:nim.spk_hist.spkhstlen));
end

pred_rate = nim.apply_spkNL(G);
LL = nim.internal_LL(pred_rate,Robs); %compute LL

%residual = LL'[r].*F'[g]
residual = nim.internal_LL_deriv(pred_rate,Robs) .* nim.apply_spkNL_deriv(G,pred_rate <= nim.min_pred_rate);

LLgrad = zeros(length(params),1); %initialize LL gradient
LLgrad(end) = sum(residual);      %Calculate derivatives with respect to constant term (theta)

pcnt = 0;
for ii = 1:Nfit_subs
    nparams = length(nim.subunits(fit_subs(ii)).NLparams);
    cur_set = (1:nparams) + pcnt;
    param_deriv = -(gint(:,ii) > params(cur_set(1)));
    LLgrad(cur_set) = residual'*param_deriv;
end

% Calculate derivative with respect to spk history filter
if fit_opts.fit_spk_hist
    LLgrad(pcnt+(1:nim.spk_hist.spkhstlen)) = residual'*Xspkhst;
end

% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

end

