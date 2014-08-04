function nim = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params,init_filts,init_kLin,spk_NL_type)
%
% nim = NIMinitialize_model(stim_params,mod_signs,NL_types,<reg_params>,<init_filts>,<init_kLin>,<spk_NL_type>)
%
% Creates a new NIM
% Inputs: 
%     stim_params: struct specifying needed stimulus parameters
%     mod_signs: Nmodsx1 vector specifying the 'sign' of each desired subunit (-1 or +1)
%     NL_types: Nmodsx1 cell array specifying the type of upstream NL associated with each subunit
%     <reg_def>: Struct specifying initial regularization parameters
%     <init_filts>: klenxNmods matrix specifying initial filter coefficients for all subunits
%     <init_kLin>: vector specifying an initial linear filter if needed
%     <spk_NL_type>: string specifying the type of spk NL (either 'exp' or 'logexp')
% OUTPUTS:
%     nim: Initialized model
% 
% NIM struct contains the fields:
%     LL_seq: Sequence of log-likelihood values (recorded after each stage of optimization)
%     penLL_seq: Same as LL_seq but for the penalized LL (or log-posterior)
%     opt_history: String array that specifies the sequence of parameter optimization stages applied
%     stim_params: Struct of stimulus parameters (see NIMcreate_stim_params for details).
%     spk_NL_type: String specifying the type of spiking NL ('logexp' or 'exp' currently supported)
%     mods: Array of structs, where each 'mod' describes a given nonlinear input (see below for more details)
%         filtK: Stimulus filter
%         sign: +1 for excitatory and -1 for suppressive
%         NLtype: String specifying the type of upstream NL ('lin','quad','threshlin',or 'nonpar')
%         mod_norm: std dev of the subunit output (in response to the stimulus it was most recently fit with)  
%         NLx: Vector of tent-basis center-locations when using a nonparametric upstream NL
%         NLy: Vector of tent-basis coefficients when using a nonparametric upstream NL.
%         reg_params: Struct of regularization parameters (see NIMcreate_reg_params)
%     kLin: Filter applied to any additional linear predictors (given in 'XLin')
%     spk_hist: Struct containing parameters of the spike history term
%         coefs: Coefficients of the spike-history filter
%         bin_edges: Bin edges (units of time lags) for the piecewise constant representation of spike history
%         spkhstlen: Number of spike-history filter coefficients
%     spk_NL_params: Vector of parameters for the spiking NL function



%%
nmods = length(mod_signs); %number of subunits

if nargin < 4 || isempty(reg_params)
    reg_params = NIMcreate_reg_params(); %use default reg_params
end
if nargin < 5 || isempty(init_filts)
    %use a random vector (this is mostly arbitrary, although initializing to zeros sometimes puts you at a local minimum)
    init_filts = randn(prod(stim_params.stim_dims),nmods)/prod(stim_params.stim_dims) * 0.1; 
end
if nargin < 6
    if stim_params.lin_dims == 0
        init_kLin = []; %no additional linear filter
    else
        init_kLin = zeros(stim_params.lin_dims,1);
    end
elseif length(init_kLin) ~= stim_params.lin_dims %check length of init_h
    error('Initial linear filter doesnt match stim_params');
end
if nargin < 7
    spk_NL_type = 'logexp'; %default to log(1+exp(x)) spiking NL
end

% COMPILE VECTOR OF SUBUNIT STRUCTS
poss_NL_types = {'lin','threshlin','quad','nonpar'};
mods = [];
for imod = 1:nmods
    tmod = []; 
    tmod.filtK = init_filts(:,imod); %store current filter coefs
    tmod.sign = mod_signs(imod); %sign of subunit (+1 for Exc, -1 for Sup)
    tmod.NLtype = NL_types{imod}; %upstream NL type
    tmod.mod_norm = nan;
    tmod.NLx = [];
    tmod.NLy = [];
    tmod.reg_params = reg_params;
    if ~strcmp(NL_types{imod},poss_NL_types)
        error('Unsupported NL type');
    end
    mods = [mods,tmod];
end

%initialize WITHOUT spike history term
spk_hist.coefs = []; 
spk_hist.bin_edges = [];
spk_hist.spkhstlen = 0;

if strcmp(spk_NL_type,'logexp')
    spk_NL_params = [0 1 1]; %[theta beta alpha]
elseif strcmp(spk_NL_type,'exp')
    spk_NL_params = [0]; %[theta]
    % elseif strcmp(spk_NL_type,'gauss')
    %     spk_nl_params = [0 1]; %[theta sigma
else
    error('Unsupported spiking NL type')
end

%create model STRUCT
nim = struct('LL_seq',[],'penLL_seq',[],'opt_history',[],'stim_params',stim_params,...
    'spk_NL_type',spk_NL_type,'mods',mods,'kLin',init_kLin,'spk_hist',...
    spk_hist,'spk_NL_params',spk_NL_params);





