function [gnm0] = createGNM_sepstim(init_kerns,init_signs,kern_types,kern_stim_inds,defmod,stim_params,spk_nl)

%init kerns: DxN matrix where D is stimulus dimension and N is number of
%filters
%kern type: "lexp","threshlin","lin","quad","rquad"

if nargin < 7
    spk_nl = 'logexp'; %default spiking NL
end

%set defaults
default.h = 1; %no psc term
default.lnl = 0;
default.lh = 0;
default.lnl2 = 0;
default.lh2 = 0;
default.nlcon = 0;
default.nlmon = 0;
default.locLambda = 0;
default.lambda_dX = 0;
default.lambda_L1x = 0;
default.lambda_dT = 0;
default.lambda_d2XT = 0;
default.lambda_d2X = 0;
default.lambda_L2x = 0;
default.lambda_L2x = 0;
default.lambda_d2T = 0;
default.kscale = 1;
defmod.nlx = linspace(-3.1,3.1,21); %internal NL x-axis
defmod.lexp_beta = 2; %if using lexp internal NLs
defmod.lexp_theta = 0; %if using lexp internal NLs

if ~isfield(stim_params,'spatial_dims')
    error('Must specify stimulus spatial dims')
end
if ~isfield(stim_params,'sdim')
    error('Must specify stimulus spatial dims (1-way)')
end
if ~isfield(stim_params,'flen')
    error('Must specify length of temporal kernels')
end
if stim_params.spatial_dims == 1
    stim_params.fsdim = stim_params.sdim;
elseif stim_params.spatial_dims == 2
    stim_params.fsdim = stim_params.sdim^2;
elseif stim_params.spatial_dims == 0
    stim_params.fsdim = stim_params.sdim;
else
    error('Unsupported stimulus dimensionality')
end

%load in specified values
use_defmod = default;
spec_fields = fieldnames(defmod);
for i = 1:length(spec_fields)
    use_defmod = setfield(use_defmod,spec_fields{i},getfield(defmod,spec_fields{i}));
end
defmod = use_defmod;

[ndims,nmods] = size(init_kerns);
mods = [];
for imod = 1:nmods
    tmod  = defmod; %initialize module param structure
    tmod.k = init_kerns(:,imod); %store current filter coefs
    tmod.w = init_signs(imod);
    
    if strcmp(kern_types{imod},'threshlin')
        tmod.nly = tmod.nlx;
        tmod.nly(tmod.nlx < 0) = 0;
        tmod.nltype = 'threshlin';
    elseif strcmp(kern_types{imod},'lin')
        tmod.nly = tmod.nlx;
        tmod.nltype = 'lin';
    elseif strcmp(kern_types{imod},'lexp')
        tmod.nly = 1/defmod.lexp_beta*log(1+exp(defmod.lexp_beta*(tmod.nlx-defmod.lexp_theta))) - ...
            1/defmod.lexp_beta*log(1+exp(defmod.lexp_beta*(tmod.nlx(1)-defmod.lexp_theta)));
        tmod.beta = defmod.lexp_beta;
        tmod.theta = defmod.lexp_theta;
        tmod.nltype = 'lexp';
    elseif strcmp(kern_types{imod},'quad')
        tmod.nly = tmod.nlx.^2;
        tmod.nltype = 'quad';
    elseif strcmp(kern_types{imod},'rquad')
        tmod.nly = tmod.nlx.^2;
        tmod.nly(tmod.nlx < 0) = 0;
        tmod.nltype = 'rquad';        
    else
        error('Unsupported NL type');
    end
    tmod.stim_ind = kern_stim_inds(imod);
    mods = [mods,tmod];  
end

% gnm0 = struct('mods',mods,'const',0,'LL',0,'LP',0,'stim_params',stim_params);
gnm0 = struct('mods',mods,'LL',0,'LP',0,'stim_params',stim_params);

gnm0.spk_nl = spk_nl;
gnm0.spk_alpha = 1;
gnm0.spk_theta = 0;
gnm0.spk_beta = 1;

%initialize likelihood vecs
gnm0.LP_seq = [];
% gnm.LP_ax = []; %iteration number
gnm0.LL_seq = [];
