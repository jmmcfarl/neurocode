
function [glm0] = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,pix_conv_mat,kern_conv_mat)
%%USAGE: [glm0] = createGLM0(ks,defmod,mname)
%   takes the filters and default module structure to create a glm

%set defaults
default.h = 1;
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

%load in specified values
use_defmod = default;
spec_fields = fieldnames(defmod);
for i = 1:length(spec_fields)
    use_defmod = setfield(use_defmod,spec_fields{i},getfield(defmod,spec_fields{i}));
end
defmod = use_defmod;

%default all lexp internal NLs
defmod.nltype = 'lexp';

if nargin < 7
    if strcmp(basis,'white')
        error('Need conv mats')
    end
    pix_conv_mat = [];
    kern_conv_mat = [];
end

[ndims,nmods] = size(init_kerns);
mods = [];
for imod = 1:nmods
    tmod  = defmod; %initialize module param structure
    if strcmp(basis,'white')
        tmod.k = init_kerns(:,imod); %store current filter coefs
        tmod.pix = tmod.k'*pix_conv_mat;
        tmod.pix = tmod.pix(:);
    elseif strcmp(basis,'pix')
        tmod.k = init_kerns(:,imod); %store current filter coefs
        tmod.pix = tmod.k(:);
    end
    %if there is no psc term absorb sign of h into w
    if length(tmod.h) == 1
        tmod.w = init_signs(imod)*tmod.h;
        tmod.h = 1;
    else
        error('Cant handle PSC terms here');
    end
    
    if strcmp(tmod.nltype,'lexp')
        tmod.nly = 1/init_betas(imod)*log(1+exp(init_betas(imod)*(tmod.nlx-init_thetas(imod)))) - ...
           1/init_betas(imod)*log(1+exp(init_betas(imod)*(tmod.nlx(1)-init_thetas(imod))));
        tmod.beta = init_betas(imod);
        tmod.theta = init_thetas(imod);
    elseif strcmp(tmod.nltype,'quad')
        tmod.nly = tmod.nlx.^2;
    end
    
    mods = [mods,tmod];
    
end

glm0   = struct('mods',mods,'const',0,'basis',basis,...
    'pix_conv_mat',pix_conv_mat,'kern_conv_mat',kern_conv_mat,...
    'LL',0,'LP',0,'lambdaW',0);

glm0.spk_nl = 'logexp';