function [glm0] = createGLM_lexp_rot(STCbvs,STCcf,init_signs,init_betas,defmod,basis,pix_conv_mat,kern_conv_mat)
%%USAGE: [glm0] = createGLM0(ks,defmod,mname)
%   takes the filters and default module structure to create a glm

if nargin < 7
    if strcmp(basis,'white')
        error('Need conv mats')
    end
    pix_conv_mat = [];
    kern_conv_mat = [];
end

[ndims,nmods] = size(STCcf);
mods = [];
for imod = 1:nmods
    tmod  = defmod; %initialize module param structure
    tmod.STCcf = STCcf(:,imod); %store filter coefs in STC space
    if strcmp(basis,'white')
        tmod.k = STCbvs*STCcf(:,imod); %store current filter coefs
        tmod.pix = tmod.k'*pix_conv_mat;
        tmod.pix = tmod.pix(:);
    elseif strcmp(basis,'pix')
        tmod.k = STCbvs*STCcf(:,imod); %store current filter coefs
        tmod.pix = tmod.k(:);
    end
    %if there is no psc term absorb sign of h into w
    if length(tmod.h) == 1
        tmod.w = init_signs(imod)*tmod.h;
        tmod.h = 1;
    else
        error('Cant handle PSC terms here');
    end
    
    tmod.nly = 1/init_betas(imod)*log(1+exp(init_betas(imod)*tmod.nlx)); 
    tmod.beta = init_betas(imod);
    tmod.nltype = 'lexp';
    mods = [mods,tmod];
    
end

glm0   = struct('mods',mods,'const',0,'basis',basis,'STCbasis',STCbvs,...
    'STCdim',ndims,'pix_conv_mat',pix_conv_mat,'kern_conv_mat',kern_conv_mat,...
    'LL',0,'LP',0,'lambdaW',0);

glm0.spk_nl = 'logexp';