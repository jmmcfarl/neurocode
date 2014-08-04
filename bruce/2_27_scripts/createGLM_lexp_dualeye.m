function [glm0] = createGLM_lexp_dualeye(init_kerns,init_signs,init_betas,defmod,basis,pix_conv_mat,kern_conv_mat)
%%USAGE: [glm0] = createGLM0(ks,defmod,mname)
%   takes the filters and default module structure to create a glm

%default all lexp internal NLs
defmod.nltype = 'lexp';

if nargin < 6
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
        tmod.nly = 1/init_betas(imod)*log(1+exp(init_betas(imod)*tmod.nlx)) - ...
            1/init_betas(imod)*log(1+exp(init_betas(imod)*tmod.nlx(1)));
        tmod.beta = init_betas(imod);
    elseif strcmp(tmod.nltype,'quad')
        tmod.nly = tmod.nlx.^2;
    end
    
    mods = [mods,tmod];
    
end

klen = size(init_kerns,1);
eyeinds = ones(klen,1);
eyeinds(floor(klen/2)+1:end) = 2;
glm0   = struct('mods',mods,'const',0,'basis',basis,...
    'pix_conv_mat',pix_conv_mat,'kern_conv_mat',kern_conv_mat,...
    'LL',0,'LP',0,'lambdaW',0,'eyeinds',eyeinds);

glm0.spk_nl = 'logexp';