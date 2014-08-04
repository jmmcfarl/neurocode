function [glm0] = createGLM1d_fullbf(STCbvs,STCcf_0,defmod,nltype,init_nls,basis,mname)
%%USAGE: [glm0] = createGLM0(ks,defmod,mname)
%   takes the filters and default module structure to create a glm


[nstc_dim,nmods] = size(STCcf_0);
mods = [];
for imod = 1:nmods
    tmod  = defmod; %initialize module param structure
    tmod.k = STCbvs*STCcf_0(:,imod); %store current filter coefs
    tmod.k = tmod.k(:);
    %if there is no psc term absorb sign of h into w
    if length(tmod.h) == 1
        tmod.w = tmod.h;
        tmod.h = 1;
    end
    
    if strcmp(init_nls{imod},'tl') %threshold linear
        tmod.nly = tmod.nlx;
        tmod.nly(tmod.nlx <= 0) = 0;
    elseif strcmp(init_nls{imod},'l') %linear
        tmod.nly = tmod.nlx;
    elseif strcmp(init_nls{imod},'pq') %positive quadratic
        tmod.nly = tmod.nlx.^2;
    elseif strcmp(init_nls{imod},'nq') %negative quadratic
        tmod.nly = -tmod.nlx.^2;
    else
        error('invalid NL type')
    end
    tmod_range = range(tmod.nly);
    tmod.nly = tmod.nly/tmod_range;
    
    mods = [mods,tmod];
end
glm0   = struct('mods',mods,'const',0,'basis',basis,...
    'nltype',nltype,'LL',0,'LP',0,'lambdaW',0,'mname',mname);

