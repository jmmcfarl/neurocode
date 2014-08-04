function [glm0] = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,mname)
%%USAGE: [glm0] = createGLM0(ks,defmod,mname)
%   takes the filters and default module structure to create a glm

[nstc_dim,nmods] = size(STCcf_0);
mods = [];
for imod = 1:nmods
	tmod  = defmod; %initialize module param structure
    tmod.STCcf = STCcf_0(:,imod); %store filter coefs in STC space
	tmod.k = STCbvs*STCcf_0(:,imod); %store current filter coefs
    if mod_signs(imod) == 1
%        tmod.hcon = 1;
%        tmod.hmon = 1;
%        tmod.nlcon = 1;
        tmod.type = 'pos';
    end
    if mod_signs(imod) == -1
%        tmod.hcon = -1;
%        tmod.hmon = -1;
%        tmod.nlcon = 1;
       tmod.h = -tmod.h;
       tmod.type = 'neg';
    end    
    
    %if there is no psc term absorb sign of h into w
    if length(tmod.h) == 1
        tmod.w = tmod.h;
        tmod.h = 1;
    end
	mods = [mods,tmod]; 
end

glm0   = struct('mods',mods,'const',0,'STCdim',nstc_dim,'STCbasis',STCbvs,'LL',0,'LP',0,'lambdaW',0,'mname',mname,...
    'mod_signs',mod_signs,'dim_signs',dim_signs);

glm0 = get_model_COM_STBF(glm0);
% glm0 = get_STCB_prior_covs(glm0);