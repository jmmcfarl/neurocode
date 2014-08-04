function [glm0] = createGLM0_stcb_connl(STCbvs,STCcf_0,defmod,init_nls,nltypes,mname)
%%USAGE: [glm0] = createGLM0(ks,defmod,mname)
%   takes the filters and default module structure to create a glm

[nstc_dim,nmods] = size(STCcf_0);
mods = [];
for imod = 1:nmods
	tmod  = defmod; %initialize module param structure
    tmod.STCcf = STCcf_0(:,imod); %store filter coefs in STC space
	tmod.k = STCbvs*STCcf_0(:,imod); %store current filter coefs
    if strcmp(init_nls{imod},'threshlin') %threshold linear
        tmod.nly = tmod.nlx;
        tmod.nly(tmod.nlx <= 0) = 0;
    elseif strcmp(init_nls{imod},'lin') %linear
        tmod.nly = tmod.nlx;
    elseif strcmp(init_nls{imod},'pquad') %positive quadratic
        tmod.nly = tmod.nlx.^2;
    elseif strcmp(init_nls{imod},'nquad') %negative quadratic
        tmod.nly = tmod.nlx.^2;
        tmod.w = -tmod.w;
    elseif strcmp(init_nls{imod},'pexp')
        tmod.nly = exp(tmod.nlx);
    elseif strcmp(init_nls{imod},'nexp')
        tmod.nly = exp(tmod.nlx);
        tmod.w = -tmod.w;
    elseif strcmp(init_nls{imod},'plexp')
        tmod.nly = log(1+exp(tmod.nlx*4));
    elseif strcmp(init_nls{imod},'nlexp')
        tmod.nly = log(1+exp(tmod.nlx*4));
        tmod.w = -tmod.w;
    else
        error('invalid NL type')
    end
        tmod.nltype = nltypes{imod};

	mods = [mods,tmod]; 
end

glm0   = struct('mods',mods,'const',0,'STCdim',nstc_dim,'STCbasis',STCbvs,'LL',0,'LP',0,'lambdaW',0,'mname',mname);

glm0 = get_model_COM_STBF(glm0);
% glm0 = get_STCB_prior_covs(glm0);