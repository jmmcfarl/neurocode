function [prate,fg] = getGLM_predrate(glmod,X)
%%*** USAGE: [nll, pnll] = getLLGLM(glmod,stim,spkbs)
% get LL of the full model

%extract matrix of STCcfs across modules
k_mat = get_k_mat(glmod);
%internal filter outputs of each module
g_mat = X*k_mat;

fsdim = glmod.mods(1).fsdim; %number of pixels
hlen = 1; %no psc term

kx    = glmod.const; %initialize kx with model constant term
nmods = length(glmod.mods);
min_x = min(glmod.mods(1).nlx);
for imod = 1:nmods %loop across NL modules
    
    mod = glmod.mods(imod);
    
    if strcmp(glmod.mods(imod).nltype,'lexp')
        fg{imod} = 1/mod.beta*log(1+exp(mod.beta*(g_mat(:,imod)-mod.theta))) - 1/mod.beta*log(1+exp(mod.beta*(min_x-mod.theta)));
    elseif strcmp(glmod.mods(imod).nltype,'quad')
        fg{imod} = g_mat(:,imod).^2;
    elseif strcmp(glmod.mods(imod).nltype,'lin')
        fg{imod} = g_mat(:,imod);
    elseif strcmp(glmod.mods(imod).nltype,'threshlin')
        fg{imod} = g_mat(:,imod);
        fg{imod}(fg{imod} < 0) = 0;
    elseif strcmp(glmod.mods(imod).nltype,'rquad')
        fg{imod} = g_mat(:,imod).^2;
        fg{imod}(g_mat(:,imod) < 0) = 0;
    elseif strcmp(mod.nltype,'uncon')
        fg{imod}    = nlin_proc_stim(g_mat(:,imod),mod.nly,mod.nlx); %pass internal output through NL
    else
        error('Invalid nl type');
    end
    kx = kx + mod.w * fg{imod};
end

%% pure ll based on model prediction (e.g. likelihood)
kx(kx > 100)    = 100; %saturate input to spiking NL
if strcmp(glmod.spk_nl,'logexp')
    prate = log(1+exp(kx)); %apply spiking NL tp get predicted rate
elseif strcmp(glmod.spk_nl,'exp')
    prate = exp(kx);
else 
    error('Not accepted spiking NL')
end

