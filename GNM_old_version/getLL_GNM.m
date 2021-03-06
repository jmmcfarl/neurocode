function [nll, pnll, lpen, prate, g, int_g] = getLL_GNM(glmod,X,spkbs,disp_type)

if nargin < 4
    disp_type = 'none';
end

k_mat = get_k_mat(glmod);
g_mat = X*k_mat;

fsdim = glmod.stim_params.fsdim; %number of pixels
flen = glmod.stim_params.flen;

% g    = glmod.const; %initialize kx with model constant term
g = glmod.spk_theta; %initialize kx with model constant term
nmods = length(glmod.mods);
min_x = min(glmod.mods(1).nlx);
for imod = 1:nmods %loop across NL modules
    
    mod = glmod.mods(imod);
    
    if strcmp(glmod.mods(imod).nltype,'lexp')
        int_g{imod} = 1/mod.lexp_beta*log(1+exp(mod.lexp_beta*(g_mat(:,imod)-mod.lexp_theta))) - ...
            1/mod.lexp_beta*log(1+exp(mod.lexp_beta*(min_x-mod.lexp_theta)));
    elseif strcmp(glmod.mods(imod).nltype,'quad')
        int_g{imod} = g_mat(:,imod).^2;
    elseif strcmp(glmod.mods(imod).nltype,'lin')
        int_g{imod} = g_mat(:,imod);
    elseif strcmp(glmod.mods(imod).nltype,'threshlin')
        int_g{imod} = g_mat(:,imod);
        int_g{imod}(int_g{imod} < 0) = 0;
    elseif strcmp(glmod.mods(imod).nltype,'rquad')
        int_g{imod} = g_mat(:,imod).^2;
        int_g{imod}(g_mat(:,imod) < 0) = 0;
    elseif strcmp(mod.nltype,'uncon')
        int_g{imod} = nlin_proc_stim(g_mat(:,imod),mod.nly,mod.nlx); %pass internal output through NL
        int_g{imod}(isnan(int_g{imod})) = 0;
    else
        error('Invalid nl type');
    end
    g = g + mod.w*int_g{imod};
end

[nll,pnll,lpen,prate] = getLL_GNM_internal(g,glmod,spkbs);

if strcmp(disp_type,'full')
    disp_level = 3;
elseif strcmp(disp_type,'tots')
    disp_level = 2;
elseif strcmp(disp_type,'min')
    disp_level = 1;
elseif strcmp(disp_type,'none')
    disp_level = 0;
else
    error('Invalid disp type');
end

if disp_level > 0
    fprintf('LL:%.6g |Lpost:%.6g |Tot Pen:%.6g\n',nll,pnll,(pnll-nll));
    if disp_level > 1
        fprintf('Totals : Grad:%.3g | KLapl:%.3g | KL1:%.3g | KL2:%.3g \n',...
            sum(lpen.grad),sum(lpen.lapl),sum(lpen.l1x),sum(lpen.l2x));
        if disp_level > 2
            fprintf('***\n');
            for imod = 1:nmods
                fprintf('Module %d : Grad:%.3g | KLapl:%.3g | KL1:%.3g | KL2:%.3g\n',...
                    imod,lpen.grad(imod),lpen.lapl(imod),lpen.l1x(imod),lpen.l2x(imod));
            end
            fprintf('***\n');
        end
        fprintf('.........\n');
    end
end