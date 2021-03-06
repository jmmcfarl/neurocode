function [nll, pnll] = getLLGLM_lexp_tempmod(glmod,X,spkbs,T,disp_type)
%%*** USAGE: [nll, pnll] = getLLGLM(glmod,stim,spkbs)
% get LL of the full model

if nargin < 4
    disp_type = 'none';
end

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
        fg = 1/mod.beta*log(1+exp(mod.beta*g_mat(:,imod))) - 1/mod.beta*log(1+exp(mod.beta*min_x));
    elseif strcmp(glmod.mods(imod).nltype,'quad')
        fg = g_mat(:,imod).^2;
    elseif strcmp(glmod.mods(imod).nltype,'lin')
        fg = g_mat(:,imod);
    elseif strcmp(glmod.mods(imod).nltype,'threshlin')
        fg = g_mat(:,imod);
        fg(fg < 0) = 0;
    else
        error('Invalid nl type');
    end
    kx = kx + mod.w * fg;
end

tx = T*glmod.sacmod';
kx = kx + tx;
kx = kx(hlen:end);




[nll,pnll,lpen] = getLL_lexp(kx,glmod,spkbs);

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
        fprintf('Totals   ::: W-sparse: %.3g | Loc:%.3g | Ksmooth:%.3g | KL1:%.3g |\n',...
            sum(lpen.w),sum(lpen.loc),sum(lpen.ksmooth),sum(lpen.l1x));
        if disp_level > 2
            fprintf('***\n');
            for imod = 1:nmods
                fprintf('Module %d ::: W-sparse: %.3g | Loc:%.3g | Ksmooth:%.3g |\n',...
                    imod,lpen.w(imod),lpen.loc(imod),lpen.ksmooth(imod));
            end
            fprintf('***\n');
        end
        fprintf('.........\n');
    end
end