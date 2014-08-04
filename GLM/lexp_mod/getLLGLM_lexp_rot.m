function [nll, pnll] = getLLGLM_lexp_rot(glmod,kern_output,spkbs,disp_type)
%%*** USAGE: [nll, pnll] = getLLGLM(glmod,stim,spkbs)
% get LL of the full model

if nargin < 4
    disp_type = 'none';
end

%extract matrix of STCcfs across modules
STCcf_mat = get_STCcf_mat(glmod);
%internal filter outputs of each module
g_mat = kern_output*STCcf_mat;

sdim = glmod.mods(1).SDIM;
fsdim = glmod.mods(1).fsdim;
hlen = 1;

kx    = glmod.const; %initialize kx with model constant term
nmods = length(glmod.mods);
min_x = min(glmod.mods(1).nlx);
for imod = 1:nmods %loop across NL modules
    
    mod = glmod.mods(imod);
    fg = g_mat(:,imod);
    if strcmp(mod.nltype,'lexp')
    fg = 1/mod.beta*log(1+exp(mod.beta*fg)) - 1/mod.beta*log(1+exp(mod.beta*min_x));
    elseif strcmp(mod.nltype,'quad')
        fg = fg.^2;
    elseif strcmp(mod.nltype,'lin')
        
    elseif strcmp(mod.nltype,'threshlin')
        fg(fg < 0) = 0;
    elseif strcmp(mod.nltype,'rquad')
        fg = fg.^2;
        fg(g_mat(:,imod) < 0) = 0;
    elseif strcmp(mod.nltype,'uncon')
     %pass through current internal NL
        fg = nlin_proc_stim(fg,glmod.mods(imod).nly,glmod.mods(imod).nlx);
    else
        error('Unsupported NL')
    end
    kx = kx + mod.w * fg;
    
end

kx = kx(hlen:end);

%  [nll,pnll] = getLL_jmm(kx,glmod,spkbs);
if strcmp(glmod.image_type,'1d')
    [nll,pnll,lpen] = getLL_lexp_rot1d(kx,glmod,spkbs);
elseif strcmp(glmod.image_type,'2d')
    [nll,pnll,lpen] = getLL_lexp_rot2d(kx,glmod,spkbs);
else
    error('Unrecognized image type')
end

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
        fprintf('Totals   ::: W-sparse: %.3g | Loc:%.3g |\n',...
            sum(lpen.w),sum(lpen.loc));
        if disp_level > 2
            fprintf('***\n');
            for imod = 1:nmods
                fprintf('Module %d ::: W-sparse: %.3g | Loc:%.3g |\n',...
                    imod, lpen.w(imod),lpen.loc(imod));
            end
            fprintf('***\n');
        end
        fprintf('.........\n');
    end
end