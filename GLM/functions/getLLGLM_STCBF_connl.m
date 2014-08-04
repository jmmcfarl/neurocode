function [nll, pnll] = getLLGLM_STCBF_connl(glmod,kern_output,spkbs,disp_type)
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
hlen = 1;

kx    = glmod.const; %initialize kx with model constant term
nmods = length(glmod.mods);
for imod = 1:nmods %loop across NL modules
    
    mod = glmod.mods(imod);
    fg = g_mat(:,imod);
    
    if strcmp(mod.nltype,'threshlin')
        fg(fg < 0) = 0;
    elseif strcmp(mod.nltype,'quad')
        fg = fg.^2;
    elseif strcmp(mod.nltype,'lexp')
        fg = log(1+exp(4*fg));
    elseif ~strcmp(mod.nltype,'lin')
        error('Invalid NL type')
    end    
    
    kx = kx + mod.w * fg;
    
end

kx = kx(hlen:end);

%  [nll,pnll] = getLL_jmm(kx,glmod,spkbs);
[nll,pnll,lpen] = getLL_jmm_STCBF(kx,glmod,spkbs);

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
        fprintf('Totals   ::: H-slope:%.3g | H-curve:%.3g | NL-slope:%.3g | NL-curve:%.3g | W-sparse: %.3g | Loc:%.3g |\n',...
            sum(lpen.h),sum(lpen.h2),sum(lpen.nl),sum(lpen.nl2),sum(lpen.w),sum(lpen.loc));
        if disp_level > 2
            fprintf('***\n');
            for imod = 1:nmods
                fprintf('Module %d ::: H-slope:%.3g | H-curve:%.3g | NL-slope:%.3g | NL-curve:%.3g | W-sparse: %.3g | Loc:%.3g |\n',...
                    imod, lpen.h(imod),lpen.h2(imod),lpen.nl(imod),lpen.nl2(imod),lpen.w(imod),lpen.loc(imod));
            end
            fprintf('***\n');
        end
        fprintf('.........\n');
    end
end