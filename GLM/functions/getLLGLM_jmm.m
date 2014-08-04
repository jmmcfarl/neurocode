function [nll, pnll] = getLLGLM_jmm(glmod,stim,spkbs,disp_type)
%%*** USAGE: [nll, pnll] = getLLGLM(glmod,stim,spkbs)
% get LL of the full model

if nargin < 4
    disp_type = 'none';
end

kx    = glmod.const; %initialize kx with model constant term
nmods = length(glmod.mods);
for imod = 1:nmods %loop across NL modules
    
    %module weighting times the module output is the total contribution of
    %each module
    kx = kx + glmod.mods(imod).w * procModul(glmod.mods(imod),stim);
end

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
    fprintf('LL:%.3g |Lpost:%.3g |Tot Pen:%.3g\n',nll,pnll,(pnll-nll));
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