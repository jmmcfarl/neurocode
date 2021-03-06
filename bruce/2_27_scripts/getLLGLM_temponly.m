function [nll, pnll] = getLLGLM_temponly(glmod,T,spkbs,disp_type)
%%*** USAGE: [nll, pnll] = getLLGLM(glmod,stim,spkbs)
% get LL of the full model

if nargin < 4
    disp_type = 'none';
end

hlen = 1; %no psc term
kx    = glmod.const; %initialize kx with model constant term
tx = T*glmod.sacmod';
kx = kx + tx;
kx = kx(hlen:end);

for i = 1:length(glmod.mods)
   glmod.mods(i).lambda_L2x=0;
   glmod.mods(i).locLambda=0;
   glmod.mods(i).lambda_dX=0;
   glmod.mods(i).lambda_dT=0;
   glmod.mods(i).lambda_L1x=0;  
end
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