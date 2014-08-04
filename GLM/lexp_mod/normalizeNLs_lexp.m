function [glmod,norm_vals] = normalizeNLs_lexp(glmod,X)

% USAGE: [glmod] = normalizeRFs(glmod,stim)
%  rescales all glmod RFs to filter stim with std=1

nmods = length(glmod.mods);
for n =1:nmods
    %renormalize kernel so the output distribution
    %has support within the domain of the NL
    %     if strcmp(glmod.mods(i).nltype,'uncon')
    if strcmp(glmod.basis,'white')
        cur_k = glmod.mods(n).k;
    elseif strcmp(glmod.basis,'pix')
        cur_k = glmod.mods(n).pix;
    end
    
    g = X*cur_k;
    if strcmp(glmod.mods(n).nltype,'lexp')
        bgint = (g-glmod.mods(n).theta)*glmod.mods(n).beta;
        too_large = bgint > 50;
        fgint = 1/glmod.mods(n).beta*log(1+exp(bgint)) - 1/glmod.mods(n).beta*log(1+exp(glmod.mods(n).beta*min_x));
        fgint(too_large) = 1/glmod.mods(n).beta*bgint(too_large) - 1/glmod.mods(n).beta*min_x;
    elseif strcmp(glmod.mods(n).nltype,'uncon')
        fgint = nlin_proc_stim(g,glmod.mods(n).nly,glmod.mods(n).nlx);
    elseif strcmp(glmod.mods(n).nltype,'quad')
        fgint = g.^2;
    elseif strcmp(glmod.mods(n).nltype,'lin')
        fgint = g;
    elseif strcmp(glmod.mods(n).nltype,'threshlin')
        fgint = g;
        fgint(fgint < 0) = 0;
    elseif strcmp(glmod.mods(n).nltype,'rquad')
        fgint = g.^2;
        fgint(g < 0) = 0;
    else
        error('Invalid internal NL');
    end
    
    norm_vals(n) = std(fgint);
    glmod.mods(n).nly = glmod.mods(n).nly/norm_vals(n);
end

