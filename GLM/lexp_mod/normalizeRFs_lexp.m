function [glmod,norm_vals] = normalizeRFs_lexp(glmod,X,kscale,nlxscale)

% USAGE: [glmod] = normalizeRFs(glmod,stim)
%  rescales all glmod RFs to filter stim with std=1
if nargin < 3
    kscale = 0;
end
if nargin < 4
    nlxscale = 0;
end

nmods = length(glmod.mods);
for i =1:nmods
    %renormalize kernel so the output distribution
    %has support within the domain of the NL
    %     if strcmp(glmod.mods(i).nltype,'uncon')
    if strcmp(glmod.basis,'white')
        cur_k = glmod.mods(i).k;
    elseif strcmp(glmod.basis,'pix')
        cur_k = glmod.mods(i).pix;
    end
    norm_vals(i) = std(X*cur_k);
    if norm_vals(i) ~= 0
        glmod.mods(i).k = glmod.mods(i).k/norm_vals(i);
        glmod.mods(i).pix = glmod.mods(i).pix/norm_vals(i);
        if kscale == 1
            glmod.mods(i).kscale = glmod.mods(i).kscale*norm_vals(i);
        end
        if nlxscale ~= 0
            glmod.mods(i).nlx = glmod.mods(i).nlx/norm_vals(i); 
        end
    end
    %     else
    %         norm_vals = ones(nmods,1);
    %     end
end

