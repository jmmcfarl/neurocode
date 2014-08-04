function [glmod,norm_vals] = normalizeRFs_STCB(glmod,kern_outputs)

% USAGE: [glmod] = normalizeRFs(glmod,stim)
%  rescales all glmod RFs to filter stim with std=1

nmods = length(glmod.mods); 
for i =1:nmods
    cur_STCbcfs = glmod.mods(i).STCcf;
    norm_vals(i) = std(kern_outputs*cur_STCbcfs);
    glmod.mods(i).STCcf = glmod.mods(i).STCcf/norm_vals(i);
    glmod.mods(i).k = glmod.mods(i).k/norm_vals(i);
    if strcmp(glmod.basis,'white')
        glmod.mods(i).pix = glmod.pix_conv_mat'*glmod.mods(i).k;
    else
        glmod.mods(i).pix = glmod.mods(i).k;
    end
end

