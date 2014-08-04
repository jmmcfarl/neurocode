function [glmod] = set_k_mat(glmod,k_mat)

nmods = length(glmod.mods);
for i = 1:nmods
    glmod.mods(i).k = k_mat(:,i);
end