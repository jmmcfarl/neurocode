function k_mat = get_k_mat(glmod)

nmods = length(glmod.mods);
k_mat = zeros(length(glmod.mods(1).k),nmods);
for i = 1:nmods
   k_mat(:,i) = glmod.mods(i).k; 
end