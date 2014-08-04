function STCcf_mat = get_STCcf_mat(glmod)

nmods = length(glmod.mods);
STCcf_mat = zeros(glmod.STCdim,nmods);
for i = 1:nmods
   STCcf_mat(:,i) = glmod.mods(i).STCcf; 
end