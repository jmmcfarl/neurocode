function Xnew = get_Xcell_tInds(Xcell,target_inds)

n_targets = length(Xcell);
for ii = 1:n_targets
   Xnew{ii} = Xcell{ii}(target_inds,:); 
end