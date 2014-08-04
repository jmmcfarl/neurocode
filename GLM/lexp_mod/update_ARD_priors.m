function mod = update_ARD_priors(mod)

D = length(mod.mods(1).k);
nmods = length(mod.mods);

for i = 1:nmods
   mod.mods(i).lambda_L2x = D/sum(mod.mods(i).k.^2);     
end