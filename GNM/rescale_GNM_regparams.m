function gnm = rescale_GNM_regparams(gnm,old_scales,new_scales)

nmods = length(gnm.mods);
for n = 1:nmods
   resc_factor = old_scales(n)/new_scales(n);
   
   %adjust L2 norm regularization
   gnm.mods(n).lambda_L2x = gnm.mods(n).lambda_L2x*resc_factor^2;
   gnm.mods(n).lambda_dX = gnm.mods(n).lambda_dX*resc_factor^2;
   gnm.mods(n).lambda_dT = gnm.mods(n).lambda_dT*resc_factor^2;
   gnm.mods(n).lambda_d2T = gnm.mods(n).lambda_d2T*resc_factor^2;
   gnm.mods(n).lambda_d2X = gnm.mods(n).lambda_d2X*resc_factor^2;
   gnm.mods(n).lambda_d2XT = gnm.mods(n).lambda_d2XT*resc_factor^2;
   
   %adjust L1 norm regularization
   gnm.mods(n).lambda_L1x = gnm.mods(n).lambda_L1x*resc_factor;
end