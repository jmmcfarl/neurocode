function LL = get_gabor_energy_dualeye_LL(K,Xmat_l,Xmat_r,Robs,XX,YY)
    
cur_mask1 = get_pgabor_mask_v2(XX,YY,K(1:6),0);
cur_mask2 = get_pgabor_mask_v2(XX,YY,K(1:6),pi/2);

mask1_out = Xmat_l*cur_mask1(:);
mask2_out = Xmat_l*cur_mask2(:);
en1 = K(7)*sqrt(mask1_out.^2+mask2_out.^2);

mask1_out = Xmat_r*cur_mask1(:);
mask2_out = Xmat_r*cur_mask2(:);
en2 = K(8)*sqrt(mask1_out.^2+mask2_out.^2);

g = en1 + en2 + K(9);
too_large = find(g > 100);
r = log(1+exp(g));
r(too_large) = g(too_large);

r(r < 1e-20) = 1e-20;

LL = sum(Robs.*log(r)-r);
Nspks = sum(Robs);
LL = -LL/Nspks;

end

