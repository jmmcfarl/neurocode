function kl_div = k_l_divergence(p,q)

kl_div = nansum(p.*log(p./q));