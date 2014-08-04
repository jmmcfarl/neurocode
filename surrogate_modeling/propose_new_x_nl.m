function new_x = propose_new_x_nl(old_x,radius,n_bvs,n_filts)

pert_x = randn(n_bvs*n_filts,1)*radius;
new_x{1} = old_x{1} + pert_x;

pert_nl = randn(n_filts*11,1)*radius;
new_x{2} = old_x{2} + pert_nl;