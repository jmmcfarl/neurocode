function new_x = propose_new_x(old_x,radius,n_bvs,n_filts)

pert_x = randn(n_bvs,n_filts)*radius;
new_x = old_x + pert_x;
% for i = 1:n_filts; new_x(:,i) = new_x(:,i)/norm(new_x(:,i)); end;
