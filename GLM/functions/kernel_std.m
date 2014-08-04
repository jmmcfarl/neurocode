function kstd = kernel_std(STCbvs,STCcfs,sdim)

[kern_l,n_bvs] = size(STCbvs);
kern_t = kern_l/sdim;
spatial_dist = var(reshape(STCbvs*STCcfs,kern_t,sdim));
spatial_dist = spatial_dist/sum(spatial_dist);

x_ax = 1:sdim;
com = x_ax*spatial_dist';
kstd = (x_ax-com).^2*spatial_dist';
