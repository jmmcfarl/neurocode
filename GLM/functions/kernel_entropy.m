function ent = kernel_entropy(STCbvs,STCcfs,sdim)

[kern_l,n_bvs] = size(STCbvs);
kern_t = kern_l/sdim;
spatial_dist = var(reshape(STCbvs*STCcfs,kern_t,sdim));
spatial_dist = spatial_dist/sum(spatial_dist);
ent = -sum(spatial_dist.*log2(spatial_dist));