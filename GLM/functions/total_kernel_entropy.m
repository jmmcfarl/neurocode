function tot_ent = total_kernel_entropy(model)

STCbvs = model.STCbasis;
[kern_l,n_bvs] = size(STCbvs);
sdim = model.mods(1).SDIM;
kern_t = kern_l/sdim;
nmods = length(model.mods);
ent = zeros(nmods,1);
for i = 1:nmods
    spatial_dist = var(reshape(STCbvs*model.mods(i).STCcf,kern_t,sdim));
    spatial_dist = spatial_dist/sum(spatial_dist);
    ent(i) = -sum(spatial_dist.*log2(spatial_dist));
end
tot_ent = sum(ent);