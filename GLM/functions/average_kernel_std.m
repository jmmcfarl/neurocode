function [avg_std,kstd] = average_kernel_std(model)

STCbvs = model.STCbasis;
[kern_l,n_bvs] = size(STCbvs);
sdim = model.mods(1).SDIM;
kern_t = kern_l/sdim;
nmods = length(model.mods);
kstd = zeros(nmods,1);
for i = 1:nmods
kstd(i) = kernel_std(model.STCbasis,model.mods(i).STCcf,sdim);
end
avg_std = mean(kstd);