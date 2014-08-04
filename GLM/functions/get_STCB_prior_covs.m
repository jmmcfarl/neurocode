function model = get_STCB_prior_covs(model)

sdim = model.mods(1).SDIM;
nmods = length(model.mods);
kern_len = length(model.mods(1).k);
kern_t = kern_len/sdim;

[Xspace,Xtime] = meshgrid(1:sdim,1:kern_t);
Xspace_vec = Xspace(:);
for i = 1:nmods
    cur_COM = model.mods(i).COM;
    dist_vec = Xspace_vec-cur_COM;
    gauss_fun = exp(-1/2*dist_vec.^2/model.mods(i).spatial_sigma - model.mods(i).prior_scale);
    model.mods(i).prior_precisionmat = diag(1./gauss_fun);    
end