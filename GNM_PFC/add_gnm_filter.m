function new_mod = add_gnm_filter(f0,init_kern,weight,nl_type)

cur_nmods = length(f0.mods);
if cur_nmods < 1
    error('Need a base model with at least one filter')
end

new_mod = f0;
new_mod.mods = [new_mod.mods new_mod.mods(end)];
new_mod.mods(end).w = weight;
new_mod.mods(end).k = init_kern(:);

if strcmp(nl_type,'lexp')
    new_mod.mods(end).nly = 1/lexp_beta*log(1+exp(lexp_beta*(new_mod.mods(end).nlx-lexp_theta))) - ...
        1/lexp_beta*log(1+exp(lexp_beta*(new_mod.mods(end).nlx(1)-lexp_theta)));
elseif strcmp(nl_type,'quad')
    new_mod.mods(end).nly = new_mod.mods(end).nlx.^2;
elseif strcmp(nl_type,'rquad')
    new_mod.mods(end).nly = new_mod.mods(end).nlx.^2;
    new_mod.mods(end).nly(new_mod.mods(end).nlx < 0) = 0;
elseif strcmp(nl_type,'lin')
    new_mod.mods(end).nly = new_mod.mods(end).nlx;
elseif strcmp(nl_type,'threshlin')
    new_mod.mods(end).nly = new_mod.mods(end).nlx;
    new_mod.mods(end).nly(new_mod.mods(end).nlx < 0) = 0;
end
new_mod.mods(end).nltype = nl_type;