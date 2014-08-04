function new_mod = add_null_filter(f0,rand_type,weight,nl_type,beta,theta)

cur_nmods = length(f0.mods);
if cur_nmods < 1
    error('Need a base model with at least one filter')
end

new_mod = f0;
new_mod.mods = [new_mod.mods new_mod.mods(end)];
new_mod.mods(end).w = weight;
if strcmp(rand_type,'zero')
    new_mod.mods(end).k = zeros(size(new_mod.mods(1).k));
elseif strcmp(rand_type,'rand')
    prev_k_mat = get_k_mat(f0);
    k_norms = sqrt(sum(prev_k_mat.^2,1));
    cur_norm = min(k_norms)/2; %initialize to have half the smallest L2 norm
    new_mod.mods(end).k = randn(size(new_mod.mods(1).k));
    new_mod.mods(end).k = new_mod.mods(end).k/norm(new_mod.mods(end).k)*cur_norm;
else
    error('Unsupported initialization method');
end
if strcmp(f0.basis,'white')
    new_mod.mods(end).pix = new_mod.mods(end).k'*f0.pix_conv_mat;
    new_mod.mods(end).pix = new_mod.mods(end).pix(:);
elseif strcmp(f0.basis,'pix')
    new_mod.mods(end).pix = new_mod.mods(end).k(:);
end

if strcmp(nl_type,'lexp')
    new_mod.mods(end).nly = 1/beta*log(1+exp(beta*(new_mod.mods(end).nlx-theta))) - ...
        1/beta*log(1+exp(beta*(new_mod.mods(end).nlx(1)-theta)));
    new_mod.mods(end).beta = beta;
    new_mod.mods(end).theta = theta;
elseif strcmp(nl_type,'quad')
    new_mod.mods(end).nly = new_mod.mods(end).nlx.^2;
elseif strcmp(nl_type,'lin')
    new_mod.mods(end).nly = new_mod.mods(end).nlx;
elseif strcmp(nl_type,'threshlin')
    new_mod.mods(end).nly = new_mod.mods(end).nlx;
    new_mod.mods(end).nly(new_mod.mods(end).nlx < 0) = 0;
end
