function LLs = get_pgabor_tempmod_LLs_lfp(gabor_params,avg_pow,X_gabor_out,scales,offsets,err_sigma,blockids)

[~,NT,K] = size(X_gabor_out);

X_energy_out = squeeze(sqrt(X_gabor_out(1,:,:).^2 + X_gabor_out(2,:,:).^2));
gmat = gabor_params(7)*X_energy_out;
n_blocks = length(unique(blockids));

for i = 1:n_blocks
    cur_set = find(blockids == i);
    gmat(cur_set,:) = scales(i)*gmat(cur_set,:)+offsets(i);
end

cur_resid = bsxfun(@minus,gmat,avg_pow');
LLs = -cur_resid.^2/(2*err_sigma^2) + log(1/2/pi);