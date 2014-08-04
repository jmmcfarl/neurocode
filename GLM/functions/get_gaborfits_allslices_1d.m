function [gabor_fits,gabor_fitvals] = get_gaborfits_allslices_1d(k,flen,sdim)

fsdim = length(k)/flen;

hyperparams = [];
% hyperparams.sigma_shape = 3; 
% hyperparams.sigma_scale = 1.25;
% hyperparams.noise_var = 0.01;
% x = linspace(0.01,25,500);
% figure
% plot(x,gampdf(x,hyperparams.sigma_shape,hyperparams.sigma_scale));

k_mat = reshape(k,flen,fsdim);

temp_pow = var(k_mat,[],2);
[~,best_slice_id] = max(temp_pow);
[init_gabor_fit,init_gabor_fitvals] = james_1d_gaborfit(k_mat(best_slice_id,:)',sdim,hyperparams);

for i = 1:flen
    fprintf('Lag %d\n',i);
   [gabor_fits(i),gabor_fitvals(i,:)] = james_1d_gaborfit(k_mat(i,:)',sdim,hyperparams,init_gabor_fit);   
end

