NT = 10000;
noise_std = 1;
scale = 2;
mean_rates = [0.01 0.1 0.5 1 2 5 10 50 100];
for ii = 1:length(mean_rates)
    mean_rate = mean_rates(ii);
%     true_rate = rand(NT,1)*mean_rate*2;
    true_rate = gamrnd(ones(NT,1)*scale,ones(NT,1)*1)*mean_rate/scale;
    obs_mean_rate(ii) = mean(true_rate);
    R_obs = poissrnd(true_rate);
    Y_obs = normrnd(true_rate,noise_std);
    
    full_LL = nansum(R_obs.*log(R_obs) - R_obs);
    mod_LL = nansum(R_obs.*log(true_rate) - true_rate);
    
    null_rate = ones(size(true_rate))*mean(R_obs);
    null_LL = nansum(R_obs.*log(null_rate) - null_rate);
    
    mod_dev = 2*(full_LL-mod_LL);
    null_dev = 2*(full_LL-null_LL);
    pseudo_R2(ii) = 1-mod_dev/null_dev;
    
    cur_r = corr(Y_obs(:),true_rate(:));
    gauss_R2(ii) = cur_r^2;
end