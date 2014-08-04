function [posterior,posterior_gradient] = smooth_stcb_posterior(k,k_obs,noise_sigma,inv_prior_cov)

likelihood = 1/2*noise_sigma^2*(k_obs-k)'*(k_obs-k);
prior = 1/2*k'*inv_prior_cov*k;

posterior = likelihood + prior;

likelihood_grad = 1/noise_sigma^2*(k_obs - k);
prior_grad = inv_prior_cov*k;
posterior_gradient = likelihood_grad + prior_grad;