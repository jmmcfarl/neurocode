function [state_means,state_t,obs_dist,obsrange,data_llik,data_mprior,dist_skew] = compute_state_means_hsmm(obs,windowSize,windowSlide,min_lik,hcf,Fs)

total_dur = length(obs)/Fs;
numWins = floor((total_dur-windowSize)/windowSlide);

state_t = (0:numWins-1)*windowSlide+windowSize/2;

mix=gmm(1,2,'diag');
gmm_options(3) = 1e-15; %tolerance
gmm_options(5) = 1;%reset cov if singular values
gmm_options(14) = 200; %max iterations

state_means = zeros(numWins,2);
data_llik = zeros(numWins,1);
data_mprior = zeros(numWins,2);
dist_skew = zeros(numWins,1);
obs_dist = nan(numWins,401);
for w = 1:numWins
    begInd = (w-1)*windowSlide*Fs+1;
    endInd = min((begInd+windowSize*Fs),length(obs));
    obs_seg = obs(begInd:endInd);
    mix=gmm(1,2,'diag');
    [mix,options] = gmmem(mix, obs_seg, gmm_options);
    data_llik(w) = options(8);
    [state_means(w,:),state_order] = sort(mix.centres);
    data_mprior(w,:) = mix.priors(state_order);
    dist_skew(w) = skewness(obs_seg);
    [obs_dist(w,:),obsrange] = gpkde(obs_seg,-3,[min(obs);max(obs)]);
end

% low_lik = find(data_llik > min_lik);
% good_lik = find(data_llik < min_lik);
low_ent = find(max(data_mprior,[],2) > 0.9);
good_ent = find(max(data_mprior,[],2) < 0.9);
avg_uds_amp = mean(state_means(good_ent,2)-state_means(good_ent,1));
low_lik_ups = low_ent(dist_skew(low_ent) < 0);
low_lik_downs = low_ent(dist_skew(low_ent) > 0);
state_means(low_lik_ups,1) = state_means(low_lik_ups,2)-avg_uds_amp;
state_means(low_lik_downs,2) = state_means(low_lik_downs,1)+avg_uds_amp;

state_t = [0 state_t total_dur];
state_means = [state_means(1,:);state_means;state_means(end,:)];