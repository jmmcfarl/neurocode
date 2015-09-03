function [obs_dist,obsrange,win_t] = compute_slide_dist(obs,windowSize,windowSlide,Fs)

total_dur = length(obs)/Fs;
numWins = floor((total_dur-windowSize)/windowSlide);

win_t = (0:numWins-1)*windowSlide+windowSize/2;

obs_dist = nan(numWins,401);
for w = 1:numWins
    begInd = round((w-1)*windowSlide*Fs)+1;
    endInd = min((begInd+round(windowSize*Fs)),length(obs));
    obs_seg = obs(begInd:endInd);
    sm_bandwidth = terrell_osbs(obs_seg); %compute terrell's oversmoothing bandwidth
    obsrange = linspace(min(obs),max(obs),401); %range of grid points to evaluate density
    obs_dist(w,:) = ksdensity(obs_seg,obsrange,'width',sm_bandwidth); %kernel density estimate
end

