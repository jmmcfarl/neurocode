function [state_means,state_t,obs_dist,obsrange,uni_times] = locate_state_means(obs,windowSize,windowSlide,Fs)

%load time data
total_dur = length(obs)/Fs;
numWins = floor((total_dur-windowSize)/windowSlide); %number of windows to use
t_axis = (0:numWins-1)*windowSlide+windowSize/2; %window center times

%more initializations
up_peaks = nan(size(t_axis));
down_peaks = nan(size(t_axis));
uni_times = [];

for w = 1:numWins

    begT = (w-1)*windowSlide;
    endT = begT + windowSize;
    begInd = round(begT*Fs+1);
    endInd = min(round(begInd+windowSize*Fs),length(obs));
    obs_seg = obs(begInd:endInd);    
    
    sm_bandwidth = terrell_osbs(obs_seg); %compute terrell's oversmoothing bandwidth
    obsrange = linspace(min(obs),max(obs),401); %range of grid points to evaluate density
    obs_dist(w,:) = ksdensity(obs_seg,obsrange,'width',sm_bandwidth); %kernel density estimate
    
    %find peaks in the density that are separated by at least 0.5z
    delta_r = mean(diff(obsrange)); 
    [peak_heights,peak_locs] = findpeaks(obs_dist(w,:),'minpeakdistance',round(.5/delta_r));

    if length(peak_locs) > 1    %if bimodal

        %sort peaks in density according to their mass
        [~,peak_order] = sort(peak_heights,'descend');
        peak_locs = peak_locs(peak_order);
        %locate peak corresponding to up state and down state
        up_p = max([peak_locs(1) peak_locs(2)]);
        down_p = min([peak_locs(1) peak_locs(2)]);
        up_peaks(w) = obsrange(up_p);
        down_peaks(w) = obsrange(down_p);

    else  %if the distribution is unimodal

        %find skewness of the distribution
        obs_skew = skewness(obs_seg);
        %if the distribution is positively skewed assume it's the down
        %state
        if obs_skew > 0
            %set threshold as first inflection point after the peak
            down_peaks(w) = obsrange(peak_locs);
            %if the distribution is negatively skewed, assume up state
        else
            up_peaks(w) = obsrange(peak_locs);
        end
        
        uni_times = [uni_times w]; %storing vector of times when the distribution was unimodal
        
    end
end

uds_amp = up_peaks-down_peaks;
avg_uds_amp = nanmean(uds_amp);

%when the distribution was unimodal, extrapolate the other mode as being
%separated by the average up-down state difference
up_peaks(isnan(up_peaks)) = down_peaks(isnan(up_peaks))+avg_uds_amp;
down_peaks(isnan(down_peaks)) = up_peaks(isnan(down_peaks))-avg_uds_amp;

%interpolate signals in case there are any missing samples remaining
up_ind = find(~isnan(up_peaks));
up_peaks_int = spline(t_axis(up_ind),up_peaks(up_ind),t_axis);
down_ind = find(~isnan(down_peaks));
down_peaks_int = spline(t_axis(down_ind),down_peaks(down_ind),t_axis);

state_means(1,:) = down_peaks_int;
state_means(2,:) = up_peaks_int;
state_t = t_axis;

