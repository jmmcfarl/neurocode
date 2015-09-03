function [peaks,t_axis,emiss_dist,emissrange,uni_wins] = hsmm_uds_init_state_means(emiss,movingwin,Fs)

% initialize state-conditional means using sliding window density estimation. 
%
% Input: 
%       emiss: column vector of signal feature ('emissions' or 'observation
%       sequence') to be used for inferring UDS
%       movingwin: [windowSize windowSlide] in seconds.  Determines the
%          window length and amount of overlap for performing sliding-window
%          density estimation
%       Fs:  sample frequency of emiss vector

% Output:
%       state_means: 2xW matrix of state-conditional means. W=number of windows used
%          (first row = down state and second row = up state) 
%       state_t: 1xW time-axis (in seconds)
%       emiss_dist: WxR matrix of density estimates in each of the W data
%           windows. R is the number of samples of the density.
%       emissrange: 1xR vector of the values of emiss at which the density
%           is estimated
%       uni_wins: vector which contains the indices of any data windows
%           where the density estimate was determined to be unimodal

windowSize = movingwin(1);
windowSlide = movingwin(2);

%set fixed parameters
min_peak_sep = 0.5; %minimum separation of modes of the density (in z)
n_density_pts = 400; %number of grid points to estimate the density
emissrange = linspace(min(emiss),max(emiss),n_density_pts); %range of grid points to evaluate density
delta_r = emissrange(2)-emissrange(1);

total_dur = length(emiss)/Fs; %total duration in seconds
numWins = floor((total_dur-windowSize)/windowSlide); %number of windows to use
t_axis = (0:numWins-1)*windowSlide+windowSize/2; %window center times

%initializations
up_peaks = nan(size(t_axis)); %store modes associated with the up state
down_peaks = nan(size(t_axis)); %store modes associated with the down state
uni_wins = []; %vector storing the window indices where the emissions distribution is determined to be unimodal

emiss_dist = nan(numWins,n_density_pts); %initialize matrix of density estimates
for w = 1:numWins
    begInd = round((w-1)*windowSlide*Fs + 1);
    endInd = min(length(emiss),round(begInd + windowSize*Fs));
    emiss_seg = emiss(begInd:endInd);    
    
    sm_bandwidth = terrell_osbs(emiss_seg); %compute terrell's oversmoothing bandwidth
    emiss_dist(w,:) = ksdensity(emiss_seg,emissrange,'width',sm_bandwidth); %kernel density estimate
    
    %find peaks in the density that are separated by at least 0.5z
    [peak_heights,peak_locs] = findpeaks(emiss_dist(w,:),'minpeakdistance',round(min_peak_sep/delta_r));

    if length(peak_locs) > 1    %if bimodal
        %sort peaks in density according to their mass
        [~,peak_order] = sort(peak_heights,'descend');
        peak_locs = peak_locs(peak_order);
        %locate peak corresponding to up state and down state
        up_p = max([peak_locs(1) peak_locs(2)]);
        down_p = min([peak_locs(1) peak_locs(2)]);
        up_peaks(w) = emissrange(up_p);
        down_peaks(w) = emissrange(down_p);
    else  %if the distribution is unimodal
        %find skewness of the distribution
        emiss_skew = skewness(emiss_seg);
        %if the distribution is positively skewed assume it's the down
        %state
        if emiss_skew > 0
            %set threshold as first inflection point after the peak
            down_peaks(w) = emissrange(peak_locs);
            %if the distribution is negatively skewed, assume up state
        else
            up_peaks(w) = emissrange(peak_locs);
        end      
        uni_wins = [uni_wins w]; %storing vector of times when the distribution was unimodal        
    end
end

peaks = [down_peaks(:) up_peaks(:)];


