function [desynch_times,desynch_ids,P,f,t] = locate_desynch_times_individual(sig)

%desynch_times (Nx2) vector containing start and stop times (in s) of each
%of N desynchonized epochs
%desynch_indicator indicator function of desynch epochs (does not include
%the buffer region around the actual desynch periods), sample rate 252Hz

dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.1/niqf 40/niqf]); %bandpass filter between 0.1 and 40 Hz before computing spectrogram.  This is because we don't want artifacts at very low frequencies dominating our slow-oscillation detection.

%multitaper spectrum parameters
params.Fs = Fsd;
params.err = 0;
params.tapers = [4 7]; %[TW K]
params.fpass = [0 40];

winlength = 15; %window length for spectrogram
winslide = 2.5; %spectrogram window spacing
movingwin = [winlength winslide];

log_so_thresh = -6; %this is a fixed threshold on the maximum relative log power in the slow-oscillation range
log_hf_thresh = -2; %this is a fixed threshold on the normalized integral relative power in the 'high-frequency' band to be identified as a deysnchronized epoch
removal_window = 1.5*round(winlength/winslide)+1; %number of window slides to remove around offending position
max_ds_dur = 30; %longest duration which we will check if the potential desynchronized epoch is just a really long down state

%% filter,downsample,zscore
sig = filtfilt(b,a,sig);
sig = downsample(sig,dsf);
sig = zscore(sig);

% compute spectrogram
[P,t,f] = mtspecgramc(sig,movingwin,params);
df = f(2)-f(1);
f_so = find(f > 0.2 & f < 1.5); %define the slow-oscillation band as 0.2-1.5Hz
f_hf = find(f > 4); %define the 'high-frequency' band as > 4Hz

%find max log power in theta and SO freq bands
max_so = max(10*log10(P(:,f_so)),[],2); %maximum log relative power in the SO range
net_hf = zscore(trapz(log10(P(:,f_hf)),2)*df); %z-score of the integral HF power

%mark points where the SO power crosses below threshold
desynch_indicator = zeros(size(t));
desynch_indicator(max_so <= log_so_thresh) = 1;
desynch_start_ids = 1+find(desynch_indicator(1:end-1) == 0 & desynch_indicator(2:end) == 1);
desynch_stop_ids = 1+find(desynch_indicator(1:end-1) == 1 & desynch_indicator(2:end) == 0);
if desynch_indicator(1) == 1 & desynch_indicator(2) == 0
    desynch_start_ids = [1 desynch_start_ids];
end

%find potential instances of super-long down states 
low_hf_indicator = find(net_hf <= log_hf_thresh);


%make sure you start and stop in desynchronized epochs in correct order
if ~isempty(desynch_start_ids)
    if isempty(desynch_stop_ids)
        desynch_stop_ids = length(t);
    else
        if desynch_start_ids(1) > desynch_stop_ids(1)
            desynch_start_ids = [1 desynch_start_ids];
        end
        if desynch_start_ids(end) > desynch_stop_ids(end)
            desynch_stop_ids = [desynch_stop_ids length(t)];
        end
    end
end

if length(desynch_start_ids) ~= length(desynch_stop_ids)
    disp('error start and stop IDs not equal!')
end

%compute the duration of each putative desynchronized epoch
desynch_durs = (desynch_stop_ids-desynch_start_ids)*winslide;
short_desynchs = find(desynch_durs < max_ds_dur); %these are the putative desynchronized epochs which could be instead long down statas

%find an remove putative desynchronized epochs which are really long down
%states
long_ds = [];
for i = 1:length(short_desynchs)
    if any(ismember(desynch_start_ids(short_desynchs(i)):desynch_stop_ids(short_desynchs(i)),low_hf_indicator))
        long_ds = [long_ds short_desynchs(i)];
    end
end
desynch_start_ids(long_ds) = [];
desynch_stop_ids(long_ds) = [];

%now make a window around desynchronized times for data exclusion
for w = 1:length(desynch_start_ids)
    if desynch_start_ids(w) <= removal_window
        desynch_start_ids(w) = 1;
    else
        desynch_start_ids(w) = desynch_start_ids(w)-removal_window;
    end
    if length(t)-desynch_stop_ids(w) <= removal_window
        desynch_stop_ids(w) = length(t);
    else
        desynch_stop_ids(w) = desynch_stop_ids(w)+removal_window;
    end
end

%now make sure there are no overlapping windows
bad_desynch_start = [];
for w = 2:length(desynch_start_ids)
    if desynch_start_ids(w) < desynch_stop_ids(w-1)
        bad_desynch_start = [bad_desynch_start w];
    end
end
desynch_start_ids(bad_desynch_start) = [];
desynch_stop_ids(bad_desynch_start-1) = [];

desynch_start_times = t(desynch_start_ids);
desynch_stop_times = t(desynch_stop_ids);
desynch_times = [desynch_start_times(:) desynch_stop_times(:)];
desynch_ids = [desynch_start_ids(:) desynch_stop_ids(:)];




