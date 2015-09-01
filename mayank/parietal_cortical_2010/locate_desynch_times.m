function [desynch_times,P_mp,P_lfp,f,t] = locate_desynch_times(mp,lfp)

%desynch_times (Nx2) vector containing start and stop times (in s) of each
%of N desynchonized epochs
%desynch_indicator indicator function of desynch epochs (does not include
%the buffer region around the actual desynch periods), sample rate 252Hz

dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.1/niqf 40/niqf]);

params.Fs = Fsd;
params.err = 0;
params.tapers = [4 7];
params.fpass = [0 40];

winlength = 15;
winslide = 2.5;
movingwin = [winlength winslide];

% log_theta_thresh = -8;
% log_so_thresh = -6;
log_so_thresh = -6;
log_hf_thresh = -2;
removal_window = round(winlength/winslide/2)+1; %number of window slides to remove around offending position
% min_between_dur = 10; %minimum duration of data between desynchronized epochs (s)
min_ds_dur = 30; %longest duration which we will check if the potential desynchronized epoch is just a really long down state

%% filter,downsample,zscore
mp = filtfilt(b,a,mp);
mp = downsample(mp,dsf);
mp = zscore(mp);

lfp = filtfilt(b,a,lfp);
lfp = downsample(lfp,dsf);
lfp = zscore(lfp);

% compute spectrogram
[P_mp,t,f] = mtspecgramc(mp,movingwin,params);
[P_lfp,t,f] = mtspecgramc(lfp,movingwin,params);
df = f(2)-f(1);
f_so = find(f > 0.2 & f < 1.5);
% f_theta = find(f > 2 & f < 6);
f_hf = find(f > 4);
%find max log power in theta and SO freq bands
% max_theta_mp = max(10*log10(P_mp(:,f_theta)),[],2);
max_so_mp = max(10*log10(P_mp(:,f_so)),[],2);
max_so_lfp = max(10*log10(P_lfp(:,f_so)),[],2);
% max_theta_lfp = max(10*log10(P_lfp(:,f_theta)),[],2);
net_hf_mp = zscore(trapz(10*log10(P_mp(:,f_hf)),2)*df);
net_hf_lfp = zscore(trapz(10*log10(P_lfp(:,f_hf)),2)*df);

%mark points where the SO power crosses threshold
desynch_indicator = zeros(size(t));
% desynch_indicator(max_theta_mp >= log_theta_thresh | max_so_mp <= log_so_thresh...
%     | max_theta_lfp >= log_theta_thresh | max_so_lfp <= log_so_thresh) = 1;
% desynch_indicator(max_so_lfp <= log_so_thresh | max_so_mp <= log_so_thresh) = 1;
desynch_indicator(max_so_lfp <= log_so_thresh) = 1;
desynch_start_ids = 1+find(desynch_indicator(1:end-1) == 0 & desynch_indicator(2:end) == 1);
desynch_stop_ids = 1+find(desynch_indicator(1:end-1) == 1 & desynch_indicator(2:end) == 0);
low_hf_indicator = find(net_hf_mp <= log_hf_thresh | net_hf_lfp <= log_hf_thresh);


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
short_desynchs = find(desynch_durs < min_ds_dur);
long_ds = [];
for i = 1:length(short_desynchs)
    if any(ismember(desynch_start_ids(short_desynchs(i)):desynch_stop_ids(short_desynchs(i)),low_hf_indicator))
        long_ds = [long_ds short_desynchs(i)];
    end
end
desynch_start_ids(long_ds) = [];
desynch_stop_ids(long_ds) = [];

%now make a window around desynchronized times
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

% %make sure that the UDS periods are long enough to use for analysis
% if length(desynch_start_ids) >= 2
%     in_between_dur = [];
%     for i = 2:length(desynch_start_ids)
%         in_between_dur = [in_between_dur desynch_start_ids(i)-desynch_stop_ids(i-1)];
%     end
%     in_between_dur = in_between_dur/movingwin(2);
%     too_short_between = find(in_between_dur < min_between_dur);
%     desynch_start_ids(too_short_between+1) = [];
%     desynch_stop_ids(too_short_between) = [];
% end
% if ~isempty(desynch_start_ids)
%     remaining_dur = (length(t)-desynch_stop_ids(end))/movingwin(2);
%     if remaining_dur < min_between_dur
%         desynch_stop_ids(end) = length(t);
%     end
% end
       

desynch_start_times = t(desynch_start_ids);
desynch_stop_times = t(desynch_stop_ids);
desynch_times = [desynch_start_times(:) desynch_stop_times(:)];





