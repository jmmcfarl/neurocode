clear all
load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
% load C:\WC_Germany\Persistent_activity\sigmoid_fit\sig_fit_data
load C:\WC_Germany\Persistent_activity\lf8_period_f_data2
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

d = 1;
dsf = 4;
Fsd = 2016/dsf;

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
sm_length = 1+round(0.4*Fsd);
sm_win = ones(1,sm_length)/sm_length;

maxlag = 2*Fsd;
lags = -maxlag:maxlag;

[b,a] = butter(2,[lcf hcf]);

lcf = 40/niqf;
hcf = 100/niqf;
[b2,a2] = butter(2,[lcf hcf]);

% lcf = 60/niqf;
% hcf = 200/niqf;
% [b3,a3] = butter(2,[lcf hcf]);
% 

cd(dir_array{d})

pwd

load used_data wcv_minus_spike lf8
wcv = wcv_minus_spike;
wcv_f_low = filtfilt(b,a,wcv);
lf8_f_low = filtfilt(b,a,lf8);

wcv_low = downsample(wcv_f_low,dsf);
lf8_low = downsample(lf8_f_low,dsf);

wcv_low = zscore(wcv_low);
lf8_low = zscore(lf8_low);

wcv_f_med = filtfilt(b2,a2,wcv);
wcv_med = downsample(wcv_f_med,dsf);
wcv_med = sqrt(wcv_med.^2);
wcv_med = conv(wcv_med,sm_win);
wcv_med(1:floor(sm_length/2)) = [];
wcv_med(end-floor(sm_length/2)+1:end) = [];
wcv_med = zscore(wcv_med);

% wcv_f_high = filtfilt(b3,a3,wcv);
% wcv_high = downsample(wcv_f_high,dsf);
% wcv_high = sqrt(wcv_high.^2);
% wcv_high = conv(wcv_high,sm_win);
% wcv_high(1:floor(sm_length/2)) = [];
% wcv_high(end-floor(sm_length/2)+1:end) = [];
% wcv_high = zscore(wcv_high);

t = (1:length(wcv_med))/Fsd;


up_trans{d} = up_trans{d}*2;
down_trans{d} = down_trans{d}*2;
up_trans8{d} = up_trans8{d}*2;
down_trans8{d} = down_trans8{d}*2;


mp_up_times = zeros(size(t));
for i = 1:length(up_trans{d})
    mp_up_times(up_trans{d}(i):down_trans{d}(i)) = 1;
end

mp_up_med = wcv_med;
mp_up_med(~logical(mp_up_times)) = nan;

mp_down_med = wcv_med;
mp_down_med(logical(mp_up_times)) = nan;


%now lfp up and down trig averages
mp_up_lf8up_avg = zeros(length(up_trans8{d}),length(lags));
mp_up_lf8down_avg = mp_up_lf8up_avg;
mp_down_lf8up_avg = mp_up_lf8up_avg;
mp_down_lf8down_avg = mp_up_lf8up_avg;

for i = 1:length(up_trans8{d})
   
    mp_up_lf8up_avg(i,:) = mp_up_med(up_trans8{d}(i)-maxlag:up_trans8{d}(i)+maxlag);
    mp_down_lf8up_avg(i,:) = mp_down_med(up_trans8{d}(i)-maxlag:up_trans8{d}(i)+maxlag);
    
end

av_mp_up_lf8up_avg = nanmean(mp_up_lf8up_avg);
av_mp_down_lf8up_avg = nanmean(mp_down_lf8up_avg);






