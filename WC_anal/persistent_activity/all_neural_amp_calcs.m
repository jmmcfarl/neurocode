clear all
close all
load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

Fs = 2016;
niqf = Fs/2;
dsf = 8;
Fsd = Fs/dsf;

lcf = 0.05/niqf;
hcf = 4/niqf;
[blow,alow] = butter(2,[lcf hcf]);
maxlag = 10*Fsd;
lags = -maxlag:maxlag;
% lcf = 2/niqf;
% hcf = 10/niqf;
% [bmid,amid] = butter(2,[lcf hcf]);
%
% lcf = 100/niqf;
% hcf = 250/niqf;
% [bhigh,ahigh] = butter(2,[lcf hcf]);

win = 20;

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf2 lf3 lf5 lf6 lf7 lf8 wcv_minus_spike
    if exist('spike_time.mat') > 0
        load spike_time
    else
        load spike_time_br
    end
    %    lf2_theta = filtfilt(bmid,amid,lf2);
    %    lf2_ripple = filtfilt(bhigh,ahigh,lf2);
    lf2_uds = filtfilt(blow,alow,lf2);

    lf3_uds = filtfilt(blow,alow,lf3);

    lf5_uds = filtfilt(blow,alow,lf5);
%     lf6_uds = filtfilt(blow,alow,lf6);
    lf7_uds = filtfilt(blow,alow,lf7);
    lf8_uds = filtfilt(blow,alow,lf8);
    wcv = filtfilt(blow,alow,wcv_minus_spike);

    %    lf2_theta = downsample(lf2_theta,dsf);
    %    lf2_ripple = downsample(lf2_ripple,dsf);
    lf2_uds = downsample(lf2_uds,dsf);
    lf3_uds = downsample(lf3_uds,dsf);
    lf5_uds = downsample(lf5_uds,dsf);
%     lf6_uds = downsample(lf6_uds,dsf);
    lf7_uds = downsample(lf7_uds,dsf);
    lf8_uds = downsample(lf8_uds,dsf);
    wcv = downsample(wcv,dsf);

    spike_id = round(spkid/dsf);

    %    lf2_theta = zscore(lf2_theta);
    %    lf2_ripple = zscore(lf2_ripple);
    lf2_uds = zscore(lf2_uds);
%     lf6_uds = zscore(lf6_uds);
    lf7_uds = zscore(lf7_uds);
    lf8_uds = zscore(lf8_uds);
    lf5_uds = zscore(lf5_uds);
    lf3_uds = zscore(lf3_uds);
    wcv = zscore(wcv);
[y2{d},x2{d}] = gpkde(lf2_uds,-3);
[y3{d},x3{d}] = gpkde(lf3_uds,-3);
[y5{d},x5{d}] = gpkde(lf5_uds,-3);
[ymp{d},xmp{d}] = gpkde(wcv,-3);
[y8{d},x8{d}] = gpkde(lf8_uds,-3);
close all

plot(xmp{d},ymp{d})
hold on
plot(x2{d},y2{d},'c')
plot(x3{d},y3{d},'k')
plot(x5{d},y5{d},'g')
plot(x8{d},y8{d},'r')
xlim([-4 4])
t_names = ['C:\WC_Germany\Persistent_activity\all_dist\narrow_filt_' f_names{d}];
print('-dpng',t_names);
close
end