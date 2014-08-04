clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);
maxlag = 6*Fsd;
lags = -maxlag:maxlag;

for d = 1:length(dir_array)

    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd

    load used_data wcv_minus_spike lf8

    %bandlimit signals
    down_w = filtfilt(b,a,wcv_minus_spike);
    down_8 = filtfilt(b,a,lf8);

    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);

    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);

    lfp_trig_mp = zeros(length(synch_ups8{d}),length(lags));
    lfp_trig_lfp = zeros(length(synch_ups8{d}),length(lags));
    for i = 1:length(synch_ups8{d})
        if up_trans8{d}(synch_ups8{d}(i)) > maxlag & ...
                length(down_8)-up_trans8{d}(synch_ups8{d}(i))> maxlag
            
        lfp_trig_mp(i,:) = down_w(up_trans8{d}(synch_ups8{d}(i))-maxlag:...
            up_trans8{d}(synch_ups8{d}(i))+maxlag);
        lfp_trig_lfp(i,:) = down_8(up_trans8{d}(synch_ups8{d}(i))-maxlag:...
            up_trans8{d}(synch_ups8{d}(i))+maxlag);

        else
           lfp_trig_mp(i,:) = nan;
           lfp_trig_lfp(i,:) = nan;
        end
        
    end

    mean_lfp_trig_mp(d,:) = nanmean(lfp_trig_mp);
    mean_lfp_trig_lfp(d,:) = nanmean(lfp_trig_lfp);
    
    
    clear down_w down_8 wcv* lf8

end


% save C:\WC_Germany\Persistent_activity\corr_data w_acorr lags l_acorr Fsd x_cor