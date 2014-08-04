clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data_new_method
load C:\WC_Germany\Persistent_activity\lf8_period_f_data_new_method
load C:\WC_Germany\Persistent_activity\UDS_synch_state_dur\UDS_synch_state_dur_data_new_method

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;
Fsd = 2016/dsf;
backlag = 5*Fsd;
forwardlag = 10*Fsd;
lags = -backlag:forwardlag;
d=14
    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8_f = filtfilt(b,a,lf8);
    wcv_f = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8_f,dsf);
    wcv_d = downsample(wcv_f,dsf);

    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);

    wcv_up_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    lf8_up_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    
    lfp_period_dur{d} = zeros(1,length(synch_ups{d}));
    cor_lfp_period_dur{d} = zeros(1,length(synch_ups{d}));
    mean_lfp_up = zeros(1,length(synch_ups{d}));
    bad_rows = [];
    
    lfp_up_lag{d} = zeros(size(synch_ups{d}));
    
    for i = 1:length(synch_ups{d})

        %find duration of MP up state in units of LFP periods
        lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i)))-lf8_period_f{d}(up_trans{d}(synch_ups{d}(i)));

        %find first lfp down state after MP up transition
        first_in_lfp_down = find(down_trans8{d} > up_trans{d}(synch_ups{d}(i)),1,'first');


        %acquire MP up related lfp transitions
        [dummy,near_lfp_up] = min(abs(up_trans{d}(synch_ups{d}(i)) - up_trans8{d}(synch_ups8{d})));
        lfp_up_lag{d}(i) = up_trans{d}(synch_ups{d}(i))-up_trans8{d}(synch_ups8{d}(near_lfp_up));
        
        prev_lfp_up = find(up_trans8{d}<up_trans{d}(synch_ups{d}(i)),1,'last');
        prev_lfp_down = find(down_trans8{d}<down_trans{d}(synch_ups{d}(i)),1,'last');
        next_lfp_up = find(up_trans8{d}>down_trans{d}(synch_ups{d}(i)),1,'first');
        if ~isempty(prev_lfp_up)
            up_lag = up_trans{d}(synch_ups{d}(i))-up_trans8{d}(prev_lfp_up);
            cor_lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i))-up_lag)-lf8_period_f{d}(up_trans{d}(synch_ups{d}(i))-up_lag);
            lfp_pup(i) = up_trans8{d}(prev_lfp_up);
            lfp_pupdur(i) = up_state_dur8{d}(prev_lfp_up);
        end
        

        %calculate mup triggered averages
        if up_trans{d}(synch_ups{d}(i)) > backlag & up_trans{d}(synch_ups{d}(i)) ...
                < length(wcv_d)-forwardlag
            wcv_up_ctrig_mat(i,:) = wcv_d(up_trans{d}(synch_ups{d}(i))-backlag: ...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
            lf8_up_ctrig_mat(i,:) = lf8_d(up_trans{d}(synch_ups{d}(i))-backlag: ...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
        else
            wcv_up_ctrig_mat(i,:) = nan;
            lf8_up_ctrig_mat(i,:) = nan;
            bad_rows = [bad_rows i];
        end

    end



save C:\WC_Germany\Persistent_activity\persistent_analysis\sing_examp_new_method ...
    wcv_up_ctrig_mat lf8_up_ctrig_mat lags Fsd bad_rows