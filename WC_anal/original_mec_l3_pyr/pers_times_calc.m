clear all

load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\correl_data
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\lag_data

dsf = 8;
Fsd = 2016/dsf;
winsize = 30;
niqf = 2016/2;
hcf = 40/niqf;
lcf = 0.05/niqf;
[b,a] = butter(2,[lcf hcf]);
cnt = 0;
maxlag = round(5*Fsd);
tlags = [-maxlag:maxlag]/Fsd;

pers_dur = 5;
back_time = 5;
forward_time = 30;
min_dur = 30;
no_pa = [];
pers_thresh = 5;

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8 = filtfilt(b,a,lf8);
    wcv = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8,dsf);
    wcv_d = downsample(wcv,dsf);

    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);

    total_dur = floor(length(lf8_d)/Fsd);

    %find all wcv up states lasting more than pers_dur seconds
    pers_ups = round(up_trans{d}(find(up_state_dur{d} > pers_thresh))/Fsd);
    pers_ups8 = round(up_trans8{d}(find(up_state_dur8{d} > pers_thresh))/Fsd);

    pers_begs = pers_ups - back_time;
    pers_ends = pers_ups + forward_time;
    pers_begs8= pers_ups8-back_time;
    pers_ends8 = pers_ups8+forward_time;

    pers_secs = [];
    pers_secs8 = [];
    for be = 1:length(pers_begs)
        pers_secs = [pers_secs pers_begs(be):pers_ends(be)];
    end
    for be = 1:length(pers_begs8)
        pers_secs8 = [pers_secs8 pers_begs8(be):pers_ends8(be)];
    end

    pers_secs = unique(pers_secs);
    pers_secs8 = unique(pers_secs8);
    
    pers_secs(pers_secs <= 0) = [];
    pers_secs(pers_secs > total_dur) = [];
    
    pers_secs8(pers_secs8 <= 0) = [];
    pers_secs8(pers_secs8 > total_dur) = [];

    pers_secs = setdiff(pers_secs,pers_secs8);
    npers_secs = setdiff([1:total_dur],pers_secs);
    npers_secs = setdiff(npers_secs,pers_secs8);
    
    pers_measure = zeros(1,total_dur);
    pers_measure(pers_secs) = 1;
    pers_diff = [0 diff(pers_measure)];
    pers_start = find(pers_diff==1);
    pers_stop = find(pers_diff == -1);
    
    npers_measure = zeros(1,total_dur);
    npers_measure(npers_secs) = 1;
    npers_diff = [0 diff(npers_measure)];
    npers_start = find(npers_diff ==1);
    npers_stop = find(npers_diff == -1);
    if npers_measure(1) == 1
        npers_start = [1 npers_start];
    elseif pers_measure(1) == 1
        pers_start = [1 pers_start];
    end
    if npers_measure(end) == 1
        npers_stop = [npers_stop total_dur];
    elseif pers_measure(end) == 1
        pers_stop = [pers_stop total_dur];
    end
    
    pers_dur = pers_stop - pers_start;
    npers_dur = npers_stop - npers_start;
    
    bad_pers = find(pers_dur < min_dur);
    bad_npers = find(npers_dur < min_dur);
    
    pers_start(bad_pers) = [];
    pers_stop(bad_pers) = [];
    npers_start(bad_npers) = [];
    npers_stop(bad_npers) = [];
   
    pers_start_times{d} = pers_start;
    pers_stop_times{d} = pers_stop;
    npers_start_times{d} = npers_start;
    npers_stop_times{d} = npers_stop;
    
    clear pers_secs npers_secs
    pers_secs = [];
    npers_secs = [];
    for p = 1:length(pers_start_times{d})
       pers_secs = [pers_secs pers_start_times{d}(p):pers_stop_times{d}(p)];
    end
    for p = 1:length(npers_start_times{d})
        npers_secs = [npers_secs npers_start_times{d}(p):npers_stop_times{d}(p)];
    end
    percent_pers(d) = length(pers_secs)/(length(pers_secs) + length(npers_secs));
    
end

save C:\WC_Germany\JMM_analysis_pyr\heka_anal\pers_start_times percent_pers pers_start_times pers_stop_times npers_start_times npers_stop_times
