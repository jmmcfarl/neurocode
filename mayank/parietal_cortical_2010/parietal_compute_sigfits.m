clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\persistent_revised')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

sess_data = sess_data(parietal);
desynch_start_times = desynch_start_times(parietal);
desynch_stop_times = desynch_stop_times(parietal);
%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_start_times(interneurons) = [];
desynch_stop_times(interneurons) = [];

n = length(sess_data);

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
[b_wide,a_wide] = butter(2,[0.05/niqf 40/niqf]);
[b_high,a_high] = butter(2,[10/niqf 40/niqf]);

for d = 1:n
    cd(sess_data(d).directory)
    pwd
    load hsmm_state_seq8_lf_pert15
    load hsmm_state_seq_lf_pert15
    wcv_state_seq = hsmm_bbstate_seq;
    lf8_state_seq = hsmm_bbstate_seq8;

    load used_data wcv_minus_spike lf8
    wcv_f = filtfilt(b_wide,a_wide,wcv_minus_spike);
    wcv_f = zscore(downsample(wcv_f,dsf));
    lf8_f = filtfilt(b_wide,a_wide,lf8);
    lf8_f = zscore(downsample(lf8_f,dsf));
    time = (1:length(lf8_f))/Fsd;
    
    up_trans = find(wcv_state_seq(1:end-1) == 1 & wcv_state_seq(2:end) == 2);
    up_trans(up_trans < round(Fsd*0.5) | up_trans > length(wcv_f) - round(0.5*Fsd)) = [];
    up_trans8 = find(lf8_state_seq(1:end-1) == 1 & lf8_state_seq(2:end) == 2);
    up_trans8(up_trans8 < round(Fsd*0.5) | up_trans8 > length(wcv_f) - round(0.5*Fsd)) = [];
    down_trans = find(wcv_state_seq(1:end-1) == 2 & wcv_state_seq(2:end) == 1);
    down_trans(down_trans < round(Fsd*0.5) | down_trans > length(wcv_f) - round(0.5*Fsd)) = [];
    down_trans8 = find(lf8_state_seq(1:end-1) == 2 & lf8_state_seq(2:end) == 1);
    down_trans8(down_trans8 < round(Fsd*0.5) | down_trans8 > length(wcv_f) - round(0.5*Fsd)) = [];
    up_trans(up_trans > down_trans(end)) = [];
    down_trans(down_trans < up_trans(1)) = [];
    up_trans8(up_trans8 > down_trans8(end)) = [];
    down_trans8(down_trans8 < up_trans8(1)) = [];


    desynch_times = [desynch_start_times{d} desynch_stop_times{d}];
    bad_ups = [];
    bad_ups8 = [];
    if ~isempty(desynch_times)
        for i = 1:size(desynch_times,1)
            bad_ups = [bad_ups find(up_trans > desynch_times(i,1)*Fsd & up_trans < desynch_times(i,2)*Fsd)];
            bad_ups8 = [bad_ups8 find(up_trans8 > desynch_times(i,1)*Fsd & up_trans8 < desynch_times(i,2)*Fsd)];
        end
        up_trans(bad_ups) = [];
        up_trans8(bad_ups8) = [];
        down_trans(bad_ups) = [];
        down_trans8(bad_ups8) = [];
    end

    [rlid_wup{d},rltime_wup{d},rlamp_wup{d},rlshift_wup{d},rltau_wup{d},rlerror_wup{d}] = ...
        get_lfp_wcv_sigmoid_fit_ut_12_4(up_trans,down_trans,wcv_f,time,Fsd);
    [rlid_lup{d},rltime_lup{d},rlamp_lup{d},rlshift_lup{d},rltau_lup{d},rlerror_lup{d}] = ...
        get_lfp_wcv_sigmoid_fit_ut_12_4(up_trans8,down_trans8,lf8_f,time,Fsd);
    good_fits_wup{d} = find(rlerror_wup{d} > 0.75);
    good_fits_lup{d} = find(rlerror_lup{d} > 0.75);
    
end

cd G:\WC_Germany\parietal_cortical_2010\
save sigmoid_fit_data rl* good*