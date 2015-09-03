clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\persistent_revised')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8
load sigmoid_fit_data

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
[b_high,a_high] = butter(2,[20/niqf 40/niqf]);

delay_range = linspace(-0.5,0.3,40);
delta_r = delay_range(2)-delay_range(1);
min_samps = 5;
for d = 1:n
    cd(sess_data(d).directory)
    pwd
    load hsmm_state_seq8_lf_pert15
    load hsmm_state_seq_lf_pert15
    wcv_state_seq = hsmm_bbstate_seq;
    lf8_state_seq = hsmm_bbstate_seq8;

    load used_data lf8
    lf8_f = filtfilt(b_high,a_high,lf8);
    lf8_f = sqrt(jmm_smooth_1d_cor(lf8_f.^2,round(Fsd*0.05)));
    lf8_f = zscore(downsample(lf8_f,dsf));
    time = (1:length(lf8_f))/Fsd;
    
    up_trans = find(wcv_state_seq(1:end-1) == 1 & wcv_state_seq(2:end) == 2);
    up_trans(up_trans < round(Fsd*0.5) | up_trans > length(lf8_f) - round(0.5*Fsd)) = [];
    up_trans8 = find(lf8_state_seq(1:end-1) == 1 & lf8_state_seq(2:end) == 2);
    up_trans8(up_trans8 < round(Fsd*0.5) | up_trans8 > length(lf8_f) - round(0.5*Fsd)) = [];
    down_trans = find(wcv_state_seq(1:end-1) == 2 & wcv_state_seq(2:end) == 1);
    down_trans(down_trans < round(Fsd*0.5) | down_trans > length(lf8_f) - round(0.5*Fsd)) = [];
    down_trans8 = find(lf8_state_seq(1:end-1) == 2 & lf8_state_seq(2:end) == 1);
    down_trans8(down_trans8 < round(Fsd*0.5) | down_trans8 > length(lf8_f) - round(0.5*Fsd)) = [];
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
   
    up_trans = up_trans(good_fits_wup{d});
    lfp_power = zeros(1,length(up_trans));
    lfp_delay = zeros(1,length(up_trans));
    wup_slope = rltau_wup{d}(good_fits_wup{d});
    for i = 1:length(up_trans)    
        [dummy,nearest_lfp_up] = min(abs(up_trans(i)-up_trans8)); 
        lfp_delay(i) = (up_trans(i)-up_trans8(nearest_lfp_up))/Fsd;  
        lfp_power(i) = mean(lf8_f(up_trans8(nearest_lfp_up)-round(Fsd*0.1):up_trans8(nearest_lfp_up)));
    end
    
    pow_delay_fun = nan(size(delay_range));
    for i = 1:length(delay_range)
        cur_ups = find(lfp_delay > delay_range(i) - delta_r/2 & lfp_delay < delay_range(i) + delta_r/2);
        if length(cur_ups) > min_samps
        pow_delay_fun(i) = mean(lfp_power(cur_ups));
        end
    end
    
%     plot(wup_slope,lfp_power,'.')
%     pause
%     clf
%     plot(wup_slope,lfp_delay,'.')
%     pause
%     clf
    plot(lfp_delay,lfp_power,'.'), hold on
    plot(delay_range,pow_delay_fun,'r','linewidth',2)
    xlim([delay_range(1) delay_range(end)])
    pause
    clf
    
    
end
    
    
