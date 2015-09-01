clear all
close all

addpath('G:\Code\smoothing\software')
addpath('G:\Code\FullBNT-1.0.4\KPMstats\')
addpath('G:\Code\FullBNT-1.0.4\netlab3.3')
addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8
load specgram_data

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_start_times(interneurons) = [];
desynch_stop_times(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

for d = 1:length(sess_data)

    cd(sess_data(d).directory)
    load hsmm_state_seq_lf_pert15

    desynch_times = [desynch_start_times{d}; desynch_stop_times{d}];
    time = (1:length(hsmm_bbstate_seq))/Fs_bb;
    is_desynch = zeros(size(time));
    for i = 1:size(desynch_times,2)
        is_desynch(time > desynch_times(1,i) & time < desynch_times(2,i)) = 1;
    end
    is_desynch = logical(is_desynch);
    time2 = (1:length(hsmm_state_seq))/Fs;
    is_desynch2 = zeros(size(time2));
    for i = 1:size(desynch_times,2)
        is_desynch2(time2 > desynch_times(1,i) & time2 < desynch_times(2,i)) = 1;
    end
    is_desynch2 = logical(is_desynch2);
    time2(1) = time(1);
    time2(end) = time(end);


    load hsmm_state_seq8_lf_pert15
    [hsmm_bbstate_seq8_shift] = shift_state_trans(hsmm_bbstate_seq8,round(0.1*Fs_bb),round(0.04*Fs_bb));
    hsmm_state_seq8 = round(interp1(time2,hsmm_state_seq8,time))';
    hsmm_state_seq = round(interp1(time2,hsmm_state_seq,time))';
    state_seq1 = hsmm_bbstate_seq;
    state_seq2 = hsmm_bbstate_seq8;
    [up_lags,down_lags,up_times1,up_times2,down_times1,down_times2] = ...
        compute_synch_translags(state_seq1,state_seq2,desynch_times,Fs_bb);
    if length(up_times1) > length(down_times1)
        up_times1(end) = [];
        up_lags(end) = [];
    end
    if length(down_times1) > length(up_times1)
        down_times1(end) = [];
        down_lags(end) = [];
    end

    [t_samples] = find_nearest_tsample(t{d},up_times1);
    tdr_8_samps = tdr_8{d}(t_samples);
    uds_pow8_samps = uds_pow_8{d}(t_samples);
    uds_comfreq8_samps = uds_comfreq_8{d}(t_samples);
    
%     plot(tdr_8_samps,up_lags,'.'), ylim([-0.6 0.3]), pause, clf
%     plot(uds_pow8_samps,up_lags,'.'), ylim([-0.6 0.3]), pause, clf
%     plot(uds_comfreq8_samps,up_lags,'.'), ylim([-0.6 0.3]), pause, clf
%     
%     plot(tdr_8_samps,down_lags,'r.'), ylim([-0.6 0.6]), pause, clf
%     plot(uds_pow8_samps,down_lags,'r.'), ylim([-0.6 0.6]), pause, clf
%     plot(uds_comfreq8_samps,down_lags,'r.'), ylim([-0.6 0.6]), pause, clf
  
    [rew,raw]=mcdcov(up_lags,'plots',0);
    mean_up_lag(d) = rew.center;
    var_up_lag(d) = rew.cov;
    [rew,raw]=mcdcov(down_lags,'plots',0);
    mean_down_lag(d) = rew.center;
    var_down_lag(d) = rew.cov;

    [rew,raw]=mcdcov(tdr_8_samps,'plots',0);
    med_tdr(d) = rew.center;
    [rew,raw]=mcdcov(uds_pow8_samps,'plots',0);
    med_pow(d) =  rew.center;
    [rew,raw]=mcdcov(uds_comfreq8_samps,'plots',0);
    med_comfreq(d) =  rew.center;

    d
end

%%
figure(1)
% plot(mean_up_lag,med_tdr,'.'), hold on
hold on, plot(mean_up_lag(parietal),med_tdr(parietal),'r.')

figure(2)
% plot(mean_up_lag,med_pow,'.')
hold on, plot(mean_up_lag(parietal),med_pow(parietal),'r.')

figure(3)
% plot(mean_up_lag,med_comfreq,'.')
hold on, plot(mean_up_lag(parietal),med_comfreq(parietal),'r.')

figure(4)
% plot(mean_down_lag,med_tdr,'.')
hold on, plot(mean_down_lag(parietal),med_tdr(parietal),'r.')

figure(5)
% plot(mean_down_lag,med_pow,'.')
hold on, plot(mean_down_lag(parietal),med_pow(parietal),'r.')

figure(6)
% plot(mean_down_lag,med_comfreq,'.')
hold on, plot(mean_down_lag(parietal),med_comfreq(parietal),'r.')

figure(7)
% plot(var_up_lag,med_tdr,'.')
hold on, plot(var_up_lag(parietal),med_tdr(parietal),'r.')

figure(8)
% plot(var_up_lag,med_pow,'.')
hold on, plot(var_up_lag(parietal),med_pow(parietal),'r.')

figure(9)
% plot(var_up_lag,med_comfreq,'.')
hold on, plot(var_up_lag(parietal),med_comfreq(parietal),'r.')

figure(10)
% plot(var_down_lag,med_tdr,'.')
hold on, plot(var_down_lag(parietal),med_tdr(parietal),'r.')

figure(11)
% plot(var_down_lag,med_pow,'.')
hold on, plot(var_down_lag(parietal),med_pow(parietal),'r.')

figure(12)
% plot(var_down_lag,med_comfreq,'.')
hold on, plot(var_down_lag(parietal),med_comfreq(parietal),'r.')


