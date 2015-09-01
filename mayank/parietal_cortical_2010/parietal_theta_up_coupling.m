clear all
% close all

load G:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_individual
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = 2016/2;
lcf = 0.05;
hcf = 40;

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];
desynch_times_mp(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_times_lf8 = desynch_times_lf8(parietal);
% desynch_times_mp = desynch_times_mp(parietal);


winsize = 1;

n = length(sess_data);

%%
freq_bins = linspace(2,6,50);
phase_bins = linspace(-pi,pi,30);
for d = 1:n
    to_dir = sess_data(d).directory;
    to_dir(1) = 'G';
    disp(sess_data(d).name)
    cd(to_dir);

    load used_data wcv_minus_spike lf8 lf4 lf3
    load hsmm_state_seq8_seg_lf_indivdesynch_dp
    lf8_state_seq = hsmm_bbstate_seq8;
    load hsmm_state_seq_seg_lf_indivdesynch_dp
    mp_state_seq = hsmm_bbstate_seq;

    [wcv_f,t_axis,Fs] = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    %     [lf8_f,t_axis,Fs] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    [b,a] = butter(2,[0.05/niqf 30/niqf]);
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);

    [b2,a2] = butter(2,[5/niqf 20/niqf]);
    lf8_hf = filtfilt(b2,a2,lf8);
    lf8_hf = downsample(lf8_hf,dsf)/sess_data(d).gains(8);


    winsize = 1;
%         [theta_freqs{d},theta_phases{d}] = ...
%             bump_det_test_ext(lf8_f,lf8_hf,lf8_state_seq,lf8_state_seq,hmm8,hmm,Fsd,winsize,1);
% 
    [theta_freqs{d},ltheta_phases{d},theta_phases{d}] = bump_det_test_ext_f...
        (lf8_f,wcv_f/5,lf8_state_seq,mp_state_seq,hmm8,hmm,Fsd,winsize,0);

    n_used_downs(d) = length(theta_freqs{d});
    theta_mean_phase(d) = circ_mean(theta_phases{d});
    theta_kappa(d) = circ_kappa(theta_phases{d});
    theta_rtest(d) = circ_rtest(theta_phases{d});
    theta_otest(d) = circ_otest(theta_phases{d});
    theta_raotest(d) = circ_raotest(theta_phases{d});
    theta_mean_lphase(d) = circ_mean(ltheta_phases{d});
    theta_lkappa(d) = circ_kappa(ltheta_phases{d});
    theta_lrtest(d) = circ_rtest(ltheta_phases{d});
    theta_lotest(d) = circ_otest(ltheta_phases{d});
    theta_lraotest(d) = circ_raotest(ltheta_phases{d});

    %     temp = mod(theta_phases{d}+pi,2*pi)-pi;
    hist_phase(d,:) = hist(theta_phases{d},phase_bins);
    hist_lphase(d,:) = hist(ltheta_phases{d},phase_bins);
    %     hist_freqs(d,:) = hist(theta_freqs{d},freq_bins);
%     figure
%     bar(phase_bins,hist_phase(d,:))
%     xlim([-pi pi])
%     figure
%     bar(phase_bins,hist_lphase(d,:))
%     xlim([-pi pi])
%         pause
%         close all

end