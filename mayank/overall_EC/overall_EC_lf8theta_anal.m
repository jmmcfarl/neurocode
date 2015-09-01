clear all
% close all

load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = 2016/2;
lcf = 0.05;
hcf = 40;

winsize = 1;

n = length(sess_data);

d = 1;

%%
freq_bins = linspace(2,6,50);
phase_bins = linspace(-pi,pi,30);
% for d = 1:n
    to_dir = sess_data(d).directory;
    to_dir(1) = 'G';
    disp(sess_data(d).name)
    cd(to_dir);

    load used_data wcv_minus_spike lf8 
%     lf8 = lf5;
    load ec_hmm_state_seq8
    lf8_state_seq = hmm_bbstate_seq8;
    load ec_hmm_state_seq
    mp_state_seq = hmm_bbstate_seq;

    [wcv_f,t_axis,Fs] = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    %     [lf8_f,t_axis,Fs] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    [b,a] = butter(2,[0.05/niqf 30/niqf]);
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);

%     [b2,a2] = butter(2,[5/niqf 20/niqf]);
%     lf8_hf = filtfilt(b2,a2,lf8);
%     lf8_hf = downsample(lf8_hf,dsf)/sess_data(d).gains(8);


    winsize = 1;
%         [theta_freqs{d},theta_phases{d}] = ...
%             bump_det_test_ext(lf8_f,lf8_hf,lf8_state_seq,lf8_state_seq,hmm8,hmm,Fsd,winsize,1);
% 
    [ltheta_phases{d}] = ov_EC_bump_det_test(lf8_f,lf8_state_seq,hmm8,Fsd,winsize,1);
    if ~isempty(ltheta_phases{d})
        n_used_downs(d) = length(ltheta_phases{d});
        theta_mean_lphase(d) = circ_mean(ltheta_phases{d});
        theta_lkappa(d) = circ_kappa(ltheta_phases{d});
        theta_lrtest(d) = circ_rtest(ltheta_phases{d});
        theta_lotest(d) = circ_otest(ltheta_phases{d});
        theta_lraotest(d) = circ_raotest(ltheta_phases{d});

        hist_lphase(d,:) = hist(ltheta_phases{d},phase_bins);
        figure
        bar(phase_bins,hist_lphase(d,:))
        xlim([-pi pi])
        pause
        close all
    end
% end