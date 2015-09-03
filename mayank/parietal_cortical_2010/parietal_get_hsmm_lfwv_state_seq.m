function [hsmm_state_seq,hmm_state_seq,hmm,Fs] = parietal_get_hsmm_lfwv_state_seq(raw_data,raw_Fs,f_names)

addpath('G:\Code\smoothing\software')
addpath('G:\Code\FullBNT-1.0.4\KPMstats\')
addpath('G:\Code\FullBNT-1.0.4\netlab3.3')
addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\hsmm_state_detection')

Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
lf_lcf = 0.05; %low cut-off for the low-freq filter
lf_hcf = 4; %high cut-off for the low-freq filter
lf_hcf_d = 1; %high cut-off for the low-freq dist estimator
hf_lcf = 5;
hf_hcf = 100;
band1 = [10 20];
band2 = [20 40];
band3 = [40 80];

windowSize = 60; %window duration to estimate time varying state means (s)
windowSlide = 5; %window spacing to estimate time varying state means (s)
hf_smooth = 0.06; % default .05 std dev of gaussian smoothing for hf RMS
num_states = 2; %number of hidden states (must be 2 in current implementation)
state_mean_hcf = 0.01; %hcf for estimating time varying state means (0 for constant mean)

min_state_dur = 1;
max_state_dur = round(Fs*10);
dur_range = (1:max_state_dur)/Fs;
dur_hist_range = linspace(0,10,100);

%% initializations
%compute filter coefficients
[b_lf,a_lf] = butter(2,[lf_lcf/niqf lf_hcf/niqf]);
[b_lfd,a_lfd] = butter(2,[lf_lcf/niqf lf_hcf_d/niqf]);
[b_hf,a_hf] = butter(2,[hf_lcf/niqf hf_hcf/niqf]);

min_freq = 10; max_freq = 80; delta_j = 0.075;
k0 = 6; %wavenumber for morlet wavelet
fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));
min_scale = 1/max_freq/fourier_factor;
max_scale = 1/min_freq/fourier_factor;
n_scales = round(log2(max_scale/min_scale)/delta_j);
obs_hf = filtfilt(b_hf,a_hf,raw_data);
dsf1 = 8;
Fsd = raw_Fs/dsf1;
dsf2 = dsf/dsf1;
obs_hf = downsample(obs_hf,dsf1);
[wav_trans,periods,scales,COI] = wavelet(obs_hf,1/Fsd,0,delta_j,min_scale,n_scales);
wfreqs = 1./periods;
inv_scales = 1./scales';
raw_scalogram = abs(wav_trans).^2;
n = size(raw_scalogram,2);
sm_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*raw_scalogram,1,delta_j,scales);

hf_b1 = mean(log(sm_scalogram(find(wfreqs > band1(1) & wfreqs < band1(2)),:)));
hf_b2 = mean(log(sm_scalogram(find(wfreqs > band2(1) & wfreqs < band2(2)),:)));
hf_b3 = mean(log(sm_scalogram(find(wfreqs > band3(1) & wfreqs < band3(2)),:)));
hf_b1(hf_b1 < 0) = 0;
hf_b2(hf_b2 < 0) = 0;
hf_b3(hf_b3 < 0) = 0;

hf_b1 = zscore(downsample(hf_b1,dsf2));
hf_b2 = zscore(downsample(hf_b2,dsf2));
hf_b3 = zscore(downsample(hf_b3,dsf2));
hf_b1 = hf_b1(:);
hf_b2 = hf_b2(:);
hf_b3 = hf_b3(:);

%compute low-freq observations
obs_lf = filtfilt(b_lf,a_lf,raw_data);
obs_lf = zscore(downsample(obs_lf,dsf));
obs_lf = obs_lf(:);

obs_lfd = filtfilt(b_lfd,a_lfd,raw_data);
obs_lfd = zscore(downsample(obs_lfd,dsf));
obs_lfd = obs_lfd(:);

T = length(obs_lf);
t_axis = (1:T)/Fs;

%% initialize hsmm
emiss = [obs_lfd hf_b1 hf_b2 hf_b3];
[hmm] = initialize_hsmm(emiss,Fs,num_states,windowSize,windowSlide,state_mean_hcf);
emiss(:,1) = obs_lf;

[hmm] = jmm_hmmtrain_hsmm(hmm,emiss,10);

[hmm_state_seq,lik_best,delta] = hsms_viterbi(hmm,emiss);

[smoothed_seq,orig_dur,reject_dur,sojourn_times] = thresh_state_smooth(hmm_state_seq,Fs,300,300);
[state_durations] = compute_state_durations(hmm_state_seq,Fs);
n_hmm(1,:) = hist(state_durations{1},dur_hist_range);
n_hmm(2,:) = hist(state_durations{2},dur_hist_range);
n_hmm(1,:) = n_hmm(1,:)/sum(n_hmm(1,:));
n_hmm(2,:) = n_hmm(2,:)/sum(n_hmm(2,:));

[state_durations] = compute_state_durations(smoothed_seq,Fs);
[alpha(1),beta(1)] = gamma_mlfit(state_durations{1});
[alpha(2),beta(2)] = gamma_mlfit(state_durations{2});

%%
hmm.state(1).alpha = alpha(1);
hmm.state(1).beta = beta(1);
hmm.state(2).alpha = alpha(2);
hmm.state(2).beta = beta(2);

hmm.dur_range = dur_range;
hmm.max_state_dur = max_state_dur;
hmm.min_state_dur = min_state_dur;

hmm.P(1,:) = gamma_pmf(dur_range,hmm.state(1).alpha,hmm.state(1).beta);
hmm.P(2,:) = gamma_pmf(dur_range,hmm.state(2).alpha,hmm.state(2).beta);

hmm.P(:,1:min_state_dur) = zeros(num_states,min_state_dur);

numIt = 1;
[hmm2]=hsmm_train_real(hmm,emiss,numIt);

[hsmm_state_seq] = hsms_viterbi_real(hmm2,emiss);
[state_durations] = compute_state_durations(hsmm_state_seq,Fs);
state_durations{1}(state_durations{1}<min_state_dur/Fs) = [];
state_durations{2}(state_durations{2}<min_state_dur/Fs) = [];

hmm = hmm2;

% figure
% subplot(2,1,1)
% n = hist(state_durations{1},dur_hist_range);
% n=  n/sum(n);
% stairs(dur_hist_range+0.05,n), hold on
% stairs(dur_hist_range,n_hmm(1,:),'r')
% plot(dur_hist_range,gamma_pmf(dur_hist_range,hmm2.state(1).alpha,hmm2.state(1).beta),'k')
% xlim([0 10])
% subplot(2,1,2)
% n = hist(state_durations{2},dur_hist_range);
% n=  n/sum(n);
% stairs(dur_hist_range+0.05,n), hold on
% stairs(dur_hist_range,n_hmm(2,:),'r')
% plot(dur_hist_range,gamma_pmf(dur_hist_range,hmm2.state(2).alpha,hmm2.state(2).beta),'k')
% xlim([0 10])
% t_names = ['G:\WC_Germany\parietal_cortical_2010\test_state_durs\hf_' f_names];
% print('-dpng',t_names)
% subplot(2,1,1)
% set(gca,'xscale','log')
% subplot(2,1,2)
% set(gca,'xscale','log')
% t_names = ['G:\WC_Germany\parietal_cortical_2010\test_state_durs\hf_log_' f_names];
% print('-dpng',t_names), close
%
%
% [obs_dist,obsrange,win_t] = compute_slide_dist(obs_lf,windowSize,windowSlide,Fs);
% figure
% pcolor(win_t,obsrange,log(obs_dist'));shading flat
% caxis([-4 0])
% hold on
% plot(t_axis,hmm2.state(1).lf_meanfun,'w')
% plot(t_axis,hmm2.state(2).lf_meanfun,'w')
% t_names = ['G:\WC_Germany\parietal_cortical_2010\test_meanfuns\' f_names];
% print('-dpng',t_names), close

% [obs_dist,obsrange,win_t] = compute_slide_dist(obs_hf,windowSize,windowSlide,Fs);
% figure
% pcolor(win_t,obsrange,log(obs_dist'));shading flat
% caxis([-4 0])
% t_names = ['G:\WC_Germany\parietal_cortical_2010\test_meanfuns\hf_' f_names];
% print('-dpng',t_names), close

