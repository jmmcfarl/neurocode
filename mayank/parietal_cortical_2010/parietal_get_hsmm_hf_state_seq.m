function [hsmm_state_seq,hmm_state_seq,hmm,Fs] = parietal_get_hsmm_hf_state_seq(raw_data,raw_Fs,f_names)

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
lf_hcf = 2; %high cut-off for the low-freq filter
lf_hcf_d = 1; %high cut-off for the low-freq dist estimator
hf_lcf = 10; %low cut-off for high freq filter
hf_hcf = 80; %high cut-off for high freq filter
windowSize = 60; %window duration to estimate time varying state means (s)
windowSlide = 5; %window spacing to estimate time varying state means (s)
hf_smooth = 0.05; %default 0.05 std dev of gaussian smoothing for hf RMS
num_states = 2; %number of hidden states (must be 2 in current implementation)
state_mean_hcf = 0.01; %hcf for estimating time varying state means (0 for constant mean)

% min_state_dur = round(Fs*0.4);
min_state_dur = 1;
max_state_dur = round(Fs*20);
dur_range = (1:max_state_dur)/Fs;
dur_hist_range = linspace(0,10,100);

%% initializations
%compute filter coefficients
[b_hf,a_hf] = butter(2,[hf_lcf/niqf hf_hcf/niqf]);
[b_lf,a_lf] = butter(2,[lf_lcf/niqf lf_hcf/niqf]);

obs_hf = filtfilt(b_hf,a_hf,raw_data);
obs_hf = sqrt(jmm_smooth_1d_cor(obs_hf.^2,round(hf_smooth*raw_Fs))); %smoothed RMS power
obs_hf = zscore(obs_hf);
obs_hf = log(obs_hf-min(obs_hf)+1); %log transform the rms power to normalize
% obs_hf = log(obs_hf); %log transform the rms power to normalize
obs_hf = zscore(downsample(obs_hf,dsf));
obs_hf = obs_hf(:);

%compute low-freq observations
obs_lf = filtfilt(b_lf,a_lf,raw_data);
obs_lf = zscore(downsample(obs_lf,dsf));
obs_lf = obs_lf(:);

T = length(obs_hf);
t_axis = (1:T)/Fs;

%% initialize hsmm
emiss = [obs_lf obs_hf];
[hmm] = initialize_hsmm(emiss,Fs,num_states,windowSize,windowSlide,state_mean_hcf);
hmm.state(1).var = hmm.state(1).var(2:end,2:end);
hmm.state(2).var = hmm.state(2).var(2:end,2:end);
hmm.state(1).mean = hmm.state(1).hf_mean;
hmm.state(2).mean = hmm.state(2).hf_mean;
hmm.p = hmm.p-1;
emiss(:,1) = [];

% [hmm] = initialize_hsmm_fixmeans(emiss,Fs,num_states);

[hmm] = jmm_hmmtrain_hsmm_fixmean(hmm,emiss,10);

[hmm_state_seq,lik_best,delta] = hsms_viterbi_fixmean(hmm,emiss);

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
[hmm2,gamma]=hsmm_train_real_fixmean(hmm,emiss,numIt);

[hsmm_state_seq,max_lik] = hsms_viterbi_real_fixmean(hmm2,emiss);
[state_durations] = compute_state_durations(hsmm_state_seq,Fs);
state_durations{1}(state_durations{1}<min_state_dur/Fs) = [];
state_durations{2}(state_durations{2}<min_state_dur/Fs) = [];

hmm = hmm2;
hmm.max_lik = max_lik;
hmm.gamma = gamma;


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

% [obs_dist,obsrange,win_t] = compute_slide_dist(obs_hf,windowSize,windowSlide,Fs);
% figure
% pcolor(win_t,obsrange,log(obs_dist'));shading flat
% line([0 win_t(end)],[hmm.state(1).mean hmm.state(1).mean],'Color','w','linewidth',2)
% line([0 win_t(end)],[hmm.state(2).mean hmm.state(2).mean],'Color','w','linewidth',2)
% caxis([-4 0])
% t_names = ['G:\WC_Germany\parietal_cortical_2010\test_meanfuns\hf_' f_names];
% print('-dpng',t_names), close

